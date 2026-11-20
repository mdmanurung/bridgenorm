# Shared-batch-effect variant, one design per marker over all populations.
#
# Ported from SpectraWeaveR R/corrector_heb_shared.R (engine "heb_shared").
#
# Model on aliquot summaries of reference samples, all populations together:
#   y = donor-by-population (fixed) + study (random, shared over populations)
#       + run/batch (random, shared over populations)
#       + batch-by-population deviation (random, own variance) + error.
# The effect of (batch b, population k) is the BLUP difference to the anchor of
#   u_study(b) + v_b + w_{b,k}.
# A population whose summaries in b were dropped by min_cells still gets u + v of b
# with a zero deviation (w has no data there, so its BLUP is 0).

# Negative REML log-likelihood (up to a constant); Zq is a list of tcrossprod(Z) matrices.
.hebs_reml <- function(log_theta, y, X, Zq) {
  th <- exp(log_theta)
  V <- th[[length(th)]] * diag(length(y))
  for (i in seq_along(Zq)) {
    V <- V + th[[i]] * Zq[[i]]
  }
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) {
    return(1e10)
  }
  Vinv_X <- backsolve(R, forwardsolve(t(R), X))
  Vinv_y <- backsolve(R, forwardsolve(t(R), y))
  Rx <- tryCatch(chol(crossprod(X, Vinv_X)), error = function(e) NULL)
  if (is.null(Rx)) {
    return(1e10)
  }
  beta <- backsolve(Rx, forwardsolve(t(Rx), crossprod(X, Vinv_y)))
  res <- y - X %*% beta
  Vinv_res <- backsolve(R, forwardsolve(t(R), res))
  0.5 * (2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) + sum(res * Vinv_res))
}

# Design for one marker: all populations together. NULL when no contrast to the anchor exists.
.hebs_design <- function(d, hier, anchor) {
  batches <- sort(unique(d$batch))
  if (length(batches) < 2 || !anchor %in% batches) {
    return(NULL)
  }
  study_of <- hier$study[match(batches, hier$batch)]
  studies <- sort(unique(study_of))
  pops <- sort(unique(d$population))
  # Deviation columns for every (batch, population), so a population with no data in b has one.
  bp <- expand.grid(
    batch = batches,
    population = pops,
    stringsAsFactors = FALSE
  )
  bp_key <- paste(bp$batch, bp$population, sep = "\r")
  Zs <- outer(hier$study[match(d$batch, hier$batch)], studies, `==`) * 1
  Zb <- outer(d$batch, batches, `==`) * 1
  Zw <- outer(paste(d$batch, d$population, sep = "\r"), bp_key, `==`) * 1
  list(
    d = d,
    batches = batches,
    studies = studies,
    study_of = study_of,
    pops = pops,
    bp = bp,
    y = d$value,
    # Donor-by-population indicators; outer() also covers one level, which model.matrix refuses.
    X = outer(
      paste(d$reference_id, d$population, sep = "\r"),
      sort(unique(paste(d$reference_id, d$population, sep = "\r"))),
      `==`
    ) *
      1,
    Z = list(Zs, Zb, Zw)
  )
}

# Variance components (study, run, deviation, error) for one marker by REML.
.hebs_theta <- function(z, variance_floor) {
  Zq <- lapply(z$Z, tcrossprod)
  scale <- max(stats::var(z$y), 1e-6)
  lo <- rep(log(variance_floor), 4)
  fit <- stats::optim(
    log(scale * c(0.1, 0.1, 0.1, 0.5)),
    function(lt) .hebs_reml(lt, z$y, z$X, Zq),
    method = "L-BFGS-B",
    lower = lo,
    upper = lo + 40
  )
  exp(fit$par)
}

# BLUPs of the effect for every (batch, emitted population) at given variance components.
.hebs_blup <- function(z, th, anchor, emit_pops) {
  y <- z$y
  X <- z$X
  Z <- do.call(cbind, z$Z)
  g_diag <- rep(th[1:3], times = vapply(z$Z, ncol, integer(1)))
  V <- Z %*% (g_diag * t(Z)) + th[[4]] * diag(length(y))
  Vinv <- chol2inv(chol(V))
  P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
  GZt <- g_diag * t(Z)
  u <- GZt %*% P %*% y
  pev <- diag(g_diag) - GZt %*% P %*% t(GZt)
  ns <- ncol(z$Z[[1]])
  nb <- ncol(z$Z[[2]])
  bp_key <- paste(z$bp$batch, z$bp$population, sep = "\r")
  col_w <- function(b, k) ns + nb + match(paste(b, k, sep = "\r"), bp_key)
  grid <- expand.grid(
    batch = z$batches,
    population = emit_pops,
    stringsAsFactors = FALSE
  )
  sa <- match(z$study_of[match(anchor, z$batches)], z$studies)
  L <- matrix(0, nrow(grid), ncol(Z))
  for (i in seq_len(nrow(grid))) {
    b <- grid$batch[i]
    k <- grid$population[i]
    sb <- match(z$study_of[match(b, z$batches)], z$studies)
    L[i, sb] <- L[i, sb] + 1
    L[i, sa] <- L[i, sa] - 1
    L[i, ns + match(b, z$batches)] <- L[i, ns + match(b, z$batches)] + 1
    L[i, ns + match(anchor, z$batches)] <- L[i, ns + match(anchor, z$batches)] -
      1
    L[i, col_w(b, k)] <- L[i, col_w(b, k)] + 1
    L[i, col_w(anchor, k)] <- L[i, col_w(anchor, k)] - 1
  }
  grid$estimate <- as.vector(L %*% u)
  grid$se <- sqrt(pmax(rowSums((L %*% pev) * L), 0))
  grid
}

# Effects of the shared model for every marker.
.hebs_effects <- function(summ, hierarchy, anchor, emit_pops, variance_floor) {
  rows <- lapply(split(seq_len(nrow(summ)), summ$marker), function(ii) {
    d <- summ[ii, , drop = FALSE]
    z <- .hebs_design(d, hierarchy, anchor)
    if (is.null(z)) {
      return(NULL)
    }
    r <- .hebs_blup(z, .hebs_theta(z, variance_floor), anchor, emit_pops)
    n_don <- mapply(
      function(b, k) {
        length(unique(d$reference_id[d$batch == b & d$population == k]))
      },
      r$batch,
      r$population
    )
    n_aliq <- mapply(
      function(b, k) {
        length(unique(d$aliquot[d$batch == b & d$population == k]))
      },
      r$batch,
      r$population
    )
    n_anc <- vapply(
      r$population,
      function(k) {
        length(unique(d$aliquot[d$batch == anchor & d$population == k]))
      },
      integer(1)
    )
    tibble::tibble(
      batch = r$batch,
      population = r$population,
      marker = d$marker[[1]],
      estimate = r$estimate,
      se_donor = r$se,
      n_donors = as.integer(n_don),
      n_aliquots = as.integer(n_aliq),
      n_aliquots_anchor = as.integer(n_anc),
      support = dplyr::case_when(
        r$batch == anchor ~ "anchor",
        n_aliq == 0 ~ "shared_borrowed",
        .default = "joint_reml_shared"
      ),
      source = "fit"
    )
  })
  dplyr::bind_rows(rows)
}
