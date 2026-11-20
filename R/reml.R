# Joint REML hierarchical model, one design per (population, marker).
#
# Ported from SpectraWeaveR R/corrector_heb.R (HEB-v0) without the internal anchor
# rows (plans A0-A15), which failed on real data and were dropped (user 2026-11-20).
#
# Model on aliquot summaries of reference samples:
#   y = donor (fixed) + study (random) + run/batch inside study (random) + error.
# The correction of batch b is the BLUP difference beta_b - beta_anchor.

# Negative REML log-likelihood (up to a constant) for log variance components.
# V = th1 * M1 + th2 * M2 + diag(th3 + D), with M1 = Z K1 Z' and M2 = Z K2 Z'.
# D holds known sampling variances (residual_model = "known_plus_extra"), else 0.
.heb_reml <- function(log_theta, y, X, M1, M2, D = 0) {
  th <- exp(log_theta)
  V <- th[[1]] * M1 + th[[2]] * M2 + diag(th[[3]] + D, length(y))
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) {
    return(1e10)
  }
  Vinv_X <- backsolve(R, forwardsolve(t(R), X))
  Vinv_y <- backsolve(R, forwardsolve(t(R), y))
  XtVX <- crossprod(X, Vinv_X)
  Rx <- tryCatch(chol(XtVX), error = function(e) NULL)
  if (is.null(Rx)) {
    return(1e10)
  }
  beta <- backsolve(Rx, forwardsolve(t(Rx), crossprod(X, Vinv_y)))
  res <- y - X %*% beta
  Vinv_res <- backsolve(R, forwardsolve(t(R), res))
  0.5 * (2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) + sum(res * Vinv_res))
}

# Design pieces for one (population, marker): NULL when it gives no contrast to the anchor.
# anchor_effect = "random": study and run effects of every batch are random and the anchor's
#   own effect is a BLUP; the correction is the BLUP difference to the anchor (derivation v0).
# anchor_effect = "fixed": the anchor-relative effects b_b = beta_b - beta_anchor are the random
#   effects, b_anchor = 0 exactly, with covariance
#   tau_s^2 (S_ij - S_ia - S_ja + 1) + tau_b^2 (delta_ij + 1), S = same-study indicator (v1).
.heb_design <- function(d, hier, anchor, anchor_effect, residual_model) {
  batches <- sort(unique(d$batch))
  if (length(batches) < 2 || !anchor %in% batches) {
    return(NULL)
  }
  study_of <- hier$study[match(batches, hier$batch)]
  studies <- sort(unique(study_of))
  nb <- length(batches)
  if (anchor_effect == "random") {
    Zs <- outer(hier$study[match(d$batch, hier$batch)], studies, `==`) * 1
    Zb <- outer(d$batch, batches, `==`) * 1
    Z <- cbind(Zs, Zb)
    q1 <- length(studies)
    K1 <- diag(c(rep(1, q1), rep(0, nb)))
    K2 <- diag(c(rep(0, q1), rep(1, nb)))
    L <- matrix(0, nb, ncol(Z))
    L[cbind(seq_len(nb), match(study_of, studies))] <- 1
    L[cbind(seq_len(nb), q1 + seq_len(nb))] <- 1
    L <- sweep(L, 2, L[match(anchor, batches), ], `-`)
  } else {
    nab <- batches[batches != anchor]
    s_n <- study_of[match(nab, batches)]
    s_a <- study_of[match(anchor, batches)]
    Z <- outer(d$batch, nab, `==`) * 1
    ia <- as.numeric(s_n == s_a)
    K1 <- outer(s_n, s_n, `==`) *
      1 -
      outer(ia, rep(1, length(nab))) -
      outer(rep(1, length(nab)), ia) +
      1
    K2 <- diag(length(nab)) + 1
    L <- matrix(0, nb, length(nab))
    L[cbind(match(nab, batches), seq_along(nab))] <- 1
  }
  # Donor indicators; outer() also covers a design with one donor, which model.matrix refuses.
  X <- outer(d$reference_id, sort(unique(d$reference_id)), `==`) * 1
  # Known part of each summary's residual variance: var(cells) / n (Fay-Herriot style).
  D <- if (residual_model == "known_plus_extra") d$var_mean else rep(0, nrow(d))
  list(
    d = d,
    batches = batches,
    y = d$value,
    X = X,
    Z = Z,
    K1 = K1,
    K2 = K2,
    L = L,
    M1 = Z %*% K1 %*% t(Z),
    M2 = Z %*% K2 %*% t(Z),
    D = D
  )
}

# Variance components by REML. Designs in one pooled set share the components:
# the objective is the sum of their REML criteria.
#
# variance_prior = "boundary_avoiding" adds Chung et al. (2013, PMID 24092484) gamma(shape 2,
# rate -> 0) priors on the study and run SDs: -log p(sigma_k) = -0.5 log theta_k, once per pooled
# set. It keeps them off zero when one donor per batch cannot separate them from the residual.
# The residual variance is not penalised. "none" is plain REML.
.heb_theta <- function(designs, variance_floor, variance_prior) {
  scale <- max(
    mean(vapply(designs, function(z) stats::var(z$y), numeric(1))),
    1e-6
  )
  penalty <- if (variance_prior == "boundary_avoiding") {
    function(lt) -0.5 * (lt[[1]] + lt[[2]])
  } else {
    function(lt) 0
  }
  obj <- function(lt) {
    sum(vapply(
      designs,
      function(z) .heb_reml(lt, z$y, z$X, z$M1, z$M2, z$D),
      numeric(1)
    )) +
      penalty(lt)
  }
  lo <- rep(log(variance_floor), 3)
  fit <- stats::optim(
    log(scale * c(0.1, 0.1, 0.5)),
    obj,
    method = "L-BFGS-B",
    lower = lo,
    upper = lo + 40
  )
  # Covariance of the log components from the inverse observed REML information.
  # A component pinned at the floor is treated as known: its rows and columns are zero.
  # A non-positive curvature direction gets no variance.
  free <- fit$par - lo > 1e-4
  cov_log <- matrix(0, 3, 3)
  if (any(free)) {
    H <- stats::optimHess(fit$par, obj)[free, free, drop = FALSE]
    e <- eigen((H + t(H)) / 2, symmetric = TRUE)
    pos <- e$values > 1e-10 * max(abs(e$values), 1)
    cov_log[free, free] <- e$vectors[, pos, drop = FALSE] %*%
      (t(e$vectors[, pos, drop = FALSE]) / e$values[pos])
  }
  list(theta = exp(fit$par), cov_log = cov_log)
}

# BLUPs and prediction error variance of each anchor contrast, from Henderson's
# mixed-model equations scaled by sigma^2:
#   C = [X'X, X'Z; Z'X, Z'Z + sigma^2 G^-1],  C [b; u] = [X'y; Z'y],  PEV(u) = sigma^2 (C^-1)_uu.
# The form G - G Z'PZ G cancels catastrophically once G / sigma^2 exceeds about 1e6.
# With known sampling variances D the residual is R = diag(th3 + D), and the unscaled
# equations are used: C = W' R^-1 W + [0, 0; 0, G^-1], PEV(u) = (C^-1)_uu.
.heb_pev <- function(z, th) {
  G <- th[[1]] * z$K1 + th[[2]] * z$K2
  p <- ncol(z$X)
  W <- cbind(z$X, z$Z)
  iu <- p + seq_len(ncol(z$Z))
  if (any(z$D > 0)) {
    rinv <- 1 / (th[[3]] + z$D)
    C <- crossprod(W * sqrt(rinv))
    C[iu, iu] <- C[iu, iu] + solve(G)
    R <- chol(C)
    sol <- backsolve(R, forwardsolve(t(R), crossprod(W, rinv * z$y)))
    Cinv_uu <- chol2inv(R)[iu, iu, drop = FALSE]
    return(list(u = sol[iu], var = pmax(rowSums((z$L %*% Cinv_uu) * z$L), 0)))
  }
  C <- crossprod(W)
  C[iu, iu] <- C[iu, iu] + th[[3]] * solve(G)
  R <- chol(C)
  sol <- backsolve(R, forwardsolve(t(R), crossprod(W, z$y)))
  Cinv_uu <- chol2inv(R)[iu, iu, drop = FALSE]
  list(u = sol[iu], var = pmax(th[[3]] * rowSums((z$L %*% Cinv_uu) * z$L), 0))
}

# BLUPs, prediction error and Satterthwaite df for one design at given components.
# df = 2 v^2 / (g' C g), v the contrast's prediction error variance, g its gradient in the
# log components (central differences) and C their covariance. cov_log = NULL gives df = NA.
.heb_blup <- function(z, th, anchor, cov_log = NULL) {
  p <- .heb_pev(z, th)
  df <- rep(NA_real_, length(z$batches))
  if (!is.null(cov_log)) {
    h <- 1e-4
    nk <- length(th)
    g <- vapply(
      seq_len(nk),
      function(k) {
        if (all(cov_log[k, ] == 0)) {
          return(rep(0, length(p$var)))
        }
        step <- replace(rep(0, nk), k, h)
        (.heb_pev(z, th * exp(step))$var - .heb_pev(z, th * exp(-step))$var) /
          (2 * h)
      },
      numeric(length(p$var))
    )
    g <- matrix(g, nrow = length(p$var))
    denom <- rowSums((g %*% cov_log) * g)
    df <- ifelse(denom > 0, 2 * p$var^2 / denom, Inf)
    df[z$batches == anchor] <- NA_real_
  }
  list(
    batch = z$batches,
    estimate = as.vector(z$L %*% p$u),
    se = sqrt(p$var),
    df = df
  )
}

# Fit pooled variance components, then the BLUP contrasts of every design (named list).
.heb_fit_designs <- function(
  designs,
  pool,
  anchor,
  variance_floor,
  variance_prior
) {
  pool_key <- switch(
    pool,
    none = names(designs),
    population = vapply(designs, function(z) z$d$marker[[1]], character(1)),
    all = rep("all", length(designs))
  )
  names(pool_key) <- names(designs)
  theta <- lapply(
    split(designs, pool_key),
    .heb_theta,
    variance_floor = variance_floor,
    variance_prior = variance_prior
  )
  lapply(stats::setNames(names(designs), names(designs)), function(nm) {
    tk <- theta[[pool_key[[nm]]]]
    z <- designs[[nm]]
    # A singular system (e.g. one reference per batch with near-zero shifts) stops with its design named.
    tryCatch(
      c(.heb_blup(z, tk$theta, anchor, tk$cov_log), list(theta = tk$theta)),
      error = function(e) {
        stop(
          sprintf(
            paste(
              "The mixed-model system for population '%s', marker '%s' is not positive definite",
              "(variance components %s): %s. Try anchor_effect = \"fixed\" or",
              "residual_model = \"known_plus_extra\"."
            ),
            z$d$population[[1]],
            z$d$marker[[1]],
            paste(signif(tk$theta, 3), collapse = ", "),
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  })
}

# scale_model = "log_sd" (plan F12): the same model on the log SD of each reference aliquot
# summary (known sampling variance 1 / (2 (n - 1))). delta_b is the BLUP log-SD contrast to the
# anchor; the correction multiplies deviations from the batch's shifted reference median by
# scale = exp(-delta_b). Batches without reference cells in a population get no scale (NA).
.heb_scale <- function(
  effects,
  xm,
  cells,
  summ,
  hierarchy,
  anchor,
  anchor_effect,
  residual_model,
  pool,
  variance_floor,
  variance_prior
) {
  ref <- cells[cells$is_reference, , drop = FALSE]
  ref_key <- paste(
    ref$reference_id,
    ref$aliquot,
    ref$batch,
    ref$population,
    sep = "\r"
  )
  rows_of <- split(ref$row, ref_key)
  summ_key <- paste(
    summ$reference_id,
    summ$aliquot,
    summ$batch,
    summ$population,
    sep = "\r"
  )
  ls <- summ
  ls$value <- vapply(
    seq_len(nrow(summ)),
    function(i) log(stats::sd(xm[rows_of[[summ_key[[i]]]], summ$marker[[i]]])),
    numeric(1)
  )
  ls$var_mean <- 1 / (2 * (summ$n_cells - 1))
  ls <- ls[is.finite(ls$value), , drop = FALSE]
  key <- paste(ls$population, ls$marker, sep = "\r")
  designs <- lapply(split(seq_len(nrow(ls)), key), function(ii) {
    .heb_design(
      ls[ii, , drop = FALSE],
      hierarchy,
      anchor,
      anchor_effect,
      residual_model
    )
  })
  designs <- designs[!vapply(designs, is.null, logical(1))]
  blups <- .heb_fit_designs(
    designs,
    pool,
    anchor,
    variance_floor,
    variance_prior
  )
  sc <- dplyr::bind_rows(lapply(names(designs), function(nm) {
    tibble::tibble(
      batch = blups[[nm]]$batch,
      population = designs[[nm]]$d$population[[1]],
      marker = designs[[nm]]$d$marker[[1]],
      log_sd_effect = blups[[nm]]$estimate,
      scale_se = blups[[nm]]$se
    )
  }))
  # Centre: median of the batch's reference cells for that population, minus the shift.
  # Only aliquots that passed min_cells count.
  ref <- ref[ref_key %in% summ_key, , drop = FALSE]
  grp_key <- paste(ref$batch, ref$population, sep = "\r")
  grp <- split(ref$row, grp_key)
  centre <- dplyr::bind_rows(lapply(names(grp), function(k) {
    tibble::tibble(
      batch = ref$batch[match(k, grp_key)],
      population = ref$population[match(k, grp_key)],
      marker = colnames(xm),
      ref_median = apply(xm[grp[[k]], , drop = FALSE], 2, stats::median)
    )
  }))
  out <- dplyr::left_join(
    effects,
    sc,
    by = c("batch", "population", "marker"),
    relationship = "one-to-one"
  )
  out <- dplyr::left_join(
    out,
    centre,
    by = c("batch", "population", "marker"),
    relationship = "one-to-one"
  )
  ok <- !is.na(out$log_sd_effect) & !is.na(out$ref_median)
  out$scale <- unname(ifelse(ok, exp(-out$log_sd_effect), NA_real_))
  out$scale_centre <- unname(ifelse(
    ok,
    out$ref_median - out$estimate,
    NA_real_
  ))
  out$scale_se <- unname(out$scale_se)
  out$scale_se[!ok] <- NA_real_
  out$log_sd_effect <- NULL
  out$ref_median <- NULL
  out
}
