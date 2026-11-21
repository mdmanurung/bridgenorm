# Moderation of batch effects across markers (plans/2026-11-20_adaptive-moderation.md).
#
# Within one batch and population, the markers' anchor contrasts y (sampling variances s2) are
# shrunk toward a centre. The prior variance tau2 is fitted by profile maximum likelihood, as in
# ComBat-refQL fit_adaptive_eb (R/shrink-batch-effects.R, MIT licence, T. Sirchi and X. Zhang).
#   standard_eb: mu + w (y - mu),           w = tau2 / (tau2 + s2)
#   adaptive:    mu + a (y - mu),           a = w + (1 - w) z2 / (1 + z2), z2 = (y - mu)^2 / (s2 + tau2)
# centre = "fitted": mu is the precision-weighted mean over markers; "zero": mu = 0.
# With tau2 = 0, standard EB returns mu for every marker, while the adaptive weight keeps
# z2 / (1 + z2) of each deviation: a heuristic departure from the point-mass posterior, kept from
# ComBat-refQL so that an outlying marker is not pulled fully to the centre.

# Moderate one vector of effects; returns estimate, weight, tau2, centre and status.
# tau2 = NULL fits the prior variance; a number uses it as given.
.heb_moderate <- function(y, s2, method, centre, tau2 = NULL) {
  given_tau2 <- tau2
  # A marker with zero (or unknown) sampling variance cannot be weighed; it is left as it is.
  valid <- is.finite(y) & is.finite(s2) & s2 > 0
  out <- list(
    estimate = y,
    weight = rep(1, length(y)),
    tau2 = NA_real_,
    centre = NA_real_,
    status = "too_few_markers"
  )
  # A fitted tau2 or centre needs two markers; a given tau2 with a zero centre needs one.
  if (sum(valid) < if (!is.null(given_tau2) && centre == "zero") 1 else 2) {
    return(out)
  }
  yv <- y[valid]
  sv <- s2[valid]
  mean_at <- function(t) {
    if (centre == "zero") 0 else sum(yv / (sv + t)) / sum(1 / (sv + t))
  }
  profile <- function(t) {
    0.5 * sum(log(sv + t) + (yv - mean_at(t))^2 / (sv + t))
  }
  tau2 <- if (is.null(given_tau2)) {
    .heb_fit_tau2(yv, sv, profile, centre)
  } else {
    given_tau2
  }
  mu <- mean_at(tau2)
  w <- if (tau2 == 0) rep(0, length(sv)) else tau2 / (tau2 + sv)
  z2 <- (yv - mu)^2 / (sv + tau2)
  weight <- if (method == "adaptive") w + (1 - w) * z2 / (1 + z2) else w
  out$estimate[valid] <- mu + weight * (yv - mu)
  out$weight[valid] <- weight
  out$tau2 <- tau2
  out$centre <- mu
  out$status <- "ok"
  out
}

# Cell-only sampling variance of each batch-minus-anchor contrast: per donor, the variance of
# its mean of aliquot means is sum(var_mean) / k^2; the contrast averages the shared donors.
.heb_cell_variance <- function(summ, anchor) {
  donor_var <- summ |>
    dplyr::summarise(
      v = sum(.data$var_mean) / dplyr::n()^2,
      .by = c("reference_id", "batch", "population", "marker")
    )
  anchor_var <- donor_var |>
    dplyr::filter(.data$batch == anchor) |>
    dplyr::select("reference_id", "population", "marker", v_anchor = "v")
  donor_var |>
    dplyr::filter(.data$batch != anchor) |>
    dplyr::inner_join(
      anchor_var,
      by = c("reference_id", "population", "marker"),
      relationship = "many-to-one"
    ) |>
    dplyr::summarise(
      cell_var = sum(.data$v + .data$v_anchor) / dplyr::n()^2,
      .by = c("batch", "population", "marker")
    )
}

# Moderate the effects table of a fit; "none" returns it unchanged.
.heb_moderate_effects <- function(
  effects,
  summ,
  anchor,
  moderation,
  centre,
  variance
) {
  if (moderation == "none") {
    return(effects)
  }
  s2 <- if (variance == "pev") {
    effects$se_donor^2
  } else {
    cv <- .heb_cell_variance(summ, anchor)
    cv$cell_var[match(
      paste(effects$batch, effects$population, effects$marker, sep = "\r"),
      paste(cv$batch, cv$population, cv$marker, sep = "\r")
    )]
  }
  effects$estimate_raw <- effects$estimate
  effects$moderation_weight <- NA_real_
  effects$moderation_tau2 <- NA_real_
  effects$moderation_centre <- NA_real_
  effects$moderation_status <- NA_character_
  # One moderation per non-anchor batch and population, across its markers.
  key <- paste(effects$batch, effects$population, sep = "\r")
  rows <- split(which(effects$batch != anchor), key[effects$batch != anchor])
  for (ii in rows) {
    r <- .heb_moderate(effects$estimate[ii], s2[ii], moderation, centre)
    effects$estimate[ii] <- r$estimate
    effects$moderation_weight[ii] <- r$weight
    effects$moderation_tau2[ii] <- r$tau2
    effects$moderation_centre[ii] <- r$centre
    effects$moderation_status[ii] <- r$status
  }
  effects
}

# Profile-ML prior variance: log-grid scan, log-scale refine, then compared with tau2 = 0.
.heb_fit_tau2 <- function(yv, sv, profile, centre) {
  upper <- max(stats::var(yv) * 100, max(sv), 1e-8)
  # A zero centre must reach a shift shared by every marker.
  if (centre == "zero") {
    upper <- max(upper, max(yv^2))
  }
  # The profile can have several basins: scan a log grid relative to the variances, then refine
  # around the best grid point and compare with tau2 = 0.
  grid <- exp(seq(log(min(sv) * 1e-6), log(upper), length.out = 200))
  nll <- vapply(grid, profile, numeric(1))
  i <- which.min(nll)
  # Refine on the log scale, so the tolerance is relative to tau2.
  opt <- stats::optimize(
    function(lt) profile(exp(lt)),
    log(c(grid[max(i - 1, 1)], grid[min(i + 1, length(grid))])),
    tol = 1e-10
  )
  best <- if (opt$objective < nll[[i]]) exp(opt$minimum) else grid[[i]]
  tau2 <- if (profile(0) <= profile(best) + sqrt(.Machine$double.eps)) {
    0
  } else {
    best
  }
  tau2
}

# Well (aliquot) variance beyond cell sampling, per population and marker: the pooled variance of
# aliquot means within donor and batch (groups with 2+ aliquots) minus their mean var(cells)/n,
# floored at 0. Batches with one aliquot borrow it (roadmap R2, user choice R2c).
.heb_well_variance <- function(summ) {
  multi <- summ |>
    dplyr::mutate(
      dev = .data$value - mean(.data$value),
      n_al = dplyr::n(),
      .by = c("reference_id", "batch", "population", "marker")
    ) |>
    dplyr::filter(.data$n_al >= 2)
  all_pm <- dplyr::distinct(summ, .data$population, .data$marker)
  if (nrow(multi) == 0) {
    return(dplyr::mutate(all_pm, well_var = 0))
  }
  pooled <- multi |>
    dplyr::summarise(
      ss = sum(.data$dev^2),
      df = dplyr::n() -
        dplyr::n_distinct(paste(.data$reference_id, .data$batch)),
      cell = mean(.data$var_mean),
      .by = c("population", "marker")
    ) |>
    dplyr::transmute(
      .data$population,
      .data$marker,
      well_var = pmax(.data$ss / .data$df - .data$cell, 0)
    )
  dplyr::left_join(
    all_pm,
    pooled,
    by = c("population", "marker"),
    relationship = "one-to-one"
  ) |>
    dplyr::mutate(well_var = dplyr::coalesce(.data$well_var, 0))
}

# Shrinkage across batches (roadmap R2): each batch's unshrunk paired contrast y (variance
# cell + well) is moderated across the batches of its group, per population and marker.
# centre "study": toward the study's precision-weighted mean; "zero": toward the anchor.
# tau2 "profile": fitted by .heb_moderate; "reml": the batch variance of the heb fit.
.heb_batch_shrink <- function(
  effects,
  summ,
  hierarchy,
  anchor,
  method,
  centre,
  group,
  tau2_source,
  reml_tau2
) {
  if (method == "blup") {
    return(effects)
  }
  donor <- .heb_donor_level(summ)
  paired <- .heb_paired_effects(donor, .heb_targets(donor, anchor), "fit") |>
    dplyr::filter(.data$batch != anchor)
  cv <- .heb_cell_variance(summ, anchor)
  wv <- .heb_well_variance(summ)
  paired <- paired |>
    dplyr::left_join(
      cv,
      by = c("batch", "population", "marker"),
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      wv,
      by = c("population", "marker"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      k_b = .data$n_aliquots / .data$n_donors,
      k_a = .data$n_aliquots_anchor / .data$n_donors,
      s2 = .data$cell_var +
        .data$well_var * (1 / .data$k_b + 1 / .data$k_a) / .data$n_donors,
      study = hierarchy$study[match(.data$batch, hierarchy$batch)],
      grp = if (group == "study") {
        paste(.data$population, .data$marker, .data$study, sep = "\r")
      } else {
        paste(.data$population, .data$marker, sep = "\r")
      }
    )
  paired$shrunk <- paired$estimate
  paired$weight <- 1
  paired$tau2 <- NA_real_
  paired$status <- NA_character_
  for (ii in split(seq_len(nrow(paired)), paired$grp)) {
    y <- paired$estimate[ii]
    s2 <- paired$s2[ii]
    t2 <- if (tau2_source == "reml") {
      unname(reml_tau2[paste(
        paired$population[ii[1]],
        paired$marker[ii[1]],
        sep = "\r"
      )])
    } else {
      NULL
    }
    # A study centre inside an all-batch group: remove each study's precision-weighted mean first.
    mu_s <- if (centre == "study" && group == "all") {
      stats::ave(
        y / s2,
        paired$study[ii],
        FUN = sum
      ) /
        stats::ave(1 / s2, paired$study[ii], FUN = sum)
    } else {
      rep(0, length(ii))
    }
    ctr <- if (centre == "study" && group == "study") "fitted" else "zero"
    r <- .heb_moderate(y - mu_s, s2, method, ctr, tau2 = t2)
    paired$shrunk[ii] <- r$estimate + mu_s
    paired$weight[ii] <- r$weight
    paired$tau2[ii] <- r$tau2
    # Too few batches to moderate: a one-batch study group keeps its paired contrast (its own centre).
    paired$status[ii] <- if (r$status == "ok") "ok" else "single_batch_group"
  }
  key_e <- paste(effects$batch, effects$population, effects$marker, sep = "\r")
  key_p <- paste(paired$batch, paired$population, paired$marker, sep = "\r")
  m <- match(key_e, key_p)
  effects$estimate_blup <- effects$estimate
  effects$batch_shrinkage_weight <- paired$weight[m]
  effects$batch_shrinkage_tau2 <- paired$tau2[m]
  effects$batch_shrinkage_status <- paired$status[m]
  hit <- !is.na(m)
  effects$estimate[hit] <- paired$shrunk[m[hit]]
  effects
}
