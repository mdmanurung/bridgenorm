# User-facing fit functions: heb_fit() and heb_fit_shared().

# Checks both fit functions share; returns the cell table and the bridge table.
.heb_fit_common <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  hierarchy,
  min_cells,
  sample_ids,
  what
) {
  transform <- .heb_check_transform(transform)
  if (
    !is.character(markers) ||
      length(markers) == 0 ||
      anyNA(markers) ||
      anyDuplicated(markers) > 0
  ) {
    stop("'markers' must be unique, non-missing strings.", call. = FALSE)
  }
  if (!is.character(anchor) || length(anchor) != 1 || is.na(anchor)) {
    stop("'anchor' must be one batch name.", call. = FALSE)
  }
  .heb_check_min_cells(min_cells)
  if (
    !is.data.frame(hierarchy) || !all(c("batch", "study") %in% names(hierarchy))
  ) {
    stop(
      what,
      " needs hierarchy = data frame with columns 'batch' and 'study'.",
      call. = FALSE
    )
  }
  d <- .heb_cells(x, samples, population, markers, sample_ids)
  batches <- sort(unique(d$cells$batch))
  if (!anchor %in% batches) {
    stop(
      sprintf(
        "Anchor batch '%s' is not in 'x'. Batches: %s.",
        anchor,
        paste(batches, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(batches) < 2) {
    stop("'x' has one batch; there is nothing to correct.", call. = FALSE)
  }
  absent <- setdiff(batches, as.character(hierarchy$batch))
  if (length(absent) > 0 || anyDuplicated(hierarchy$batch) > 0) {
    stop(
      "'hierarchy' needs one row per batch in 'x'; missing: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  # Bridges from the samples present in x.
  present <- samples[
    as.character(samples$sample_id) %in% d$cells$sample_id,
    ,
    drop = FALSE
  ]
  bridges <- heb_bridges(present, anchor)
  list(d = d, transform = transform, batches = batches, bridges = bridges)
}

#' Fit a frozen HEB batch-correction model
#'
#' For each (population, marker), fits a linear mixed model to aliquot
#' summaries of the reference samples: donor effects fixed, study and run
#' (batch) effects random, variance components by REML. The effect of batch
#' \eqn{b} is the difference of the best linear unbiased predictions of the
#' batch effects of \eqn{b} and of the anchor, so the anchor's effect is 0.
#' \code{\link{heb_integrate}} subtracts it: \code{corrected = x - effect}.
#'
#' @details
#' \strong{Data.} \code{x} must already be on the scale named in
#' \code{transform}; the transform is recorded, never applied. A sample is a
#' reference when its \code{reference_id} in \code{samples} is not \code{NA}.
#' Only reference cells are evidence. The replicate unit is the reference
#' aliquot (\code{samples$aliquot}, or \code{reference_id} when that column is
#' absent), never the cell. Aliquot summaries with fewer than
#' \code{min_cells} cells are left out and listed in \code{excluded}.
#'
#' \strong{Settings.} All are required; none has a default.
#' \describe{
#'   \item{\code{variance_floor}}{Lower bound for every variance component.}
#'   \item{\code{pool}}{Which (population, marker) designs share one set of
#'     variance components: \code{"none"}, \code{"population"} (one set per
#'     marker, across populations) or \code{"all"}.}
#'   \item{\code{anchor_effect}}{\code{"random"}: the anchor's own effect is
#'     a BLUP and the correction is the BLUP difference. \code{"fixed"}: the
#'     anchor-relative effects are the random effects and the anchor is 0
#'     exactly.}
#'   \item{\code{variance_prior}}{\code{"none"} for plain REML, or
#'     \code{"boundary_avoiding"} for the gamma(2) prior of Chung et al.
#'     (2013) on the study and run SDs, which keeps them off zero when one
#'     donor per batch cannot separate them from the residual.}
#'   \item{\code{residual_model}}{\code{"estimated"} for one residual
#'     variance, or \code{"known_plus_extra"}: each summary's residual is its
#'     known sampling variance, var(cells) / n, plus one estimated extra
#'     variance.}
#'   \item{\code{moderation}}{\code{"none"}, or shrink each batch's
#'     marker effects (one population at a time) toward a centre:
#'     \code{"standard_eb"} with weight \eqn{w = \tau^2/(\tau^2 + s^2)}, or
#'     \code{"adaptive"} with \eqn{a = w + (1 - w) z^2/(1 + z^2)},
#'     \eqn{z^2 = (y - \mu)^2/(s^2 + \tau^2)}, which shrinks outlying markers
#'     less (ComBat-refQL, Sirchi and Zhang). \eqn{\tau^2} is fitted by profile
#'     maximum likelihood over markers. Moderated fits keep \code{estimate_raw}
#'     and add \code{moderation_weight}, \code{moderation_tau2},
#'     \code{moderation_centre} and \code{moderation_status}; \code{se_donor}
#'     and \code{df} stay those of the raw estimate. A batch with fewer than two
#'     usable markers is left as it is.}
#'   \item{\code{batch_shrinkage}}{\code{"blup"} keeps the mixed-model
#'     BLUP. \code{"standard_eb"} or \code{"adaptive"} replace it by the
#'     unshrunk paired contrast moderated across batches (the formulas of
#'     \code{moderation}), with sampling variance = cell variance plus pooled
#'     well variance (from batches with 2+ aliquots) over the aliquot counts.
#'     The BLUP is kept as \code{estimate_blup}; \code{batch_shrinkage_weight}
#'     and \code{batch_shrinkage_tau2} are added. Runs before
#'     \code{moderation}.}
#'   \item{\code{batch_shrinkage_centre}}{\code{"study"}: the study's
#'     precision-weighted mean contrast; \code{"zero"}: the anchor.}
#'   \item{\code{batch_shrinkage_group}}{\code{"study"}: batches of one
#'     study share a prior; \code{"all"}: every non-anchor batch does.}
#'   \item{\code{batch_shrinkage_tau2}}{\code{"profile"}: fitted by profile
#'     maximum likelihood over the group; \code{"reml"}: the batch variance
#'     component of the mixed model.}
#'   \item{\code{moderation_centre}}{\code{"fitted"}: the precision-weighted
#'     mean over markers, which keeps a shift shared by all markers;
#'     \code{"zero"}: no batch effect.}
#'   \item{\code{moderation_variance}}{\eqn{s^2}: \code{"pev"}, the squared
#'     \code{se_donor} of the fit, or \code{"cell"}, the cell-only variance of
#'     the contrast (aliquot var/n over the batch and anchor references).}
#'   \item{\code{scale_model}}{\code{"none"} for a pure shift, or
#'     \code{"log_sd"}: the same model is fitted to the log SD of each
#'     reference aliquot (known variance \eqn{1/(2(n-1))}), and the cells of
#'     batch \eqn{b} become \eqn{c + e^{-\delta_b}(x - \mathrm{shift} - c)},
#'     with \eqn{\delta_b} the log-SD contrast to the anchor and \eqn{c} the
#'     batch's shifted reference median. Effects gain \code{scale},
#'     \code{scale_se} and \code{scale_centre}. New batches at integrate time
#'     get the shift only.}
#'   \item{\code{unbridged}}{A batch is bridged when it shares a reference id
#'     with the anchor batch (\code{\link{heb_bridges}}). \code{"stop"}
#'     refuses the fit when any batch is unbridged after \code{min_cells}. \code{"abstain"} fits
#'     the data and drops the effects of unbridged batches, so
#'     \code{\link{heb_integrate}} leaves their cells raw and flags them.}
#' }
#'
#' \strong{Uncertainty.} \code{se_donor} is the posterior standard error of
#' the effect, from the prediction error covariance of Henderson's
#' mixed-model equations. \code{df} is its Satterthwaite degrees of freedom
#' from the inverse REML information, with components at the floor treated
#' as known (\code{Inf} when all are).
#'
#' \strong{Population gate.} A (batch, population) whose reference donors
#' share no donor with the anchor batch in that population gets no effect, so
#' \code{\link{heb_integrate}} abstains on it. \code{min_cells} must be at
#' least 2, so every summary has a sampling variance.
#'
#' \strong{Limits.} An additive shift on asinh data equals a gain change only
#' for bright populations. A population with no reference summary in the
#' anchor batch gets no effect and is abstained.
#'
#' @param x A numeric matrix (cells by markers) or a data frame with one row
#'   per cell, already transformed. A data frame needs a \code{sample_id}
#'   column unless \code{sample_ids} is given.
#' @param samples Data frame, one row per sample: \code{sample_id},
#'   \code{batch}, \code{reference_id} (\code{NA} when not a reference) and
#'   optionally \code{aliquot}.
#' @param population Per-cell population labels, one per row of \code{x},
#'   frozen before correction.
#' @param markers Marker columns to correct.
#' @param transform \code{list(name = "asinh", cofactor = <number>)}: the
#'   scale \code{x} is on.
#' @param anchor The batch that defines the target scale.
#' @param hierarchy Data frame with columns \code{batch} and \code{study},
#'   one row per batch.
#' @param min_cells Minimum cells per reference aliquot summary.
#' @param unbridged \code{"stop"} or \code{"abstain"}; see Details.
#' @param variance_floor,pool,anchor_effect,variance_prior,residual_model,scale_model
#'   Model settings; see Details.
#' @param moderation,moderation_centre,moderation_variance Moderation of each
#'   batch's effects across markers; see Details.
#' @param batch_shrinkage,batch_shrinkage_centre,batch_shrinkage_group,batch_shrinkage_tau2
#'   Shrinkage across batches (roadmap R2); see Details.
#' @param sample_ids Optional sample id per row; required for a matrix
#'   \code{x}.
#'
#' @return A \code{heb_model}: a list with \code{schema_version},
#'   \code{method} (\code{"heb"}), \code{settings}, \code{markers},
#'   \code{transform}, \code{population_levels}, \code{anchor},
#'   \code{min_cells}, \code{hierarchy}, \code{known_batches},
#'   \code{bridges}, \code{effects} (batch, population, marker, estimate,
#'   se_donor, df, n_donors, n_aliquots, n_aliquots_anchor, support, source),
#'   \code{targets}, \code{excluded}, \code{aliquot_source},
#'   \code{software}, \code{created}, \code{parent_digest} and
#'   \code{digest}.
#'
#' @seealso \code{\link{heb_fit_shared}}, \code{\link{heb_integrate}},
#'   \code{\link{heb_save}}
#'
#' @examples
#' set.seed(1)
#' batches <- c("b1", "b2", "b3")
#' shift <- c(b1 = 0, b2 = 0.5, b3 = -0.3)
#' samples <- expand.grid(donor = c("D1", "D2"), batch = batches, stringsAsFactors = FALSE)
#' samples$sample_id <- paste(samples$donor, samples$batch, sep = "_")
#' samples$reference_id <- samples$donor
#' x <- data.frame(
#'   sample_id = rep(samples$sample_id, each = 50),
#'   CD3 = rnorm(300, 2) + rep(shift[samples$batch], each = 50)
#' )
#' model <- heb_fit(
#'   x, samples, population = rep("T", 300), markers = "CD3",
#'   transform = list(name = "asinh", cofactor = 6000), anchor = "b1",
#'   hierarchy = data.frame(batch = batches, study = "s1"), min_cells = 10,
#'   unbridged = "stop", variance_floor = 1e-8, pool = "none",
#'   anchor_effect = "random", variance_prior = "boundary_avoiding",
#'   residual_model = "known_plus_extra", scale_model = "none",
#'   moderation = "none", moderation_centre = "fitted", moderation_variance = "pev",
#'   batch_shrinkage = "blup", batch_shrinkage_centre = "study",
#'   batch_shrinkage_group = "study", batch_shrinkage_tau2 = "profile"
#' )
#' model$effects
#'
#' @export
heb_fit <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  hierarchy,
  min_cells,
  unbridged,
  variance_floor,
  pool,
  anchor_effect,
  variance_prior,
  residual_model,
  scale_model,
  moderation,
  moderation_centre,
  moderation_variance,
  batch_shrinkage,
  batch_shrinkage_centre,
  batch_shrinkage_group,
  batch_shrinkage_tau2,
  sample_ids = NULL
) {
  .heb_stop_missing(c(
    x = missing(x),
    samples = missing(samples),
    population = missing(population),
    markers = missing(markers),
    transform = missing(transform),
    anchor = missing(anchor),
    hierarchy = missing(hierarchy),
    min_cells = missing(min_cells),
    unbridged = missing(unbridged),
    variance_floor = missing(variance_floor),
    pool = missing(pool),
    anchor_effect = missing(anchor_effect),
    variance_prior = missing(variance_prior),
    residual_model = missing(residual_model),
    scale_model = missing(scale_model),
    moderation = missing(moderation),
    moderation_centre = missing(moderation_centre),
    moderation_variance = missing(moderation_variance),
    batch_shrinkage = missing(batch_shrinkage),
    batch_shrinkage_centre = missing(batch_shrinkage_centre),
    batch_shrinkage_group = missing(batch_shrinkage_group),
    batch_shrinkage_tau2 = missing(batch_shrinkage_tau2)
  ))
  if (
    !is.numeric(variance_floor) ||
      length(variance_floor) != 1 ||
      !(variance_floor > 0)
  ) {
    stop("heb needs variance_floor > 0 (no default).", call. = FALSE)
  }
  .heb_check_choice(unbridged, "unbridged", c("stop", "abstain"))
  .heb_check_choice(pool, "pool", c("none", "population", "all"))
  .heb_check_choice(anchor_effect, "anchor_effect", c("random", "fixed"))
  .heb_check_choice(
    variance_prior,
    "variance_prior",
    c("none", "boundary_avoiding")
  )
  .heb_check_choice(
    residual_model,
    "residual_model",
    c("estimated", "known_plus_extra")
  )
  .heb_check_choice(scale_model, "scale_model", c("none", "log_sd"))
  .heb_check_choice(
    moderation,
    "moderation",
    c("none", "standard_eb", "adaptive")
  )
  .heb_check_choice(moderation_centre, "moderation_centre", c("fitted", "zero"))
  .heb_check_choice(
    moderation_variance,
    "moderation_variance",
    c("pev", "cell")
  )
  .heb_check_choice(
    batch_shrinkage,
    "batch_shrinkage",
    c("blup", "standard_eb", "adaptive")
  )
  .heb_check_choice(
    batch_shrinkage_centre,
    "batch_shrinkage_centre",
    c("study", "zero")
  )
  .heb_check_choice(
    batch_shrinkage_group,
    "batch_shrinkage_group",
    c("study", "all")
  )
  .heb_check_choice(
    batch_shrinkage_tau2,
    "batch_shrinkage_tau2",
    c("profile", "reml")
  )
  # A summary's sampling variance needs two cells.
  if (
    is.numeric(min_cells) &&
      length(min_cells) == 1 &&
      !is.na(min_cells) &&
      min_cells < 2
  ) {
    stop("heb_fit needs min_cells of at least 2.", call. = FALSE)
  }
  com <- .heb_fit_common(
    x,
    samples,
    population,
    markers,
    transform,
    anchor,
    hierarchy,
    min_cells,
    sample_ids,
    "heb"
  )
  s <- .heb_summaries(com$d$xm, com$d$cells, as.integer(min_cells))
  if (nrow(s$summ) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }
  com$bridges <- .heb_bridges_after_min_cells(com$batches, s$summ, anchor)
  unbridged_batches <- .heb_gate(com$bridges, unbridged)
  targets <- .heb_targets(.heb_donor_level(s$summ), anchor)
  key <- paste(s$summ$population, s$summ$marker, sep = "\r")
  designs <- lapply(split(seq_len(nrow(s$summ)), key), function(ii) {
    .heb_design(
      s$summ[ii, , drop = FALSE],
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
  effects <- dplyr::bind_rows(lapply(names(designs), function(nm) {
    r <- blups[[nm]]
    d <- designs[[nm]]$d
    n_don <- vapply(
      r$batch,
      function(b) length(unique(d$reference_id[d$batch == b])),
      integer(1)
    )
    n_aliq <- vapply(
      r$batch,
      function(b) length(unique(d$aliquot[d$batch == b])),
      integer(1)
    )
    tibble::tibble(
      batch = r$batch,
      population = d$population[[1]],
      marker = d$marker[[1]],
      estimate = r$estimate,
      se_donor = r$se,
      df = r$df,
      n_donors = unname(n_don),
      n_aliquots = unname(n_aliq),
      n_aliquots_anchor = unname(n_aliq[match(anchor, r$batch)]),
      support = ifelse(r$batch == anchor, "anchor", "joint_reml"),
      source = "fit"
    )
  }))
  # A batch-population whose donors share no reference with the anchor in that population has no
  # contrast to the anchor in the data; it gets no effect (abstained at integrate time).
  ids <- unique(s$summ[, c("batch", "population", "reference_id")])
  anchor_ids <- ids[ids$batch == anchor, c("population", "reference_id")]
  shared <- unique(dplyr::semi_join(
    ids,
    anchor_ids,
    by = c("population", "reference_id")
  )[, c("batch", "population")])
  effects <- dplyr::semi_join(effects, shared, by = c("batch", "population"))
  # Shrinkage across batches (R2) replaces the BLUP shift; REML batch variance per design.
  reml_tau2 <- vapply(blups, function(r) r$theta[[2]], numeric(1))
  effects <- .heb_batch_shrink(
    effects,
    s$summ,
    hierarchy,
    anchor,
    batch_shrinkage,
    batch_shrinkage_centre,
    batch_shrinkage_group,
    batch_shrinkage_tau2,
    reml_tau2
  )
  # Moderation across markers acts on the shift before the scale centre is computed.
  effects <- .heb_moderate_effects(
    effects,
    s$summ,
    anchor,
    moderation,
    moderation_centre,
    moderation_variance
  )
  if (scale_model == "log_sd") {
    effects <- .heb_scale(
      effects,
      com$d$xm,
      com$d$cells,
      s$summ,
      hierarchy,
      anchor,
      anchor_effect,
      residual_model,
      pool,
      variance_floor,
      variance_prior
    )
  }
  effects <- effects[!effects$batch %in% unbridged_batches, , drop = FALSE]
  settings <- list(
    unbridged = unbridged,
    variance_floor = variance_floor,
    pool = pool,
    anchor_effect = anchor_effect,
    variance_prior = variance_prior,
    residual_model = residual_model,
    scale_model = scale_model,
    moderation = moderation,
    moderation_centre = moderation_centre,
    moderation_variance = moderation_variance,
    batch_shrinkage = batch_shrinkage,
    batch_shrinkage_centre = batch_shrinkage_centre,
    batch_shrinkage_group = batch_shrinkage_group,
    batch_shrinkage_tau2 = batch_shrinkage_tau2
  )
  .heb_freeze(
    "heb",
    settings,
    com,
    anchor,
    min_cells,
    hierarchy,
    effects,
    targets,
    s$excluded
  )
}

#' Fit a frozen HEB model with a batch effect shared across populations
#'
#' For each marker, fits one mixed model over all populations:
#' donor-by-population effects fixed; study and run effects random and
#' shared by every population; a batch-by-population deviation random with
#' its own variance. The effect of (batch, population) is the BLUP
#' difference to the anchor of the study, run and deviation parts.
#'
#' @details
#' \strong{Borrowing.} A population with no summary in a batch (every
#' aliquot below \code{min_cells}) still gets the shared study and run effect
#' of that batch, with zero deviation. Its rows carry
#' \code{support = "shared_borrowed"} and \code{n_donors = 0}.
#'
#' \strong{Abstention.} A population with no summary in the anchor batch has
#' no donor-by-population target. It gets no effect, so
#' \code{\link{heb_integrate}} abstains on it. Its data still inform the
#' shared effects.
#'
#' Every batch must share a reference id with the anchor batch
#' (\code{\link{heb_bridges}}); the fit stops otherwise. This model has no
#' \code{df} column, no prior, no known-variance residual and no scale term.
#'
#' @inheritParams heb_fit
#' @param variance_floor Lower bound for every variance component.
#'   **Required.**
#' @return A \code{heb_model} with \code{method = "heb_shared"}; see
#'   \code{\link{heb_fit}}.
#' @seealso \code{\link{heb_fit}}, \code{\link{heb_integrate}}
#' @export
heb_fit_shared <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  hierarchy,
  min_cells,
  variance_floor,
  sample_ids = NULL
) {
  .heb_stop_missing(c(
    x = missing(x),
    samples = missing(samples),
    population = missing(population),
    markers = missing(markers),
    transform = missing(transform),
    anchor = missing(anchor),
    hierarchy = missing(hierarchy),
    min_cells = missing(min_cells),
    variance_floor = missing(variance_floor)
  ))
  if (
    !is.numeric(variance_floor) ||
      length(variance_floor) != 1 ||
      !(variance_floor > 0)
  ) {
    stop("heb_shared needs variance_floor > 0 (no default).", call. = FALSE)
  }
  com <- .heb_fit_common(
    x,
    samples,
    population,
    markers,
    transform,
    anchor,
    hierarchy,
    min_cells,
    sample_ids,
    "heb_shared"
  )
  s <- .heb_summaries(com$d$xm, com$d$cells, as.integer(min_cells))
  if (nrow(s$summ) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }
  com$bridges <- .heb_bridges_after_min_cells(com$batches, s$summ, anchor)
  .heb_gate(com$bridges, "stop")
  targets <- .heb_targets(.heb_donor_level(s$summ), anchor)
  # Populations that may be corrected: a donor-by-population target exists in the anchor.
  emit_pops <- sort(unique(targets$population))
  effects <- .hebs_effects(s$summ, hierarchy, anchor, emit_pops, variance_floor)
  settings <- list(variance_floor = variance_floor)
  .heb_freeze(
    "heb_shared",
    settings,
    com,
    anchor,
    min_cells,
    hierarchy,
    effects,
    targets,
    s$excluded
  )
}

# Stop on unbridged batches, or return them for abstention.
.heb_gate <- function(bridges, unbridged) {
  out <- bridges$batch[!bridges$bridged]
  if (length(out) > 0 && unbridged == "stop") {
    stop(
      "Batch(es) share no reference with the anchor batch: ",
      paste(out, collapse = ", "),
      ". Use unbridged = \"abstain\" to fit the others and leave these raw.",
      call. = FALSE
    )
  }
  out
}

# min_cells must be one whole number of at least 1.
.heb_check_min_cells <- function(min_cells) {
  if (
    !is.numeric(min_cells) ||
      length(min_cells) != 1 ||
      is.na(min_cells) ||
      min_cells < 1 ||
      min_cells != round(min_cells)
  ) {
    stop("'min_cells' must be one whole number of at least 1.", call. = FALSE)
  }
}

# Bridges recomputed from the summaries that passed min_cells: a batch whose shared references all fall
# below min_cells is unbridged. Batches with no summary keep a row with no shared id.
.heb_bridges_after_min_cells <- function(batches, summ, anchor) {
  .heb_bridge_table(
    c(as.character(summ$batch), batches),
    c(as.character(summ$reference_id), rep(NA_character_, length(batches))),
    anchor
  )
}
