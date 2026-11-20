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
  if (
    !is.numeric(min_cells) ||
      length(min_cells) != 1 ||
      is.na(min_cells) ||
      min_cells < 1 ||
      min_cells != round(min_cells)
  ) {
    stop("'min_cells' must be one whole number of at least 1.", call. = FALSE)
  }
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
#'     refuses the fit when any batch is unbridged. \code{"abstain"} fits
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
#'   residual_model = "known_plus_extra", scale_model = "none"
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
    scale_model = missing(scale_model)
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
  unbridged_batches <- .heb_gate(com$bridges, unbridged)

  s <- .heb_summaries(com$d$xm, com$d$cells, as.integer(min_cells))
  if (nrow(s$summ) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }
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
    scale_model = scale_model
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
  .heb_gate(com$bridges, "stop")
  s <- .heb_summaries(com$d$xm, com$d$cells, as.integer(min_cells))
  if (nrow(s$summ) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }
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
