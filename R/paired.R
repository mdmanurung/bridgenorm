# User-facing paired-shift fit: the baseline every bridgenorm model must beat.

#' Fit a frozen paired-shift model
#'
#' For each batch, population and marker, the effect is the mean over the
#' references shared with the anchor batch of (reference summary in the batch
#' minus the same reference's summary in the anchor). Reference summaries are
#' equal-weight means of aliquot means. Nothing is shrunk. This is the
#' baseline every bridgenorm model must beat, the estimator \code{\link{bn_fit}}
#' uses in each SOM node, and the one \code{\link{bn_apply}} uses for batches
#' unseen at fit time.
#'
#' @details
#' \code{se_donor} is the SD of the paired differences over references divided
#' by the square root of their number (\code{NA} with one reference), with
#' \code{df} = number of references - 1. Anchor rows have estimate 0 and
#' \code{support = "anchor"}. With \code{target = "virtual"} each reference
#' summary of the anchor is replaced by its mean over the bridged batches, so
#' the anchor batch gets an effect too.
#'
#' With \code{unbridged = "abstain"}, a batch that shares no reference with the
#' anchor (after \code{min_cells}) gets no effect, so \code{\link{bn_apply}} leaves its cells raw
#' and flags them; \code{"stop"} refuses the fit.
#'
#' @inheritParams bn_fit
#' @param target \code{"anchor"} corrects every batch to the anchor batch.
#'   \code{"virtual"} replaces each anchor reference summary by its mean over
#'   the bridged batches, so the anchor batch gets an effect too. The anchor
#'   still sets the bridges.
#' @return A \code{bn_model}: a list with \code{method = "paired_shift"},
#'   \code{settings}, \code{markers}, \code{transform}, \code{anchor},
#'   \code{effects} (one row per batch, population and marker: estimate,
#'   \code{se_donor}, \code{df}, donor and aliquot counts, \code{support}),
#'   \code{targets} (anchor summaries), \code{bridges}, \code{excluded},
#'   \code{known_batches} and \code{digest}.
#' @seealso \code{\link{bn_fit}}, \code{\link{bn_apply}}
#' @examples
#' set.seed(1)
#' batches <- c("b1", "b2", "b3")
#' samples <- expand.grid(donor = c("D1", "D2"), batch = batches, stringsAsFactors = FALSE)
#' samples$sample_id <- paste(samples$donor, samples$batch, sep = "_")
#' samples$reference_id <- samples$donor
#' shift <- c(b1 = 0, b2 = 0.5, b3 = -0.3)
#' x <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
#'   data.frame(
#'     sample_id = samples$sample_id[i],
#'     Lin = rnorm(100, 2 + shift[[samples$batch[i]]], 0.3)
#'   )
#' }))
#' model <- bn_fit_paired(
#'   x, samples, population = rep("T", nrow(x)), markers = "Lin",
#'   transform = list(name = "asinh", cofactor = 5), anchor = "b1",
#'   min_cells = 20
#' )
#' model$effects[, c("batch", "estimate", "se_donor")]
#' @export
bn_fit_paired <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  min_cells = 50,
  unbridged = "abstain",
  target = "anchor",
  goal_batches = NULL,
  sample_ids = NULL
) {
  .bn_fit_paired(
    x,
    samples,
    population,
    markers,
    transform,
    anchor,
    min_cells,
    unbridged,
    target,
    goal_batches,
    sample_ids
  )
}

# bn_fit_paired() with the rule for goal summaries the anchor has and the goal lacks: "stop" for user fits,
# "drop" for the node fits inside bn_fit(), whose dropped nodes take the global fallback.
.bn_fit_paired <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  min_cells = 50,
  unbridged = "abstain",
  target = "anchor",
  goal_batches = NULL,
  sample_ids = NULL,
  goal_lacking = "stop"
) {
  .bn_stop_missing(c(
    x = missing(x),
    samples = missing(samples),
    population = missing(population),
    markers = missing(markers),
    transform = missing(transform),
    anchor = missing(anchor)
  ))
  .bn_check_choice(
    unbridged,
    "unbridged",
    c("stop", "abstain"),
    what = "bn_fit_paired"
  )
  .bn_check_target(target)
  samples <- .bn_check_samples(samples)
  goal_batches <- .bn_check_goal_batches(goal_batches, target, samples)
  # Paired shift has no hierarchy; every batch sits in one study.
  hierarchy <- data.frame(
    batch = sort(unique(as.character(samples$batch))),
    study = "all"
  )
  com <- .bn_fit_common(
    x,
    samples,
    population,
    markers,
    transform,
    anchor,
    hierarchy,
    min_cells,
    sample_ids,
    "bn_fit_paired"
  )
  s <- .bn_summaries(com$d$xm, com$d$cells, as.integer(min_cells))
  if (nrow(s$summ) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }
  com$bridges <- .bn_bridges_after_min_cells(com$batches, s$summ, anchor)
  unbridged_batches <- .bn_gate(com$bridges, unbridged)
  donor <- .bn_donor_level(s$summ)
  # Unbridged batches never shape a virtual target; goal batches that are unbridged are dropped.
  targets <- .bn_targets(
    donor[!donor$batch %in% unbridged_batches, , drop = FALSE],
    anchor,
    target,
    goal = goal_batches,
    lacking = goal_lacking
  )
  effects <- .bn_paired_effects(donor, targets, "fit")
  if (target == "anchor") {
    effects$support[effects$batch == anchor] <- "anchor"
  }
  effects <- effects[!effects$batch %in% unbridged_batches, , drop = FALSE]
  .bn_freeze(
    "paired_shift",
    c(
      list(unbridged = unbridged, target = target),
      if (!is.null(goal_batches)) list(goal_batches = goal_batches)
    ),
    com,
    anchor,
    min_cells,
    effects,
    targets,
    s$excluded
  )
}
