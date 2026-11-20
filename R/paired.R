# User-facing paired-shift fit: the baseline every heb variant must beat.

#' Fit a frozen paired-shift model
#'
#' For each batch, population and marker, the effect is the mean over the
#' references shared with the anchor batch of (reference summary in the batch
#' minus the same reference's summary in the anchor). Reference summaries are
#' equal-weight means of aliquot means, as in \code{\link{heb_fit}}. Nothing is
#' shrunk. This is the estimator \code{\link{heb_integrate}} uses for batches
#' unseen at fit time, exposed so that it can be fitted and compared directly.
#'
#' @details
#' \code{se_donor} is the SD of the paired differences over references divided
#' by the square root of their number (\code{NA} with one reference), with
#' \code{df} = number of references - 1. Anchor rows have estimate 0 and
#' \code{support = "anchor"}.
#'
#' With \code{unbridged = "abstain"}, a batch that shares no reference with the
#' anchor (after \code{min_cells}) gets no effect, so \code{\link{heb_integrate}} leaves its cells raw
#' and flags them; \code{"stop"} refuses the fit.
#'
#' @inheritParams heb_fit
#' @return A \code{heb_model} with \code{method = "paired_shift"}; see
#'   \code{\link{heb_fit}}. Its \code{hierarchy} puts every batch in one study.
#' @seealso \code{\link{heb_fit}}, \code{\link{heb_integrate}}
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
#'     CD3 = rnorm(100, 2 + shift[[samples$batch[i]]], 0.3)
#'   )
#' }))
#' model <- heb_fit_paired(
#'   x, samples, population = rep("T", nrow(x)), markers = "CD3",
#'   transform = list(name = "asinh", cofactor = 5), anchor = "b1",
#'   min_cells = 20, unbridged = "abstain"
#' )
#' model$effects[, c("batch", "estimate", "se_donor")]
#' @export
heb_fit_paired <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  min_cells,
  unbridged,
  sample_ids = NULL
) {
  .heb_stop_missing(c(
    x = missing(x),
    samples = missing(samples),
    population = missing(population),
    markers = missing(markers),
    transform = missing(transform),
    anchor = missing(anchor),
    min_cells = missing(min_cells),
    unbridged = missing(unbridged)
  ))
  .heb_check_choice(unbridged, "unbridged", c("stop", "abstain"))
  samples <- .heb_check_samples(samples)
  # Paired shift has no hierarchy; every batch sits in one study.
  hierarchy <- data.frame(
    batch = sort(unique(as.character(samples$batch))),
    study = "all"
  )
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
    "heb_fit_paired"
  )
  s <- .heb_summaries(com$d$xm, com$d$cells, as.integer(min_cells))
  if (nrow(s$summ) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }
  com$bridges <- .heb_bridges_after_min_cells(com$batches, s$summ, anchor)
  unbridged_batches <- .heb_gate(com$bridges, unbridged)
  donor <- .heb_donor_level(s$summ)
  targets <- .heb_targets(donor, anchor)
  effects <- .heb_paired_effects(donor, targets, "fit")
  effects$support[effects$batch == anchor] <- "anchor"
  effects <- effects[!effects$batch %in% unbridged_batches, , drop = FALSE]
  .heb_freeze(
    "paired_shift",
    list(unbridged = unbridged),
    com,
    anchor,
    min_cells,
    hierarchy[hierarchy$batch %in% com$batches, , drop = FALSE],
    effects,
    targets,
    s$excluded
  )
}
