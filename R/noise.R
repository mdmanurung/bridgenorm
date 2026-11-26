# Reference-noise report per batch (roadmap R7). Report only: no estimate depends on it.

#' Report the reference noise of each batch
#'
#' For each batch, population and marker, summarises how noisy and how
#' unusual the reference samples are, so the user can choose the anchor batch.
#' It never chooses the anchor itself.
#'
#' @details
#' Columns:
#' \describe{
#'   \item{\code{n_references}, \code{n_aliquots}}{Reference donors and
#'     aliquots with at least \code{min_cells} cells.}
#'   \item{\code{replicate_sd}}{Pooled within-donor SD of aliquot means:
#'     \eqn{\sqrt{\sum (m_a - \bar m_r)^2 / (n_{aliquots} - n_{references})}}.
#'     \code{NA} when no donor has two aliquots in the batch.}
#'   \item{\code{cell_sd}}{Mean over aliquots of the SD of their cells.}
#'   \item{\code{loo_z}}{Leave-one-batch-out z of the batch's donor mean,
#'     averaged over donors: (donor mean in the batch - its mean over the other
#'     batches) / its SD over the other batches. A donor needs two other
#'     batches; \code{NA} when none has.}
#' }
#' A donor mean is the equal-weight mean of its aliquot means, as in
#' \code{\link{bn_fit_paired}}.
#'
#' @inheritParams bn_fit
#' @param transform Optional \code{list(name = "asinh", cofactor = <number>)}:
#'   the scale \code{x} is on. It does not change the report; when given, it
#'   is checked and stored as the \code{transform} attribute of the result.
#' @return A tibble with one row per batch, population and marker.
#' @example inst/examples/bn-toy.R
#' @examples
#' bn_reference_noise(x, samples, rep("T", nrow(x)), markers = "M")
#' @seealso \code{\link{bn_fit}}, \code{\link{bn_bridges}}
#' @export
bn_reference_noise <- function(
  x,
  samples,
  population,
  markers,
  transform = NULL,
  min_cells = 50,
  sample_ids = NULL
) {
  .bn_stop_missing(c(
    x = missing(x),
    samples = missing(samples),
    population = missing(population),
    markers = missing(markers)
  ))
  if (!is.null(transform)) {
    transform <- .bn_check_transform(transform)
  }
  .bn_check_min_cells(min_cells)
  d <- .bn_cells(x, samples, population, markers, sample_ids)
  s <- .bn_summaries(d$xm, d$cells, as.integer(min_cells))$summ
  if (nrow(s) == 0) {
    stop("No reference summary has at least min_cells cells.", call. = FALSE)
  }

  # Aliquot level: replicate deviations within donor and cell SDs.
  aliquot_level <- s |>
    dplyr::mutate(
      dev = .data$value - mean(.data$value),
      .by = c("reference_id", "batch", "population", "marker")
    ) |>
    dplyr::summarise(
      n_references = dplyr::n_distinct(.data$reference_id),
      n_aliquots = dplyr::n(),
      replicate_sd = if (dplyr::n() > dplyr::n_distinct(.data$reference_id)) {
        sqrt(
          sum(.data$dev^2) /
            (dplyr::n() - dplyr::n_distinct(.data$reference_id))
        )
      } else {
        NA_real_
      },
      cell_sd = mean(sqrt(.data$var_mean * .data$n_cells)),
      .by = c("batch", "population", "marker")
    )

  # Donor level: leave-one-batch-out z per donor, then the mean over donors.
  donor <- .bn_donor_level(s)
  loo <- donor |>
    dplyr::mutate(
      z = vapply(
        seq_along(.data$value),
        function(i) {
          rest <- .data$value[-i]
          if (length(rest) < 2) {
            NA_real_
          } else {
            (.data$value[[i]] - mean(rest)) / stats::sd(rest)
          }
        },
        numeric(1)
      ),
      .by = c("reference_id", "population", "marker")
    ) |>
    dplyr::summarise(
      loo_z = if (all(is.na(.data$z))) {
        NA_real_
      } else {
        mean(.data$z, na.rm = TRUE)
      },
      .by = c("batch", "population", "marker")
    )

  out <- dplyr::left_join(
    aliquot_level,
    loo,
    by = c("batch", "population", "marker"),
    relationship = "one-to-one"
  ) |>
    dplyr::arrange(.data$batch, .data$population, .data$marker)
  attr(out, "transform") <- transform
  out
}
