# Input checks, the per-cell table and the anchor bridge.
#
# Ported from SpectraWeaveR R/corrector.R (.sw_stop_missing, .sw_check_transform,
# .sw_corrector_cells). The sample table replaces sw_sample_manifest(): a sample is a
# reference when its reference_id is not NA.

# Stop naming every required argument the caller left out.
.heb_stop_missing <- function(flags) {
  absent <- names(flags)[flags]
  if (length(absent) > 0) {
    stop(
      "Required argument(s) with no default: ",
      paste0("'", absent, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# One string out of a fixed set, or stop with the allowed values.
.heb_check_choice <- function(value, name, choices, what = "heb") {
  if (!is.character(value) || length(value) != 1 || !value %in% choices) {
    stop(
      sprintf(
        "%s needs %s = %s (no default).",
        what,
        name,
        paste0("'", choices, "'", collapse = " or ")
      ),
      call. = FALSE
    )
  }
  value
}

# Normalise and check list(name = "asinh", cofactor = <number>).
.heb_check_transform <- function(transform) {
  if (
    !is.list(transform) ||
      !is.character(transform[["name"]]) ||
      length(transform[["name"]]) != 1 ||
      is.na(transform[["name"]])
  ) {
    stop(
      "'transform' must be a list with a 'name' string, ",
      "for example list(name = \"asinh\", cofactor = 6000).",
      call. = FALSE
    )
  }
  out <- list(name = transform[["name"]])
  if (transform[["name"]] == "asinh") {
    cf <- transform[["cofactor"]]
    if (!is.numeric(cf) || length(cf) != 1 || is.na(cf) || cf <= 0) {
      stop(
        "An asinh 'transform' needs a positive 'cofactor'; it has no default.",
        call. = FALSE
      )
    }
    out$cofactor <- as.double(cf)
  } else if (!is.null(transform[["cofactor"]])) {
    out$cofactor <- as.double(transform[["cofactor"]])
  }
  out
}

# Short text for a transform, used in messages.
.heb_transform_label <- function(transform) {
  if (is.null(transform[["cofactor"]])) {
    transform[["name"]]
  } else {
    sprintf(
      "%s (cofactor %s)",
      transform[["name"]],
      format(transform[["cofactor"]])
    )
  }
}

# Check the sample table: one row per sample, the columns the model reads.
.heb_check_samples <- function(samples) {
  need <- c("sample_id", "batch", "reference_id")
  if (!is.data.frame(samples) || !all(need %in% names(samples))) {
    stop(
      "'samples' must be a data frame with columns sample_id, batch and ",
      "reference_id (NA for samples that are not references).",
      call. = FALSE
    )
  }
  ids <- as.character(samples$sample_id)
  if (anyNA(ids) || anyDuplicated(ids) > 0) {
    stop(
      "'samples' must have one row per sample_id, none missing.",
      call. = FALSE
    )
  }
  if (anyNA(samples$batch)) {
    stop("'samples' has a missing batch.", call. = FALSE)
  }
  samples
}

# Marker matrix and per-cell table from x, the sample table and the labels.
.heb_cells <- function(x, samples, population, markers, sample_ids) {
  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop("A matrix 'x' must be numeric.", call. = FALSE)
    }
    if (is.null(sample_ids)) {
      stop("A matrix 'x' needs 'sample_ids', one per row.", call. = FALSE)
    }
    ids <- as.character(sample_ids)
  } else if (is.data.frame(x)) {
    if (!is.null(sample_ids)) {
      ids <- as.character(sample_ids)
    } else if ("sample_id" %in% names(x)) {
      ids <- as.character(x[["sample_id"]])
    } else {
      stop(
        "'x' needs a 'sample_id' column, or pass 'sample_ids'.",
        call. = FALSE
      )
    }
  } else {
    stop("'x' must be a numeric matrix or a data frame.", call. = FALSE)
  }
  if (length(ids) != nrow(x) || anyNA(ids)) {
    stop(
      "Sample ids must be one non-missing value per row of 'x'.",
      call. = FALSE
    )
  }
  absent <- setdiff(markers, colnames(x))
  if (length(absent) > 0) {
    stop(
      "Marker(s) not in 'x': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  xm <- if (is.matrix(x)) {
    x[, markers, drop = FALSE]
  } else {
    as.matrix(as.data.frame(x)[markers])
  }
  if (!is.numeric(xm)) {
    stop("Marker columns of 'x' must be numeric.", call. = FALSE)
  }
  storage.mode(xm) <- "double"
  if (any(!is.finite(xm))) {
    stop("Marker values must be finite; found NA, NaN or Inf.", call. = FALSE)
  }
  if (
    !is.atomic(population) || length(population) != nrow(x) || anyNA(population)
  ) {
    stop(
      "'population' must be one non-missing label per row of 'x'.",
      call. = FALSE
    )
  }

  samples <- .heb_check_samples(samples)
  uncovered <- setdiff(unique(ids), as.character(samples$sample_id))
  if (length(uncovered) > 0) {
    stop(
      "Sample(s) in 'x' missing from 'samples': ",
      paste(utils::head(uncovered, 10), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  m <- match(ids, as.character(samples$sample_id))
  ref <- as.character(samples$reference_id[m])
  is_ref <- !is.na(ref)

  has_aliquot <- "aliquot" %in% names(samples)
  aliquot <- if (has_aliquot) as.character(samples$aliquot[m]) else ref
  if (has_aliquot && any(is_ref & is.na(aliquot))) {
    stop(
      "Reference sample(s) have a missing 'aliquot': ",
      paste(unique(ids[is_ref & is.na(aliquot)]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  aliquot_source <- if (has_aliquot) {
    "samples aliquot column"
  } else {
    paste0(
      "defaulted to reference_id: 'samples' has no aliquot column, so all ",
      "cells of one reference in one batch form one summary"
    )
  }

  cells <- tibble::tibble(
    row = seq_len(nrow(x)),
    sample_id = ids,
    batch = as.character(samples$batch[m]),
    population = as.character(population),
    is_reference = is_ref,
    reference_id = ifelse(is_ref, ref, NA_character_),
    aliquot = ifelse(is_ref, aliquot, NA_character_)
  )
  list(xm = xm, cells = cells, aliquot_source = aliquot_source)
}

#' Which batches share a reference with the anchor batch
#'
#' A batch is bridged when at least one of its reference samples has a
#' \code{reference_id} that is also measured in the anchor batch. Connection
#' through other batches does not count. This is the rule
#' \code{\link{heb_fit}} uses for \code{unbridged}.
#'
#' @param samples Data frame with columns \code{sample_id}, \code{batch} and
#'   \code{reference_id} (\code{NA} for samples that are not references).
#' @param anchor The anchor batch. **Required.**
#' @return A tibble with one row per non-anchor batch: \code{batch},
#'   \code{n_shared_with_anchor}, \code{shared_ids} (comma-separated) and
#'   \code{bridged}.
#'
#' @examples
#' samples <- data.frame(
#'   sample_id = paste0("s", 1:4), batch = c("b1", "b2", "b3", "b3"),
#'   reference_id = c("R1", "R1", "R2", NA)
#' )
#' heb_bridges(samples, anchor = "b1")
#'
#' @export
heb_bridges <- function(samples, anchor) {
  .heb_stop_missing(c(samples = missing(samples), anchor = missing(anchor)))
  samples <- .heb_check_samples(samples)
  batch <- as.character(samples$batch)
  ref <- as.character(samples$reference_id)
  if (!anchor %in% batch) {
    stop(
      sprintf("Anchor batch '%s' is not in 'samples'.", anchor),
      call. = FALSE
    )
  }
  .heb_bridge_table(batch, ref, anchor)
}

# Bridge table from per-row batch and reference ids (NA = not a reference).
.heb_bridge_table <- function(batch, ref, anchor) {
  anchor_ids <- unique(ref[batch == anchor & !is.na(ref)])
  others <- sort(setdiff(unique(batch), anchor))
  shared <- lapply(others, function(b) {
    sort(intersect(unique(ref[batch == b & !is.na(ref)]), anchor_ids))
  })
  tibble::tibble(
    batch = others,
    n_shared_with_anchor = lengths(shared),
    shared_ids = vapply(shared, paste, character(1), collapse = ", "),
    bridged = lengths(shared) > 0
  )
}
