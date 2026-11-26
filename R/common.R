# Checks the fit functions share: inputs, bridges and the unbridged gate.

# Checks every fit shares; returns the cell table and the bridge table.
.bn_fit_common <- function(
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
  transform <- .bn_check_transform(transform)
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
  .bn_check_min_cells(min_cells)
  if (
    !is.data.frame(hierarchy) || !all(c("batch", "study") %in% names(hierarchy))
  ) {
    stop(
      what,
      " needs hierarchy = data frame with columns 'batch' and 'study'.",
      call. = FALSE
    )
  }
  if (anyNA(hierarchy$study) || anyNA(hierarchy$batch)) {
    stop(
      what,
      ": 'hierarchy' has a missing batch or study label; every batch needs a study.",
      call. = FALSE
    )
  }
  d <- .bn_cells(x, samples, population, markers, sample_ids)
  batches <- sort(unique(d$cells$batch), method = "radix")
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
  bridges <- bn_bridges(present, anchor)
  list(d = d, transform = transform, batches = batches, bridges = bridges)
}


# Stop on unbridged batches, or return them for abstention.
.bn_gate <- function(bridges, unbridged) {
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
.bn_check_min_cells <- function(min_cells) {
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
.bn_bridges_after_min_cells <- function(batches, summ, anchor) {
  .bn_bridge_table(
    c(as.character(summ$batch), batches),
    c(as.character(summ$reference_id), rep(NA_character_, length(batches))),
    anchor
  )
}
