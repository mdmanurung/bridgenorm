# Apply a frozen model: bn_apply().
#
# Ported from SpectraWeaveR R/corrector.R (sw_corrector_integrate, .sw_apply_effects).

#' Apply a frozen bridgenorm model to new or old data
#'
#' Corrects \code{x} with the frozen effects and reports one action per
#' (batch, population, marker). Earlier effects never change, so values
#' corrected before stay the same when new batches arrive.
#'
#' Arrival cases:
#' \describe{
#'   \item{Known batch}{Use the frozen effect. Reference samples in
#'     \code{x} are not used again.}
#'   \item{New batch with a bridging reference}{A reference in the new batch
#'     also has a frozen summary in the anchor batch. Its effect is the mean
#'     paired difference of its own reference cells to the frozen anchor
#'     summaries, without shrinkage. The effect is added to a copy of the
#'     model (returned as \code{model}); old effects are untouched.}
#'   \item{New batch without a bridging reference}{Abstain.}
#'   \item{Known batch or population without an effect}{Abstain. This
#'     covers batches dropped by \code{unbridged = "abstain"}.}
#' }
#' Abstain keeps the original value exactly.
#'
#' @param model A model from \code{\link{bn_fit}} or \code{\link{bn_fit_paired}}.
#' @param x Data in the same form as for \code{\link{bn_fit}}.
#' @param samples Sample table covering every sample in \code{x}.
#' @param population Per-cell population labels for \code{x}.
#' @param transform The scale \code{x} is on; it must equal the model's.
#'   Defaults to the model's own transform.
#' @param unsupported \code{"abstain"} keeps uncorrectable values and
#'   reports them; \code{"error"} stops when any action is an abstention.
#' @param sample_ids Optional sample id per row; required for a matrix
#'   \code{x}.
#'
#' @return A \code{bn_result}: a list with \code{x} (corrected data, same
#'   shape, row order and columns as the input), \code{actions} (one row per
#'   batch, population and marker: \code{action}, \code{reason_code},
#'   \code{reason}, estimate and support columns, \code{n_cells}),
#'   \code{model} (the input model, or a copy extended with new-batch
#'   effects, with a new digest and \code{parent_digest}), \code{excluded}
#'   and \code{aliquot_source}.
#' @example inst/examples/bn-toy.R
#' @examples
#' model <- bn_fit(
#'   x, samples, population = rep("T", nrow(x)), markers = c("Lin", "M"),
#'   transform = list(name = "asinh", cofactor = 5), anchor = "b1",
#'   grid = c(2, 1), mapping_markers = "Lin", seed = 1, rlen = 20
#' )
#' res <- bn_apply(model, x, samples, rep("T", nrow(x)))
#' res
#' # Mean of M in the study samples per batch, before and after.
#' tapply(x$M[startsWith(x$sample_id, "s")], x$sample_id[startsWith(x$sample_id, "s")], mean)
#' tapply(res$x$M[startsWith(x$sample_id, "s")], x$sample_id[startsWith(x$sample_id, "s")], mean)
#' @seealso \code{\link{bn_fit}}
#' @export
bn_apply <- function(
  model,
  x,
  samples,
  population,
  transform = model$transform,
  unsupported = "abstain",
  sample_ids = NULL
) {
  res <- .bn_apply(
    model,
    x,
    samples,
    population,
    transform,
    unsupported,
    sample_ids
  )
  # Say once which batches were left entirely raw: actions and print() alone are easy to miss.
  raw <- tapply(res$actions$action == "abstain", res$actions$batch, all)
  raw <- names(raw)[raw]
  if (length(raw) > 0) {
    message(
      "bn_apply(): no correction applied to batch(es) ",
      paste(raw, collapse = ", "),
      "; their values are unchanged. See result$actions$reason."
    )
  }
  res
}

# bn_apply() without the user-facing message; the fit passes and component models call this.
.bn_apply <- function(
  model,
  x,
  samples,
  population,
  transform = model$transform,
  unsupported = "abstain",
  sample_ids = NULL
) {
  .bn_stop_missing(c(
    model = missing(model),
    x = missing(x),
    samples = missing(samples),
    population = missing(population)
  ))
  .bn_check_choice(
    unsupported,
    "unsupported",
    c("abstain", "error"),
    what = "bn_apply"
  )
  # A two-pass quantile model (bn_fit) applies its node passes and quantile maps.
  if (inherits(model, "bn_node_quantile_model")) {
    return(.bn_apply_node_quantile(
      model,
      x,
      samples,
      population,
      transform,
      unsupported,
      sample_ids
    ))
  }
  # A node-level model (.bn_fit_nodes) assigns cells to nodes, then applies its component models.
  if (inherits(model, "bn_node_model")) {
    return(.bn_apply_nodes(
      model,
      x,
      samples,
      population,
      transform,
      unsupported,
      sample_ids
    ))
  }
  .bn_validate(model)
  transform <- .bn_check_transform(transform)
  if (!identical(transform, model$transform)) {
    stop(
      "Transform mismatch: the model was fitted on ",
      .bn_transform_label(model$transform),
      " but 'x' is declared as ",
      .bn_transform_label(transform),
      ".",
      call. = FALSE
    )
  }
  absent <- setdiff(model$markers, colnames(x))
  if (length(absent) > 0) {
    stop(
      "Marker set mismatch: model marker(s) not in 'x': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  d <- .bn_cells(x, samples, population, model$markers, sample_ids)

  # Only batches the model has never seen may get new effects.
  new_batches <- setdiff(
    sort(unique(d$cells$batch), method = "radix"),
    model$known_batches
  )
  new_eff <- .bn_empty_effects()
  new_exc <- .bn_empty_excluded()
  if (length(new_batches) > 0) {
    bridging <- d$cells$batch %in%
      new_batches &
      d$cells$is_reference &
      d$cells$reference_id %in% model$targets$reference_id
    s <- .bn_summaries(
      d$xm,
      d$cells[bridging, , drop = FALSE],
      model$min_cells
    )
    new_exc <- s$excluded
    if (nrow(s$summ) > 0) {
      new_eff <- .bn_paired_effects(
        .bn_donor_level(s$summ),
        model$targets,
        "integrate"
      )
    }
  }

  # Extend a copy; the input model and its effects stay as they were.
  updated <- model
  if (nrow(new_eff) > 0) {
    added <- sort(unique(new_eff$batch), method = "radix")
    updated$effects <- dplyr::bind_rows(model$effects, new_eff)
    updated$known_batches <- sort(
      c(model$known_batches, added),
      method = "radix"
    )
    updated$excluded <- dplyr::bind_rows(
      model$excluded,
      new_exc[new_exc$batch %in% added, , drop = FALSE]
    )
    updated$parent_digest <- model$digest
    updated$digest <- .bn_digest(updated)
  }

  applied <- .bn_apply_effects(
    updated,
    d$xm,
    d$cells,
    new_batches = new_batches,
    estimated_now = unique(new_eff$batch)
  )
  if (unsupported == "error" && any(applied$actions$action == "abstain")) {
    ab <- applied$actions[applied$actions$action == "abstain", ]
    stop(
      sprintf(
        "%d (batch, population, marker) combination(s) cannot be corrected, e.g. %s / %s / %s: %s",
        nrow(ab),
        ab$batch[1],
        ab$population[1],
        ab$marker[1],
        ab$reason[1]
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      x = .bn_restore(x, applied$xm, model$markers),
      actions = applied$actions,
      model = updated,
      excluded = new_exc,
      aliquot_source = d$aliquot_source
    ),
    class = "bn_result"
  )
}

# Write corrected marker columns back into the input's shape.
.bn_restore <- function(x, xm, markers) {
  out <- x
  if (is.matrix(x)) {
    out[, markers] <- xm
  } else {
    for (m in markers) {
      out[[m]] <- unname(xm[, m])
    }
  }
  out
}

# Subtract each (batch, population, marker) effect; abstain where none exists.
.bn_apply_effects <- function(model, xm, cells, new_batches, estimated_now) {
  markers <- model$markers
  groups <- unique(cells[, c("batch", "population")])
  groups <- groups[order(groups$batch, groups$population), , drop = FALSE]
  group_key <- paste(groups$batch, groups$population, sep = "\r")
  n_cells <- as.integer(table(factor(
    paste(cells$batch, cells$population, sep = "\r"),
    levels = group_key
  )))

  # One action row per (batch, population, marker) present in x.
  idx <- rep(seq_len(nrow(groups)), each = length(markers))
  actions <- tibble::tibble(
    batch = groups$batch[idx],
    population = groups$population[idx],
    marker = rep(markers, nrow(groups)),
    n_cells = n_cells[idx]
  )
  eff <- model$effects[, .bn_effect_cols]
  actions <- dplyr::left_join(
    actions,
    eff,
    by = c("batch", "population", "marker"),
    relationship = "many-to-one"
  )

  has_eff <- !is.na(actions$estimate)
  pop_unknown <- !actions$population %in% model$population_levels
  is_new <- actions$batch %in% new_batches
  now <- actions$batch %in% estimated_now
  unbridged <- actions$batch %in% model$bridges$batch[!model$bridges$bridged]
  actions$action <- ifelse(has_eff & !pop_unknown, "apply", "abstain")
  actions$reason_code <- dplyr::case_when(
    pop_unknown ~ "population_unknown",
    has_eff & now ~ "estimated_new_batch",
    has_eff ~ "frozen_effect",
    is_new ~ "no_bridging_reference",
    unbridged ~ "unbridged_at_fit",
    .default = "known_batch_no_effect"
  )
  actions$reason <- dplyr::case_when(
    pop_unknown ~ "Population label is not in the frozen model.",
    has_eff & now ~ sprintf(
      "New batch: effect estimated from %d reference donor(s) in this batch against frozen anchor summaries.",
      actions$n_donors
    ),
    has_eff ~ sprintf(
      "Frozen effect from %d reference donor(s) (%s).",
      actions$n_donors,
      actions$support
    ),
    is_new ~ "New batch has no reference shared with the anchor (after min_cells); original values kept.",
    unbridged ~ "Batch shared no reference with the anchor at fit time (unbridged = \"abstain\"); original values kept.",
    .default = "Batch is in the model but has no effect for this population; original values kept. Refit to change this."
  )
  actions <- actions[, c(
    "batch",
    "population",
    "marker",
    "action",
    "reason_code",
    "reason",
    "estimate",
    "se_donor",
    "n_donors",
    "n_aliquots",
    "n_aliquots_anchor",
    "support",
    "n_cells"
  )]

  # Subtract per cell; cells under an abstain action keep their exact value.
  out <- xm
  act_key <- paste(
    actions$batch,
    actions$population,
    actions$marker,
    sep = "\r"
  )
  shift <- ifelse(actions$action == "apply", actions$estimate, NA_real_)
  for (m in markers) {
    i <- match(paste(cells$batch, cells$population, m, sep = "\r"), act_key)
    s <- shift[i]
    ok <- !is.na(s)
    out[ok, m] <- xm[ok, m] - s[ok]
  }
  list(xm = out, actions = actions)
}

#' @export
print.bn_result <- function(x, ...) {
  n <- table(factor(x$actions$action, levels = c("apply", "abstain")))
  cat(
    "<bn_result>",
    n[["apply"]],
    "applied,",
    n[["abstain"]],
    "abstained",
    "(batch x population x marker) | model digest:",
    x$model$digest,
    "\n"
  )
  invisible(x)
}
