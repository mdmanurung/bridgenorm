# Two-pass node correction with per-node quantile maps: bn_fit() and its apply path.
#

.bn_node_quantile_schema_version <- "bn_node_quantile_model/1"

#' Fit a two-pass node model with per-node quantile maps
#'
#' Batch effects can change the shape of a marker's distribution within a
#' cell state, not only its location. \code{bn_fit()} runs a
#' node shift (\code{.bn_fit_nodes()}), then refits a second map on the
#' shifted cells and, in each second-pass node, maps every batch's marker
#' quantiles onto the anchor's.
#'
#' @details
#' Steps, all on reference cells only:
#' \enumerate{
#'   \item Pass 1: \code{.bn_fit_nodes()} with the given settings; its
#'     correction is applied to the reference cells.
#'   \item Pass 2: \code{.bn_fit_nodes()} again, with the same settings,
#'     on the pass-1 cells. Its node shifts cover the cells pass 2 cannot
#'     map by quantiles.
#'   \item For each (batch, pass-2 node) with at least \code{min_cells}
#'     reference cells pooled over its reference samples, in a batch whose
#'     reference donors (\code{reference_id}) are exactly the anchor's, the
#'     quantiles at
#'     \code{quantile_probs} of every marker (\code{stats::quantile()},
#'     type 7) are stored.
#'   \item Shrinkage: in each such (batch, node) of a non-anchor batch whose
#'     node also has anchor quantiles, the offset curve (batch minus anchor
#'     quantiles) is pulled toward the batch-wide offset curve (quantiles of
#'     all the batch's pass-1 reference cells minus the anchor's), with weight
#'     \eqn{n / (n + shrink_k)} on the node's own curve, \eqn{n} being the
#'     batch's reference cells in the node. Each column of a shrunk table is
#'     then sorted (rearrangement; Chernozhukov, Fernandez-Val and Galichon
#'     2010), because the mix can cross itself. \code{shrink_k = 0} keeps the
#'     fitted tables.
#' }
#' At integrate time, pass 1 is applied, then each cell is assigned to its
#' nearest pass-2 node. If its batch and the anchor both have quantiles in
#' that node, each marker value is mapped by linear interpolation from the
#' batch's quantiles to the anchor's; beyond the outer quantiles it moves by
#' the difference of the outer quantiles (\code{tail = "shift"}). Other cells
#' take the pass-2 node shift. Anchor cells are mapped onto their own
#' quantiles, which leaves values between the outer quantiles unchanged.
#'
#' A batch with other reference donors than the anchor has no quantiles:
#' pooled distributions of different donors differ by biology, not only by
#' batch. Such a batch takes the pass-2 node shift if it is bridged, and is
#' left raw if it is not.
#'
#' The quantile rule pools wells, while the node shifts count cells per
#' reference sample; a (batch, node) can therefore have quantiles without a
#' node shift.
#'
#' @param x A numeric matrix (cells by markers) or a data frame with one row
#'   per cell, already transformed. A data frame needs a \code{sample_id}
#'   column unless \code{sample_ids} is given.
#' @param samples Data frame, one row per sample: \code{sample_id},
#'   \code{batch}, \code{reference_id} (the bridge donor; \code{NA} for study
#'   samples) and optionally \code{aliquot}.
#' @param population Per-cell population labels, one per row of \code{x},
#'   frozen before correction.
#' @param markers Marker columns to correct.
#' @param transform \code{list(name = "asinh", cofactor = <number>)}: the
#'   scale \code{x} is on. It is recorded, never applied.
#' @param anchor The batch that defines the target scale.
#' @param grid SOM grid as \code{c(xdim, ydim)}.
#' @param mapping_markers Markers that define cell state for the map; each
#'   must be one of \code{markers}.
#' @param seed Seed for the training-cell draw and the SOM.
#' @param min_cells Minimum cells per reference aliquot summary and per
#'   (batch, node) quantile table.
#' @param unbridged \code{"abstain"} leaves a batch that shares no bridge
#'   donor with the anchor raw and flags it; \code{"stop"} refuses the fit.
#' @param method Node estimator; only \code{"paired"}.
#' @param eb_settings Only \code{NULL}; kept for future estimators.
#' @param train_cells Maximum training cells per reference sample for the SOM.
#' @param rlen Number of SOM training passes (\code{kohonen::som}).
#' @param assignment How cells are assigned to nodes; only \code{"nearest"}.
#' @param fallback For a (batch, node, marker) without a node effect:
#'   \code{"global"} takes the batch's global shift, \code{"abstain"} keeps
#'   the raw value.
#' @param quantile_probs Increasing probabilities in (0, 1), at least two,
#'   at which batch and anchor quantiles are matched.
#' @param tail Rule beyond the outer quantiles; only \code{"shift"}.
#' @param shrink_k Non-negative pseudo-count of cells for the shrinkage toward
#'   the batch-wide offset curve. Larger values borrow more; \code{0} turns
#'   shrinkage off. \code{"eb"} estimates one pseudo-count per marker by
#'   empirical Bayes from how much node offsets vary beyond their sampling
#'   noise, never below 300; the values are stored in
#'   \code{settings$shrink_k_marker}. It is kept for study; in benchmarks
#'   it did not beat the scalar default.
#' @param target \code{"anchor"} corrects every batch to the anchor batch.
#'   \code{"virtual"} corrects every batch, the anchor included, to a virtual
#'   reference: per reference summary and per SOM node, the equal-weight mean
#'   over batches (for quantiles, the mean of the batches' quantile curves).
#'   It covers the same summaries and nodes as the anchor, which still sets
#'   the bridges and the map. Not combined with \code{shrink_k = "eb"}.
#' @param sample_ids Optional sample id per row; required for a matrix
#'   \code{x}.
#'
#' @section Defaults:
#' The defaults come from benchmarking: \code{min_cells = 50},
#' \code{train_cells = 500}, \code{rlen = 100}, 21 quantiles from 0.025 to
#' 0.975, \code{shrink_k = 300}. Check them on new data. \code{grid},
#' \code{mapping_markers}, \code{anchor} and \code{seed} have no default.
#'
#' @return A \code{bn_node_quantile_model}: a list with
#'   \code{schema_version}, \code{settings}, \code{markers},
#'   \code{transform}, \code{anchor}, \code{pass1} and \code{pass2} (frozen
#'   node models), \code{quantiles} (one probabilities-by-markers matrix per
#'   "batch\\rnode"), \code{software}, \code{created},
#'   \code{parent_digest} and \code{digest}.
#' @seealso \code{\link{bn_apply}}, \code{\link{bn_fit_paired}}
#' @export
bn_fit <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  grid,
  mapping_markers,
  seed,
  min_cells = 50,
  unbridged = "abstain",
  method = "paired",
  eb_settings = NULL,
  train_cells = 500,
  rlen = 100,
  assignment = "nearest",
  fallback = "global",
  quantile_probs = seq(0.025, 0.975, length.out = 21),
  tail = "shift",
  shrink_k = 300,
  target = "anchor",
  sample_ids = NULL
) {
  .bn_stop_missing(c(
    x = missing(x),
    samples = missing(samples),
    population = missing(population),
    markers = missing(markers),
    transform = missing(transform),
    anchor = missing(anchor),
    grid = missing(grid),
    mapping_markers = missing(mapping_markers),
    seed = missing(seed)
  ))
  .bn_check_choice(tail, "tail", "shift", what = "bn_fit")
  if (
    !is.numeric(quantile_probs) ||
      length(quantile_probs) < 2 ||
      any(!is.finite(quantile_probs)) ||
      any(quantile_probs <= 0 | quantile_probs >= 1) ||
      any(diff(quantile_probs) <= 0)
  ) {
    stop(
      "'quantile_probs' must be at least two increasing probabilities in (0, 1).",
      call. = FALSE
    )
  }
  if (
    !identical(shrink_k, "eb") &&
      (!is.numeric(shrink_k) ||
        length(shrink_k) != 1 ||
        !is.finite(shrink_k) ||
        shrink_k < 0)
  ) {
    stop("'shrink_k' must be one finite number >= 0, or \"eb\".", call. = FALSE)
  }
  .bn_check_target(target)
  if (target == "virtual" && identical(shrink_k, "eb")) {
    stop("target = \"virtual\" does not take shrink_k = \"eb\".", call. = FALSE)
  }
  if (target == "virtual" && .bn_virtual %in% samples$batch) {
    stop(
      sprintf(
        "A batch is named '%s', the virtual reference's key.",
        .bn_virtual
      ),
      call. = FALSE
    )
  }

  # Reference rows only, as .bn_fit_nodes() uses.
  d <- .bn_cells(x, samples, population, markers, sample_ids)
  ref_rows <- which(d$cells$is_reference)
  x_ref <- x[ref_rows, , drop = FALSE]
  ids_ref <- d$cells$sample_id[ref_rows]
  pop_ref <- population[ref_rows]
  fit_pass <- function(xp) {
    .bn_fit_nodes(
      x = xp,
      samples = samples,
      population = pop_ref,
      markers = markers,
      transform = transform,
      anchor = anchor,
      min_cells = min_cells,
      unbridged = unbridged,
      method = method,
      eb_settings = eb_settings,
      grid = grid,
      mapping_markers = mapping_markers,
      train_cells = train_cells,
      rlen = rlen,
      seed = seed,
      assignment = assignment,
      fallback = fallback,
      sample_ids = ids_ref,
      target = target
    )
  }
  integrate_pass <- function(m, xp) {
    .bn_apply(
      m,
      xp,
      samples,
      pop_ref,
      transform = transform,
      unsupported = "abstain",
      sample_ids = ids_ref
    )
  }

  # Pass 1 node shift, then pass 2 on the shifted reference cells.
  pass1 <- fit_pass(x_ref)
  x1 <- integrate_pass(pass1, x_ref)$x
  pass2 <- fit_pass(x1)
  nodes <- integrate_pass(pass2, x1)$nodes

  # Quantiles per (batch, pass-2 node) with enough pooled reference cells.
  x1m <- .bn_marker_matrix(x1, markers)
  groups <- split(
    seq_along(nodes),
    paste(d$cells$batch[ref_rows], nodes, sep = "\r")
  )
  n_cells <- lengths(groups)
  # A quantile map matches pooled reference distributions, so only batches whose reference donors are exactly
  # the anchor's get tables; an unbridged or partly bridged batch keeps the pass-2 shift (heb 0.7.2).
  # Donors count only through aliquots with at least min_cells cells, the gate the shift passes use.
  ref_cells <- d$cells[ref_rows, ]
  aliquot_key <- paste(ref_cells$batch, ref_cells$aliquot, sep = "\r")
  big <- aliquot_key %in% names(which(table(aliquot_key) >= min_cells))
  donors <- lapply(
    split(
      ref_cells$reference_id[big],
      factor(ref_cells$batch[big], levels = unique(ref_cells$batch))
    ),
    \(r) {
      sort(unique(r), method = "radix")
    }
  )
  same_donors <- names(donors)[vapply(
    donors,
    identical,
    logical(1),
    donors[[anchor]]
  )]
  groups <- groups[
    n_cells >= min_cells & sub("\r.*", "", names(groups)) %in% same_donors
  ]
  no_tables <- setdiff(unique(d$cells$batch), c(same_donors, anchor))
  if (length(no_tables) > 0) {
    warning(
      sprintf(
        "No quantile maps for batch(es) %s: their reference donors passing min_cells differ from the anchor's (%s). They get node shifts only, or stay raw if unbridged; see the fit's actions.",
        paste(no_tables, collapse = ", "),
        paste(donors[[anchor]], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  # Pooled tables weight donors by their cell counts; unequal shares between a batch and the anchor bias the
  # map (documented limit).
  if (length(donors[[anchor]]) > 1) {
    share <- prop.table(
      table(ref_cells$batch[big], ref_cells$reference_id[big]),
      1
    )
    off <- same_donors[vapply(
      same_donors,
      \(b) max(abs(share[b, ] - share[anchor, ])) > 0.1,
      logical(1)
    )]
    if (length(off) > 0) {
      warning(
        sprintf(
          "Bridge donors' cell shares differ from the anchor's by more than 10 percentage points in batch(es) %s. Quantile maps pool donors, so unequal shares bias the correction.",
          paste(off, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  q_table <- function(i) {
    q <- apply(
      x1m[i, , drop = FALSE],
      2,
      stats::quantile,
      probs = quantile_probs,
      names = FALSE,
      type = 7
    )
    matrix(q, nrow = length(quantile_probs), dimnames = list(NULL, markers))
  }
  quantiles <- lapply(groups, q_table)

  # Target curves per node and batch-wide: the anchor's, or with target = "virtual" the equal-weight mean of the
  # batches' curves (VR3) in each node the anchor has a table for.
  ref_batch <- d$cells$batch[ref_rows]
  batch_q <- lapply(split(seq_along(ref_batch), ref_batch), q_table)
  node_of <- \(g) sub(".*\r", "", g)
  covered <- Filter(
    \(g) !is.null(quantiles[[paste(anchor, node_of(g), sep = "\r")]]),
    names(quantiles)
  )
  if (target == "virtual") {
    by_node <- split(covered, node_of(covered))
    target_q <- lapply(by_node, \(gs) Reduce(`+`, quantiles[gs]) / length(gs))
    names(target_q) <- paste(.bn_virtual, names(by_node), sep = "\r")
    target_batch_q <- Reduce(`+`, batch_q[same_donors]) / length(same_donors)
    target_key <- \(g) paste(.bn_virtual, node_of(g), sep = "\r")
    todo <- covered
  } else {
    target_q <- quantiles
    target_batch_q <- batch_q[[anchor]]
    target_key <- \(g) paste(anchor, node_of(g), sep = "\r")
    todo <- covered[sub("\r.*", "", covered) != anchor]
  }

  # Shrink each node offset curve toward its batch-wide offset curve (TR3k); the anchor's own curves stay as
  # fitted when it is the target.
  k_marker <- if (identical(shrink_k, "eb")) {
    .bn_eb_shrink_k(
      quantiles,
      groups,
      x1m,
      todo,
      anchor,
      batch_q,
      quantile_probs
    )
  } else {
    stats::setNames(rep(shrink_k, length(markers)), markers)
  }
  if (any(k_marker > 0)) {
    shrunk <- quantiles
    for (g in todo) {
      b <- sub("\r.*", "", g)
      qa <- target_q[[target_key(g)]]
      # One weight per marker column; with a number shrink_k every column has the same weight.
      w <- matrix(
        n_cells[[g]] / (n_cells[[g]] + k_marker),
        nrow(qa),
        length(markers),
        byrow = TRUE
      )
      q <- qa +
        w * (quantiles[[g]] - qa) +
        (1 - w) * (batch_q[[b]] - target_batch_q)
      # Sort each column (rearrangement): the mix of three monotone curves can cross itself, and a
      # non-monotone table would swap the order of cells (heb 0.7.1).
      shrunk[[g]] <- apply(q, 2, sort)
    }
    quantiles <- shrunk
  }
  if (target == "virtual") {
    quantiles <- c(quantiles, target_q)
  }

  model <- structure(
    list(
      schema_version = .bn_node_quantile_schema_version,
      settings = list(
        quantile_probs = quantile_probs,
        tail = tail,
        min_cells = min_cells,
        quantile_type = 7L,
        shrink_k = shrink_k,
        shrink_k_marker = k_marker,
        target = target
      ),
      markers = markers,
      transform = pass1$transform,
      anchor = anchor,
      pass1 = pass1,
      pass2 = pass2,
      quantiles = quantiles,
      software = list(R = R.version.string, bridgenorm = .bn_pkg_version()),
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      parent_digest = NA_character_,
      digest = NA_character_
    ),
    class = "bn_node_quantile_model"
  )
  model$digest <- .bn_digest(model)
  model
}

# Marker columns of x as a double matrix.
.bn_marker_matrix <- function(x, markers) {
  out <- if (is.matrix(x)) {
    x[, markers, drop = FALSE]
  } else {
    as.matrix(as.data.frame(x)[markers])
  }
  storage.mode(out) <- "double"
  out
}

# Stop unless `model` is an intact bn_node_quantile_model with intact passes.
.bn_validate_node_quantile <- function(model) {
  if (!identical(model$schema_version, .bn_node_quantile_schema_version)) {
    stop(
      sprintf(
        "Unsupported schema_version '%s'; this package reads '%s'.",
        as.character(model$schema_version),
        .bn_node_quantile_schema_version
      ),
      call. = FALSE
    )
  }
  if (!identical(model$digest, .bn_digest(model))) {
    stop(
      "Digest mismatch: the model changed after it was frozen. Refit, or load the original file.",
      call. = FALSE
    )
  }
  .bn_validate_node(model$pass1)
  .bn_validate_node(model$pass2)
  invisible(model)
}

# bn_apply() for a bn_node_quantile_model: pass 1, then per-node quantile maps or the pass-2 shift.
.bn_apply_node_quantile <- function(
  model,
  x,
  samples,
  population,
  transform,
  unsupported,
  sample_ids
) {
  .bn_validate_node_quantile(model)
  d <- .bn_cells(x, samples, population, model$markers, sample_ids)
  ids <- d$cells$sample_id
  r1 <- .bn_apply(
    model$pass1,
    x,
    samples,
    population,
    transform,
    "abstain",
    sample_ids = ids
  )
  r2 <- .bn_apply(
    model$pass2,
    r1$x,
    samples,
    population,
    transform,
    "abstain",
    sample_ids = ids
  )
  v <- .bn_marker_matrix(r1$x, model$markers)
  out <- .bn_marker_matrix(r2$x, model$markers)

  # Quantile map where the cell's batch and the target both have quantiles in its pass-2 node. Models from
  # 0.9.0 have no settings$target and map to the anchor.
  target <- if (identical(model$settings$target, "virtual")) {
    .bn_virtual
  } else {
    model$anchor
  }
  key <- paste(d$cells$batch, r2$nodes, sep = "\r")
  mapped_keys <- character(0)
  tied <- list()
  for (g in intersect(unique(key), names(model$quantiles))) {
    qa <- model$quantiles[[paste(target, sub(".*\r", "", g), sep = "\r")]]
    if (is.null(qa)) {
      next
    }
    qb <- model$quantiles[[g]]
    i <- which(key == g)
    for (m in model$markers) {
      # A marker with one value in the batch's quantile table has no map; it keeps the pass-2 node shift.
      if (length(unique(qb[, m])) < 2) {
        tied[[g]] <- c(tied[[g]], m)
        next
      }
      vi <- v[i, m]
      res <- stats::approx(qb[, m], qa[, m], xout = vi, ties = mean, rule = 2)$y
      low <- vi < qb[1, m]
      high <- vi > qb[nrow(qb), m]
      res[low] <- vi[low] + qa[1, m] - qb[1, m]
      res[high] <- vi[high] + qa[nrow(qa), m] - qb[nrow(qb), m]
      out[i, m] <- res
    }
    mapped_keys <- c(mapped_keys, g)
  }

  # Actions: pass 1 and pass 2 rows; pass-2 rows of quantile-mapped (batch, node) read node_quantile.
  a1 <- r1$actions
  a1$pass <- 1L
  a2 <- r2$actions
  a2$pass <- 2L
  a2_key <- paste(a2$batch, a2$population, sep = "\r")
  # Tied rows (one quantile value) keep the pass-2 node shift and its action.
  t_row <- vapply(
    seq_along(a2_key),
    function(r) a2$marker[[r]] %in% tied[[a2_key[[r]]]],
    logical(1)
  )
  q_row <- a2_key %in% mapped_keys & !t_row
  a2$action[q_row] <- "apply"
  a2$reason_code[q_row] <- "node_quantile"
  a2$reason[
    q_row
  ] <- "Batch and anchor quantiles in this node; the quantile map was applied."
  a2$reason_code[t_row] <- "node_quantile_tied"
  a2$reason[
    t_row
  ] <- "One value in this batch's quantiles for this node and marker; the node shift was kept."
  actions <- rbind(a1, a2)

  # Keep the pass models as bn_apply() returned them: they carry effects of new batches.
  updated <- model
  updated$pass1 <- r1$model
  updated$pass2 <- r2$model
  if (!identical(updated, model)) {
    updated$parent_digest <- model$digest
    updated$digest <- .bn_digest(updated)
  }

  if (unsupported == "error" && any(actions$action == "abstain")) {
    ab <- actions[actions$action == "abstain", ]
    stop(
      sprintf(
        "%d (pass, batch, node, marker) row(s) cannot be corrected, e.g. pass %d / %s / %s / %s: %s",
        nrow(ab),
        ab$pass[1],
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
      x = .bn_restore(x, out, model$markers),
      actions = actions,
      model = updated,
      nodes = r2$nodes,
      excluded = r2$excluded,
      aliquot_source = d$aliquot_source
    ),
    class = "bn_result"
  )
}

#' @export
print.bn_node_quantile_model <- function(x, ...) {
  cat(
    "<bn_node_quantile_model>",
    "| anchor:",
    x$anchor,
    "| nodes per pass:",
    nrow(x$pass1$map$nodes),
    "| quantile tables:",
    length(x$quantiles),
    "| digest:",
    x$digest,
    "\n"
  )
  invisible(x)
}

# Key of the virtual reference's quantile tables in model$quantiles ("(virtual)\rnode"); no batch may take it.
.bn_virtual <- "(virtual)"

# Empirical-Bayes shrink_k per marker (shrink_k = "eb"). The median offset
# of a (batch, node) from its batch-wide target, d, has variance tau^2 (true node deviation) plus sampling noise
# v = (pi/2) (s2_g / n_g + s2_a / n_a). tau^2 = mean(d^2) - mean(v), floored at 1e-6; k = (pi/2) median(s2_g) / tau^2,
# the normal-means shrinkage weight n / (n + k). k never goes below 300, the tested floor.
.bn_eb_shrink_k <- function(
  quantiles,
  groups,
  x1m,
  todo,
  anchor,
  batch_q,
  quantile_probs,
  floor = 300
) {
  mid <- which.min(abs(quantile_probs - 0.5))
  s2 <- t(vapply(
    groups,
    \(i) apply(x1m[i, , drop = FALSE], 2, stats::var),
    numeric(ncol(x1m))
  ))
  n <- lengths(groups)
  parts <- lapply(todo, \(g) {
    a <- paste(anchor, sub(".*\r", "", g), sep = "\r")
    target <- batch_q[[sub("\r.*", "", g)]] - batch_q[[anchor]]
    d <- quantiles[[g]][mid, ] - quantiles[[a]][mid, ] - target[mid, ]
    v <- (pi / 2) * (s2[g, ] / n[[g]] + s2[a, ] / n[[a]])
    rbind(d2 = d^2, v = v)
  })
  tau2 <- pmax(
    Reduce(`+`, lapply(parts, \(p) p["d2", ])) /
      length(todo) -
      Reduce(`+`, lapply(parts, \(p) p["v", ])) / length(todo),
    1e-6
  )
  k <- (pi / 2) * apply(s2[todo, , drop = FALSE], 2, stats::median) / tau2
  stats::setNames(pmax(k, floor), colnames(x1m))
}
