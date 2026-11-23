# Two-pass node correction with per-node quantile maps: bn_fit_node_quantile() and its integrate path.
#
# Ported from the tusoskill run bn_nodes_20261121, candidate c016 (decisions T1-T3).

.bn_node_quantile_schema_version <- "bn_node_quantile_model/1"

#' Fit a two-pass node model with per-node quantile maps
#'
#' Batch effects can change the shape of a marker's distribution within a
#' cell state, not only its location. \code{bn_fit_node_quantile()} runs a
#' node shift (\code{\link{bn_fit_nodes}}), then refits a second map on the
#' shifted cells and, in each second-pass node, maps every batch's marker
#' quantiles onto the anchor's.
#'
#' @details
#' Steps, all on reference cells only:
#' \enumerate{
#'   \item Pass 1: \code{\link{bn_fit_nodes}} with the given settings; its
#'     correction is applied to the reference cells.
#'   \item Pass 2: \code{\link{bn_fit_nodes}} again, with the same settings,
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
#' @inheritParams bn_fit_nodes
#' @param quantile_probs Increasing probabilities in (0, 1), at least two,
#'   at which batch and anchor quantiles are matched.
#' @param tail Rule beyond the outer quantiles; only \code{"shift"}.
#' @param shrink_k Non-negative number of cells: the pseudo-count of the
#'   shrinkage toward the batch-wide offset curve. Larger values borrow more;
#'   \code{0} turns shrinkage off. No default: it changes the estimates
#'   (decision TR3k used 300).
#' @return A \code{bn_node_quantile_model}: a list with
#'   \code{schema_version}, \code{settings}, \code{markers},
#'   \code{transform}, \code{anchor}, \code{pass1} and \code{pass2} (frozen
#'   \code{bn_node_model}s), \code{quantiles} (one probabilities-by-markers
#'   matrix per "batch\\rnode"), \code{software}, \code{created},
#'   \code{parent_digest} and \code{digest}.
#' @seealso \code{\link{bn_fit_nodes}}, \code{\link{bn_integrate}}
#' @export
bn_fit_node_quantile <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  min_cells,
  unbridged,
  method,
  eb_settings,
  grid,
  mapping_markers,
  train_cells,
  rlen,
  seed,
  assignment,
  fallback,
  quantile_probs,
  tail,
  shrink_k,
  sample_ids = NULL
) {
  .bn_stop_missing(c(
    x = missing(x),
    samples = missing(samples),
    population = missing(population),
    markers = missing(markers),
    transform = missing(transform),
    anchor = missing(anchor),
    min_cells = missing(min_cells),
    unbridged = missing(unbridged),
    method = missing(method),
    eb_settings = missing(eb_settings),
    grid = missing(grid),
    mapping_markers = missing(mapping_markers),
    train_cells = missing(train_cells),
    rlen = missing(rlen),
    seed = missing(seed),
    assignment = missing(assignment),
    fallback = missing(fallback),
    quantile_probs = missing(quantile_probs),
    tail = missing(tail),
    shrink_k = missing(shrink_k)
  ))
  .bn_check_choice(tail, "tail", "shift", what = "bn_fit_node_quantile")
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
    !is.numeric(shrink_k) ||
      length(shrink_k) != 1 ||
      !is.finite(shrink_k) ||
      shrink_k < 0
  ) {
    stop("'shrink_k' must be one finite number >= 0.", call. = FALSE)
  }

  # Reference rows only, as bn_fit_nodes() uses.
  d <- .bn_cells(x, samples, population, markers, sample_ids)
  ref_rows <- which(d$cells$is_reference)
  x_ref <- x[ref_rows, , drop = FALSE]
  ids_ref <- d$cells$sample_id[ref_rows]
  pop_ref <- population[ref_rows]
  fit_pass <- function(xp) {
    bn_fit_nodes(
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
      sample_ids = ids_ref
    )
  }
  integrate_pass <- function(m, xp) {
    bn_integrate(
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
  ref_cells <- d$cells[ref_rows, ]
  donors <- lapply(split(ref_cells$reference_id, ref_cells$batch), \(r) sort(unique(r)))
  same_donors <- names(donors)[vapply(donors, identical, logical(1), donors[[anchor]])]
  groups <- groups[n_cells >= min_cells & sub("\r.*", "", names(groups)) %in% same_donors]
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

  # Shrink each non-anchor node offset curve toward its batch-wide offset curve (TR3k).
  if (shrink_k > 0) {
    ref_batch <- d$cells$batch[ref_rows]
    batch_q <- lapply(split(seq_along(ref_batch), ref_batch), q_table)
    shrunk <- quantiles
    for (g in names(quantiles)) {
      b <- sub("\r.*", "", g)
      qa <- quantiles[[paste(anchor, sub(".*\r", "", g), sep = "\r")]]
      if (b == anchor || is.null(qa)) {
        next
      }
      w <- n_cells[[g]] / (n_cells[[g]] + shrink_k)
      q <- qa +
        w * (quantiles[[g]] - qa) +
        (1 - w) * (batch_q[[b]] - batch_q[[anchor]])
      # Sort each column (rearrangement): the mix of three monotone curves can cross itself, and a
      # non-monotone table would swap the order of cells (heb 0.7.1).
      shrunk[[g]] <- apply(q, 2, sort)
    }
    quantiles <- shrunk
  }

  model <- structure(
    list(
      schema_version = .bn_node_quantile_schema_version,
      settings = list(
        quantile_probs = quantile_probs,
        tail = tail,
        min_cells = min_cells,
        quantile_type = 7L,
        shrink_k = shrink_k
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

# bn_integrate() for a bn_node_quantile_model: pass 1, then per-node quantile maps or the pass-2 shift.
.bn_integrate_node_quantile <- function(
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
  r1 <- bn_integrate(
    model$pass1,
    x,
    samples,
    population,
    transform,
    "abstain",
    sample_ids = ids
  )
  r2 <- bn_integrate(
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

  # Quantile map where the cell's batch and the anchor both have quantiles in its pass-2 node.
  key <- paste(d$cells$batch, r2$nodes, sep = "\r")
  mapped_keys <- character(0)
  for (g in intersect(unique(key), names(model$quantiles))) {
    qa <- model$quantiles[[paste(model$anchor, sub(".*\r", "", g), sep = "\r")]]
    if (is.null(qa)) {
      next
    }
    qb <- model$quantiles[[g]]
    i <- which(key == g)
    for (m in model$markers) {
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
  q_row <- paste(a2$batch, a2$population, sep = "\r") %in% mapped_keys
  a2$action[q_row] <- "apply"
  a2$reason_code[q_row] <- "node_quantile"
  a2$reason[
    q_row
  ] <- "Batch and anchor quantiles in this node; the quantile map was applied."
  actions <- rbind(a1, a2)

  # Keep the pass models as bn_integrate() returned them: they carry effects of new batches.
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
