# Node-level correction: heb_fit_nodes() and its integrate path.
#
# Ported from experiments/data set A/27_som-grid.R (SOM grid sweep, decisions SG1-SG4, N1-N3).

.heb_node_schema_version <- "heb_node_model/1"

#' Fit a node-level batch-correction model
#'
#' Batch effects often differ between cell states. \code{heb_fit_nodes()}
#' learns a self-organising map (SOM) on the reference cells and fits the
#' chosen estimator (paired shift or \code{\link{heb_fit}}) separately in each
#' map node, with the node as the population. \code{\link{heb_integrate}} then
#' assigns every cell to its nearest node and subtracts that node's effect.
#'
#' @details
#' Steps, all on reference cells only:
#' \enumerate{
#'   \item A global model is fitted with the given \code{population} labels:
#'     paired shift always (for mapping), plus the chosen \code{method} (for
#'     the fallback).
#'   \item Mapping markers are pre-shifted by the global paired shift, so the
#'     map learns cell states rather than batches.
#'   \item Up to \code{train_cells} cells per reference sample, drawn with
#'     \code{seed}, train a \code{kohonen::som()} on the standardised mapping
#'     markers (\code{grid}, \code{rlen}).
#'   \item Every reference cell is assigned to its nearest node, and the chosen
#'     \code{method} is fitted with the node labels as populations.
#' }
#' At integrate time a (batch, node, marker) without an effect, for example a
#' node with fewer than \code{min_cells} reference cells in that batch, takes
#' the global effect (\code{fallback = "global"}) or keeps the raw value
#' (\code{"abstain"}). Both are reported in \code{actions}.
#'
#' The global pre-shift assumes the reference donor has the same cell-state
#' mix in every batch. If a batch's references miss a cell state, the global
#' shift absorbs that mix difference and can move cells to the wrong node.
#'
#' Every argument that changes an estimate is required, with no default.
#' Needs the \pkg{kohonen} package.
#'
#' @inheritParams heb_fit
#' @param method \code{"paired"} (\code{\link{heb_fit_paired}}) or
#'   \code{"heb"} (\code{\link{heb_fit}}).
#' @param heb_settings For \code{method = "heb"}: a named list with
#'   \code{hierarchy} and every model setting of \code{\link{heb_fit}}
#'   (\code{variance_floor} to \code{batch_shrinkage_tau2}). \code{NULL} for
#'   \code{"paired"}.
#' @param unbridged \code{"stop"} or \code{"abstain"}, as in
#'   \code{\link{heb_fit}}.
#' @param grid SOM grid as \code{c(xdim, ydim)}.
#' @param mapping_markers Markers that define cell state for the map.
#' @param train_cells Maximum training cells per reference sample.
#' @param rlen Number of SOM training passes (\code{kohonen::som}).
#' @param seed Seed for the training-cell draw and the SOM.
#' @param assignment How cells are assigned to nodes; only \code{"nearest"}.
#' @param fallback \code{"global"} or \code{"abstain"}; see Details.
#' @return A \code{heb_node_model}: a list with \code{schema_version},
#'   \code{method}, \code{settings}, \code{markers}, \code{transform},
#'   \code{anchor}, \code{map} (\code{centre}, \code{scale}, \code{nodes},
#'   \code{labels}, \code{mapping_markers}), \code{global_paired},
#'   \code{global} and \code{local} (frozen \code{heb_model}s),
#'   \code{software}, \code{created} and \code{digest}.
#' @seealso \code{\link{heb_fit}}, \code{\link{heb_fit_paired}},
#'   \code{\link{heb_integrate}}
#' @export
heb_fit_nodes <- function(
  x,
  samples,
  population,
  markers,
  transform,
  anchor,
  min_cells,
  unbridged,
  method,
  heb_settings,
  grid,
  mapping_markers,
  train_cells,
  rlen,
  seed,
  assignment,
  fallback,
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
    unbridged = missing(unbridged),
    method = missing(method),
    heb_settings = missing(heb_settings),
    grid = missing(grid),
    mapping_markers = missing(mapping_markers),
    train_cells = missing(train_cells),
    rlen = missing(rlen),
    seed = missing(seed),
    assignment = missing(assignment),
    fallback = missing(fallback)
  ))
  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop(
      "heb_fit_nodes() needs the 'kohonen' package: install.packages(\"kohonen\").",
      call. = FALSE
    )
  }
  .heb_check_choice(
    method,
    "method",
    c("paired", "heb"),
    what = "heb_fit_nodes"
  )
  .heb_check_choice(assignment, "assignment", "nearest", what = "heb_fit_nodes")
  .heb_check_choice(
    fallback,
    "fallback",
    c("global", "abstain"),
    what = "heb_fit_nodes"
  )
  if (method == "heb" && !is.list(heb_settings)) {
    stop(
      "method = \"heb\" needs 'heb_settings', a list of heb_fit() settings.",
      call. = FALSE
    )
  }
  if (method == "paired" && !is.null(heb_settings)) {
    stop("method = \"paired\" takes heb_settings = NULL.", call. = FALSE)
  }
  if (
    !is.numeric(grid) ||
      length(grid) != 2 ||
      any(grid < 1) ||
      any(grid != round(grid))
  ) {
    stop(
      "'grid' must be two positive whole numbers, c(xdim, ydim).",
      call. = FALSE
    )
  }
  if (prod(grid) < 2) {
    stop("'grid' must have at least two nodes.", call. = FALSE)
  }
  for (arg in c("train_cells", "rlen")) {
    v <- get(arg)
    if (!is.numeric(v) || length(v) != 1 || v < 1 || v != round(v)) {
      stop(
        sprintf("'%s' must be one positive whole number.", arg),
        call. = FALSE
      )
    }
  }
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
    stop("'seed' must be one number.", call. = FALSE)
  }
  absent <- setdiff(mapping_markers, colnames(x))
  if (length(absent) > 0) {
    stop(
      "Mapping marker(s) not in 'x': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Reference rows only: every fit and the map use them, as heb_fit() does.
  d <- .heb_cells(x, samples, population, markers, sample_ids)
  ref_rows <- which(d$cells$is_reference)
  x_ref <- x[ref_rows, , drop = FALSE]
  ids_ref <- d$cells$sample_id[ref_rows]
  pop_ref <- population[ref_rows]
  fit_one <- function(fit_method, pop, unbridged_rule) {
    args <- list(
      x = x_ref,
      samples = samples,
      population = pop,
      markers = markers,
      transform = transform,
      anchor = anchor,
      min_cells = min_cells,
      unbridged = unbridged_rule,
      sample_ids = ids_ref
    )
    if (fit_method == "paired") {
      do.call(heb_fit_paired, args)
    } else {
      do.call(heb_fit, c(args, heb_settings))
    }
  }
  global_paired <- fit_one("paired", pop_ref, unbridged)
  global <- if (method == "paired") global_paired else fit_one("heb", pop_ref, unbridged)

  # Map: trained on reference cells after the global paired shift.
  mapped_ref <- .heb_node_mapping_data(
    global_paired,
    x_ref,
    samples,
    pop_ref,
    mapping_markers,
    transform,
    ids_ref
  )
  set.seed(seed)
  train_rows <- unlist(
    lapply(split(seq_along(ids_ref), ids_ref), function(i) {
      i[sample.int(length(i), min(as.integer(train_cells), length(i)))]
    }),
    use.names = FALSE
  )
  scaled <- scale(mapped_ref[train_rows, , drop = FALSE])
  if (!all(is.finite(scaled))) {
    stop(
      "A mapping marker has zero variance in the training cells.",
      call. = FALSE
    )
  }
  som <- kohonen::som(
    scaled,
    grid = kohonen::somgrid(grid[[1]], grid[[2]]),
    rlen = rlen,
    keep.data = FALSE
  )
  map <- list(
    centre = attr(scaled, "scaled:center"),
    scale = attr(scaled, "scaled:scale"),
    nodes = som$codes[[1]],
    labels = paste0("SOM", seq_len(nrow(som$codes[[1]]))),
    mapping_markers = mapping_markers
  )

  # Node fits: every reference cell labelled by its nearest node.
  node_ref <- map$labels[.heb_node_nearest(mapped_ref, map)]
  # Node fits abstain on batches without node summaries; the fallback covers them.
  local <- fit_one(method, node_ref, "abstain")

  model <- structure(
    list(
      schema_version = .heb_node_schema_version,
      method = method,
      settings = list(
        grid = as.integer(grid),
        train_cells = as.integer(train_cells),
        rlen = as.integer(rlen),
        seed = seed,
        assignment = assignment,
        fallback = fallback,
        heb_settings = heb_settings
      ),
      markers = markers,
      transform = global$transform,
      anchor = anchor,
      map = map,
      global_paired = global_paired,
      global = global,
      local = local,
      software = list(R = R.version.string, heb = .heb_pkg_version()),
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      digest = NA_character_
    ),
    class = "heb_node_model"
  )
  model$digest <- .heb_digest(model)
  model
}

# Mapping markers after the global paired shift; cells the global model abstains on stay raw.
.heb_node_mapping_data <- function(
  global_paired,
  x,
  samples,
  population,
  mapping_markers,
  transform,
  sample_ids
) {
  res <- heb_integrate(
    global_paired,
    x,
    samples,
    population,
    transform = transform,
    unsupported = "abstain",
    sample_ids = sample_ids
  )
  out <- if (is.matrix(res$x)) {
    res$x[, mapping_markers, drop = FALSE]
  } else {
    as.matrix(as.data.frame(res$x)[mapping_markers])
  }
  storage.mode(out) <- "double"
  out
}

# Index of the nearest node for each row of `mapped` (squared Euclidean distance, ties to the first node).
.heb_node_nearest <- function(mapped, map) {
  z <- sweep(sweep(mapped, 2, map$centre, "-"), 2, map$scale, "/")
  d <- outer(rowSums(z^2), rowSums(map$nodes^2), "+") -
    2 * tcrossprod(z, map$nodes)
  d[d < 0] <- 0
  max.col(-d, ties.method = "first")
}

# Stop unless `model` is an intact heb_node_model with intact component models.
.heb_validate_node <- function(model) {
  if (!identical(model$schema_version, .heb_node_schema_version)) {
    stop(
      sprintf(
        "Unsupported schema_version '%s'; this package reads '%s'.",
        as.character(model$schema_version),
        .heb_node_schema_version
      ),
      call. = FALSE
    )
  }
  if (!identical(model$digest, .heb_digest(model))) {
    stop(
      "Digest mismatch: the model changed after it was frozen. Refit, or load the original file.",
      call. = FALSE
    )
  }
  .heb_validate(model$global_paired)
  .heb_validate(model$global)
  .heb_validate(model$local)
  invisible(model)
}

# heb_integrate() for a heb_node_model: nearest-node effects, with the global model as fallback.
.heb_integrate_nodes <- function(
  model,
  x,
  samples,
  population,
  transform,
  unsupported,
  sample_ids
) {
  .heb_validate_node(model)
  d <- .heb_cells(x, samples, population, model$markers, sample_ids)
  ids <- d$cells$sample_id
  mapped <- .heb_node_mapping_data(
    model$global_paired,
    x,
    samples,
    population,
    model$map$mapping_markers,
    transform,
    ids
  )
  nodes <- model$map$labels[.heb_node_nearest(mapped, model$map)]

  local <- heb_integrate(
    model$local,
    x,
    samples,
    nodes,
    transform,
    "abstain",
    sample_ids = ids
  )
  global <- heb_integrate(
    model$global,
    x,
    samples,
    population,
    transform,
    "abstain",
    sample_ids = ids
  )
  get_xm <- function(r) {
    if (is.matrix(r$x)) {
      r$x[, model$markers, drop = FALSE]
    } else {
      as.matrix(as.data.frame(r$x)[model$markers])
    }
  }
  out <- get_xm(local)
  global_xm <- get_xm(global)

  # Cells whose (batch, node, marker) has no local effect take the global value (fallback "global").
  actions <- local$actions
  abstained <- actions$action == "abstain"
  act_key <- paste(
    actions$batch,
    actions$population,
    actions$marker,
    sep = "\r"
  )
  if (model$settings$fallback == "global") {
    global_key <- paste(
      global$actions$batch,
      global$actions$population,
      global$actions$marker,
      sep = "\r"
    )
    for (m in model$markers) {
      i <- match(paste(d$cells$batch, nodes, m, sep = "\r"), act_key)
      g <- match(
        paste(d$cells$batch, d$cells$population, m, sep = "\r"),
        global_key
      )
      use_global <- abstained[i] & global$actions$action[g] == "apply"
      out[use_global, m] <- global_xm[use_global, m]
    }
    # A node without an effect is reported as a fallback when the batch has a global effect for that marker.
    has_global <- paste(actions$batch, actions$marker) %in%
      paste(
        global$actions$batch[global$actions$action == "apply"],
        global$actions$marker[global$actions$action == "apply"]
      )
    fb <- abstained & has_global
    actions$action[fb] <- "apply"
    actions$reason_code[fb] <- "node_fallback_global"
    actions$reason[
      fb
    ] <- "Node has no effect in this batch; the global effect was applied."
  }

  if (unsupported == "error" && any(actions$action == "abstain")) {
    ab <- actions[actions$action == "abstain", ]
    stop(
      sprintf(
        "%d (batch, node, marker) combination(s) cannot be corrected, e.g. %s / %s / %s: %s",
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
      x = .heb_restore(x, out, model$markers),
      actions = actions,
      model = model,
      nodes = nodes,
      excluded = local$excluded,
      aliquot_source = d$aliquot_source
    ),
    class = "heb_result"
  )
}

#' @export
print.heb_node_model <- function(x, ...) {
  cat(
    "<heb_node_model>",
    x$method,
    "| anchor:",
    x$anchor,
    "| nodes:",
    nrow(x$map$nodes),
    sprintf("(%d x %d)", x$settings$grid[[1]], x$settings$grid[[2]]),
    "| fallback:",
    x$settings$fallback,
    "| digest:",
    x$digest,
    "\n"
  )
  invisible(x)
}
