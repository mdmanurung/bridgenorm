# Node-level correction (internal): .bn_fit_nodes(), the pass-1 and pass-2 model of bn_fit(), and its apply path.
#

.bn_node_schema_version <- "bn_node_model/1"

#' Fit a node-level batch-correction model
#'
#' Batch effects often differ between cell states. \code{.bn_fit_nodes()}
#' learns a self-organising map (SOM) on the reference cells and fits the
#' paired shift separately in each
#' map node, with the node as the population. \code{\link{bn_apply}} then
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
#' The training cells are drawn by position within each reference sample, so
#' for a given \code{seed} the map depends on the row order of \code{x}
#' within a sample. \code{.bn_fit_nodes()} restores the caller's random
#' number stream on exit.
#'
#' Every argument that changes an estimate is required, with no default; the
#' \pkg{kohonen} settings left at their defaults (topology, learning rate,
#' radius, distance) and the \pkg{kohonen} version are stored in
#' \code{settings$som}.
#' Needs the \pkg{kohonen} package.
#'
#' @inheritParams bn_fit
#' @param method Only \code{"paired"} (\code{\link{bn_fit_paired}}).
#' @param eb_settings Only \code{NULL}.
#' @param grid SOM grid as \code{c(xdim, ydim)}.
#' @param mapping_markers Markers that define cell state for the map; each
#'   must be one of \code{markers}, so it is pre-shifted.
#' @param train_cells Maximum training cells per reference sample.
#' @param rlen Number of SOM training passes (\code{kohonen::som}).
#' @param seed Seed for the training-cell draw and the SOM.
#' @param assignment How cells are assigned to nodes; only \code{"nearest"}.
#' @param fallback \code{"global"} or \code{"abstain"}; see Details.
#' @return A \code{bn_node_model}: a list with \code{schema_version},
#'   \code{method}, \code{settings}, \code{markers}, \code{transform},
#'   \code{anchor}, \code{map} (\code{centre}, \code{scale}, \code{nodes},
#'   \code{labels}, \code{mapping_markers}), \code{global_paired},
#'   \code{global} and \code{local} (frozen \code{bn_model}s),
#'   \code{software}, \code{created} and \code{digest}.
#' @noRd
.bn_fit_nodes <- function(
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
  sample_ids = NULL,
  target = "anchor"
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
    fallback = missing(fallback)
  ))
  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop(
      ".bn_fit_nodes() needs the 'kohonen' package: install.packages(\"kohonen\").",
      call. = FALSE
    )
  }
  .bn_check_choice(
    method,
    "method",
    "paired",
    what = ".bn_fit_nodes"
  )
  .bn_check_choice(assignment, "assignment", "nearest", what = ".bn_fit_nodes")
  .bn_check_choice(
    fallback,
    "fallback",
    c("global", "abstain"),
    what = ".bn_fit_nodes"
  )
  if (!is.null(eb_settings)) {
    stop("method = \"paired\" takes eb_settings = NULL.", call. = FALSE)
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
  # The pre-shift corrects `markers` only, so a mapping marker outside them would be mapped raw.
  uncorrected <- setdiff(mapping_markers, markers)
  if (length(uncorrected) > 0) {
    stop(
      "Every one of 'mapping_markers' must be in 'markers', so the map sees pre-shifted values; not in 'markers': ",
      paste(uncorrected, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Reference rows only: every fit and the map use them.
  d <- .bn_cells(x, samples, population, markers, sample_ids)
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
      target = target,
      sample_ids = ids_ref
    )
    do.call(bn_fit_paired, args)
  }
  global_paired <- fit_one("paired", pop_ref, unbridged)
  global <- global_paired

  # Map: trained on reference cells after the global paired shift.
  mapped_ref <- .bn_node_mapping_data(
    global_paired,
    x_ref,
    samples,
    pop_ref,
    mapping_markers,
    transform,
    ids_ref
  )$x
  # Seed the draw and the SOM, then give the caller back their random number stream.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit(
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    },
    add = TRUE
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
  node_ref <- map$labels[.bn_node_nearest(mapped_ref, map)]
  # Node fits abstain on batches without node summaries; the fallback covers them.
  local <- fit_one(method, node_ref, "abstain")

  model <- structure(
    list(
      schema_version = .bn_node_schema_version,
      method = method,
      settings = list(
        grid = as.integer(grid),
        train_cells = as.integer(train_cells),
        rlen = as.integer(rlen),
        seed = seed,
        assignment = assignment,
        fallback = fallback,
        eb_settings = eb_settings,
        target = target,
        # kohonen settings left at their defaults, recorded so the map can be rebuilt.
        som = list(
          kohonen = as.character(utils::packageVersion("kohonen")),
          topo = som$grid$topo,
          alpha = som$alpha,
          radius = som$radius,
          dist_fcts = som$dist.fcts
        )
      ),
      markers = markers,
      transform = global$transform,
      anchor = anchor,
      map = map,
      global_paired = global_paired,
      global = global,
      local = local,
      software = list(R = R.version.string, bridgenorm = .bn_pkg_version()),
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      parent_digest = NA_character_,
      digest = NA_character_
    ),
    class = "bn_node_model"
  )
  model$digest <- .bn_digest(model)
  model
}

# Mapping markers after the global paired shift; cells the global model abstains on stay raw.
# Returns the mapped matrix (x) and the global paired model as bn_apply() returned it (model).
.bn_node_mapping_data <- function(
  global_paired,
  x,
  samples,
  population,
  mapping_markers,
  transform,
  sample_ids
) {
  res <- bn_apply(
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
  # The global model comes back too: it may carry effects of a batch estimated just now.
  list(x = out, model = res$model)
}

# Index of the nearest node for each row of `mapped` (squared Euclidean distance, ties to the first node).
.bn_node_nearest <- function(mapped, map) {
  z <- sweep(sweep(mapped, 2, map$centre, "-"), 2, map$scale, "/")
  d <- outer(rowSums(z^2), rowSums(map$nodes^2), "+") -
    2 * tcrossprod(z, map$nodes)
  d[d < 0] <- 0
  max.col(-d, ties.method = "first")
}

# Stop unless `model` is an intact bn_node_model with intact component models.
.bn_validate_node <- function(model) {
  if (!identical(model$schema_version, .bn_node_schema_version)) {
    stop(
      sprintf(
        "Unsupported schema_version '%s'; this package reads '%s'.",
        as.character(model$schema_version),
        .bn_node_schema_version
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
  .bn_validate(model$global_paired)
  .bn_validate(model$global)
  .bn_validate(model$local)
  invisible(model)
}

# bn_apply() for a bn_node_model: nearest-node effects, with the global model as fallback.
.bn_apply_nodes <- function(
  model,
  x,
  samples,
  population,
  transform,
  unsupported,
  sample_ids
) {
  .bn_validate_node(model)
  d <- .bn_cells(x, samples, population, model$markers, sample_ids)
  ids <- d$cells$sample_id
  mapping <- .bn_node_mapping_data(
    model$global_paired,
    x,
    samples,
    population,
    model$map$mapping_markers,
    transform,
    ids
  )
  nodes <- model$map$labels[.bn_node_nearest(mapping$x, model$map)]

  local <- bn_apply(
    model$local,
    x,
    samples,
    nodes,
    transform,
    "abstain",
    sample_ids = ids
  )
  global <- bn_apply(
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
    # Count, per action row, its cells and those that took the global effect (decided per cell population).
    n_row <- integer(nrow(actions))
    n_global <- integer(nrow(actions))
    for (m in model$markers) {
      i <- match(paste(d$cells$batch, nodes, m, sep = "\r"), act_key)
      g <- match(
        paste(d$cells$batch, d$cells$population, m, sep = "\r"),
        global_key
      )
      use_global <- abstained[i] &
        !is.na(g) &
        global$actions$action[g] == "apply"
      use_global[is.na(use_global)] <- FALSE
      out[use_global, m] <- global_xm[use_global, m]
      n_row <- n_row + tabulate(i[!is.na(i)], nbins = nrow(actions))
      n_global <- n_global + tabulate(i[use_global], nbins = nrow(actions))
    }
    # A row with no cells is a fallback when its batch has a global effect for that marker.
    has_global <- paste(actions$batch, actions$marker) %in%
      paste(
        global$actions$batch[global$actions$action == "apply"],
        global$actions$marker[global$actions$action == "apply"]
      )
    fb <- abstained &
      ((n_row > 0 & n_global == n_row) | (n_row == 0 & has_global))
    partial <- abstained & n_global > 0 & n_global < n_row
    actions$action[fb] <- "apply"
    actions$reason_code[fb] <- "node_fallback_global"
    actions$reason[
      fb
    ] <- "Node has no effect in this batch; the global effect was applied."
    actions$reason_code[partial] <- "node_fallback_partial"
    actions$reason[partial] <- sprintf(
      "Node has no effect in this batch; %d of %d cells took the global effect, the rest (populations without a global effect) kept raw values.",
      n_global[partial],
      n_row[partial]
    )
  }

  # Keep the component models as bn_apply() returned them: they carry effects of new batches.
  updated <- model
  updated$global_paired <- mapping$model
  updated$global <- if (model$method == "paired") {
    mapping$model
  } else {
    global$model
  }
  updated$local <- local$model
  if (!identical(updated, model)) {
    updated$parent_digest <- model$digest
    updated$digest <- .bn_digest(updated)
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
      x = .bn_restore(x, out, model$markers),
      actions = actions,
      model = updated,
      nodes = nodes,
      excluded = local$excluded,
      aliquot_source = d$aliquot_source
    ),
    class = "bn_result"
  )
}

#' @export
print.bn_node_model <- function(x, ...) {
  cat(
    "<bn_node_model>",
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
