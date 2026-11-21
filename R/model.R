# Frozen model: fields, digest, validation, print, save and load.
#
# Ported from SpectraWeaveR R/corrector.R (frozen model, digest, save and load).

# Version tag of the frozen model layout; bump when fields change.
.heb_schema_version <- "heb_model/1"

# Fields every frozen model carries.
.heb_model_fields <- c(
  "schema_version",
  "method",
  "settings",
  "markers",
  "transform",
  "population_levels",
  "anchor",
  "min_cells",
  "hierarchy",
  "known_batches",
  "bridges",
  "effects",
  "targets",
  "excluded",
  "aliquot_source",
  "software",
  "created",
  "parent_digest",
  "digest"
)

# Package version, or NA when the code is sourced rather than installed.
.heb_pkg_version <- function() {
  tryCatch(as.character(utils::packageVersion("heb")), error = function(e) {
    NA_character_
  })
}

# MD5 of a canonical text form of every field except the digest itself.
# hexNumeric deparse keeps doubles exact, so the hash moves with any change.
.heb_digest <- function(model) {
  body <- unclass(model)
  body$digest <- NULL
  txt <- deparse(
    body,
    control = c(
      "keepNA",
      "keepInteger",
      "niceNames",
      "showAttributes",
      "hexNumeric"
    )
  )
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f))
  writeLines(txt, f, useBytes = TRUE)
  unname(tools::md5sum(f))
}

# Build and digest the frozen model.
.heb_freeze <- function(
  method,
  settings,
  com,
  anchor,
  min_cells,
  hierarchy,
  effects,
  targets,
  excluded
) {
  model <- structure(
    list(
      schema_version = .heb_schema_version,
      method = method,
      settings = settings,
      markers = colnames(com$d$xm),
      transform = com$transform,
      population_levels = sort(unique(com$d$cells$population)),
      anchor = anchor,
      min_cells = as.integer(min_cells),
      hierarchy = hierarchy,
      known_batches = com$batches,
      bridges = com$bridges,
      effects = effects,
      targets = targets,
      excluded = excluded,
      aliquot_source = com$d$aliquot_source,
      software = list(R = R.version.string, heb = .heb_pkg_version()),
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      parent_digest = NA_character_,
      digest = NA_character_
    ),
    class = "heb_model"
  )
  model$digest <- .heb_digest(model)
  model
}

# Stop unless `model` has the current schema and an intact digest.
.heb_validate <- function(model) {
  if (!inherits(model, "heb_model")) {
    stop("'model' must be a heb_model.", call. = FALSE)
  }
  absent <- setdiff(.heb_model_fields, names(model))
  if (length(absent) > 0) {
    stop(
      "Model is missing field(s): ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!identical(model$schema_version, .heb_schema_version)) {
    stop(
      sprintf(
        "Unsupported schema_version '%s'; this package reads '%s'.",
        as.character(model$schema_version),
        .heb_schema_version
      ),
      call. = FALSE
    )
  }
  if (!all(.heb_effect_cols %in% names(model$effects))) {
    stop("Model 'effects' lacks required columns.", call. = FALSE)
  }
  if (!identical(model$digest, .heb_digest(model))) {
    stop(
      "Digest mismatch: the model changed after it was frozen. ",
      "Refit, or load the original file.",
      call. = FALSE
    )
  }
  invisible(model)
}

# Validate either model type.
.heb_validate_any <- function(model) {
  if (inherits(model, "heb_node_quantile_model")) {
    .heb_validate_node_quantile(model)
  } else if (inherits(model, "heb_node_model")) {
    .heb_validate_node(model)
  } else {
    .heb_validate(model)
  }
}

#' @export
print.heb_model <- function(x, ...) {
  cat(
    "<heb_model>",
    x$method,
    "| anchor:",
    x$anchor,
    "| batches:",
    length(x$known_batches),
    "| markers:",
    length(x$markers),
    "| populations:",
    length(x$population_levels),
    "\n"
  )
  cat(
    "Transform:",
    .heb_transform_label(x$transform),
    "| min_cells:",
    x$min_cells,
    "| digest:",
    x$digest,
    "\n"
  )
  n_unb <- sum(!x$bridges$bridged)
  cat(
    nrow(x$effects),
    "effect(s);",
    nrow(x$excluded),
    "reference summary(ies) excluded;",
    n_unb,
    "unbridged batch(es)\n"
  )
  invisible(x)
}

#' Save and load a frozen HEB model
#'
#' \code{heb_save()} checks the digest and writes the model with
#' \code{saveRDS()}. \code{heb_load()} reads it and checks the schema and
#' digest again, so a changed file fails to load.
#'
#' @param model A \code{heb_model}, \code{heb_node_model} or
#'   \code{heb_node_quantile_model}.
#' @param path File path.
#' @return \code{heb_save()} returns \code{path} invisibly; \code{heb_load()}
#'   returns the model.
#' @export
heb_save <- function(model, path) {
  .heb_stop_missing(c(model = missing(model), path = missing(path)))
  .heb_validate_any(model)
  saveRDS(model, path)
  invisible(path)
}

#' @rdname heb_save
#' @export
heb_load <- function(path) {
  .heb_stop_missing(c(path = missing(path)))
  model <- readRDS(path)
  .heb_validate_any(model)
  model
}
