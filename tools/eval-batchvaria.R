# Purpose: Judge batch correction on the public simulations with BatchVaria 1.0.0. Decisions BV1-BV4. Two steps:
#          "simulate" (R 4.5) simulates, corrects with every arm and writes per-sample feature tables;
#          "evaluate" (R 4.6) builds one SummarizedExperiment per scenario and seed (features x study samples,
#          one assay per arm) and runs BatchVaria. Three feature sets are reported side by side: per-sample
#          medians of each marker in each simulated population, in each node of ONE shared SOM (fit on the
#          uncorrected cells of all samples; every arm's cells are assigned to it), and per-sample population
#          frequencies (true labels, so no-op under correction: a control).
#          Cross-sectional: the 4 scenarios of tools/sim-compare.R, design ~ batch + condition.
#          Longitudinal: the simulation of tools/sim-longitudinal.R (staggered, blocked, blocked_bridge_drift,
#          calibrated), design ~ batch + timepoint + subject (BV3).
# Run:     R 4.5 (R_LIBS = bridgenorm100 + 4.5 library), package root, as a job array over tasks 1-21:
#            Rscript tools/eval-batchvaria.R simulate <task> <run_dir>
#          R 4.6 (conda env R46_bioc324 first on PATH, LD_LIBRARY_PATH to its lib, R_LIBS unset):
#            Rscript tools/eval-batchvaria.R evaluate <run_dir>
#            Rscript tools/eval-batchvaria.R evaluate <run_dir> batch_timepoint   longitudinal only, ~ batch + timepoint,
#                                                                                 outputs in run_dir/batch_timepoint/
# Inputs:  tools/sim-compare.R and tools/sim-longitudinal.R (parsed, their command-line blocks skipped)
# Outputs (run_dir/):
#   features/task_NN.rds   task row, study-sample colData, and per feature set a list of arm = feature x sample matrix
#                          (node_median_<side>x<side> for each shared-SOM grid)
#   node_grid_check.csv   per task and grid: features in any arm and features complete in all arms
#   feature_filter.csv     features kept after dropping any with an NA in some arm (empty SOM node in a sample)
#   variance_percent.csv   anova percent per term, scenario, seed, feature set, arm (varianceTable)
#   variance_delta.csv     the same as delta versus uncorrected (varianceTable)
#   assay_variance.csv     total variance per arm and its ratio to uncorrected (assayVariance)
#   basis_retention.csv    retention on the uncorrected principal axes (basisRetention)
#   summary.csv            mean over seeds per scenario, feature set, arm
#   effect_feature_*.csv   the same anova percent on the one feature with the true group effect (P3|M05), per seed / mean
#   figures/               plotVarianceDelta (mean delta over seeds) and plotVarianceComposition (seed 1), PNG

cli <- commandArgs(trailingOnly = TRUE)

# Parameters -------------------------------------------------------------------
EVAL_SEEDS <- 1:3 # simulation seeds
SOM_SEED <- 473 # SOM seed of every method and of the shared SOM
# Side lengths of the square shared-SOM grids saved (8 = the methods' grid, GRID in tools/sim-compare.R). A sample with
# no cell in a node has no median there, and BatchVaria needs a full features x samples matrix, so a node-marker
# feature is kept only if every study sample has cells in that node in every arm. The evaluate step uses the finest
# 5x5 grid (BV5, user choice); node_grid_check.csv reports how many node-marker features each grid keeps.
SHARED_SOM_SIDES <- c(8, 6, 5, 4, 3)
NODE_GRID_SIDE <- 5 # BV5: the shared-SOM grid used for node_median
SHARED_SOM_TRAIN_CELLS <- 500 # cells per sample, bridgenorm's default train_cells
SHARED_SOM_RLEN <- 100 # SOM passes, bridgenorm's default rlen
CROSS_FORMULA <- ~ batch + condition
LONG_FORMULA <- ~ batch + timepoint + subject # BV3; fall back to ~ batch + timepoint if BatchVaria rejects it
BASELINE <- "uncorrected"
EFFECT_FEATURE <- "P3|M05" # the one feature with a true group effect (EFFECT_POP, EFFECT_MARKER in tools/sim-compare.R)
FEATURE_SETS <- c("pop_median", "node_median", "pop_freq") # node_median = the chosen grid

eval_tasks <- rbind(
  data.frame(
    expand.grid(
      seed = EVAL_SEEDS,
      design = c("balanced", "confounded"),
      batch_effect = c("global", "calibrated"),
      stringsAsFactors = FALSE
    ),
    kind = "cross"
  ),
  data.frame(
    expand.grid(
      seed = EVAL_SEEDS,
      design = c("staggered", "blocked", "blocked_bridge_drift"),
      batch_effect = "calibrated",
      stringsAsFactors = FALSE
    ),
    kind = "long"
  )
)
eval_tasks$scenario <- ifelse(
  eval_tasks$kind == "cross",
  paste(eval_tasks$batch_effect, eval_tasks$design, sep = "_"),
  paste0("long_", eval_tasks$design)
)

# ==============================================================================
# simulate (R 4.5)
# ==============================================================================

if (cli[[1]] == "simulate") {
  suppressPackageStartupMessages({
    library(bridgenorm)
    library(dplyr)
    library(tidyr)
  })
  task_id <- as.integer(cli[[2]])
  run_dir <- cli[[3]]
  tk <- eval_tasks[task_id, ]
  out_file <- file.path(run_dir, "features", sprintf("task_%02d.rds", task_id))
  stopifnot("task output is new" = !file.exists(out_file))

  # Design constants and simulate() -------------------------------------------
  # tools/sim-compare.R without its command-line blocks (the same parse as tools/sim-longitudinal.R).
  # `args` stays the base function, so the "repro" test inside sim-compare.R is FALSE, as in the original runs.
  for (e in parse("tools/sim-compare.R")) {
    txt <- paste(deparse(e), collapse = "")
    if (grepl("^args <- commandArgs", txt) || grepl("^if \\(args", txt)) {
      next
    }
    eval(e, globalenv())
  }

  # Simulation ---------------------------------------------------------------
  # Cross-sectional: simulate(). Longitudinal: the run block of tools/sim-longitudinal.R up to the methods, run
  # twice with the same random stream; the second run replaces every batch effect by the identity, which gives
  # the simulated values without batch effect (the oracle), cell for cell.
  if (tk$kind == "cross") {
    sim <- simulate(tk$seed, tk$batch_effect, tk$design)
    study_ids <- sim$samples$sample_id[!is.na(sim$samples$condition)]
    coldata_cols <- "condition"
  } else {
    lg_exprs <- parse("tools/sim-longitudinal.R")
    lg_first_line <- vapply(lg_exprs, \(e) deparse(e)[[1]], "")
    lg_top <- lg_exprs[
      !grepl(
        "^(if \\(args|for \\(e in parse|args <- commandArgs)",
        lg_first_line
      )
    ]
    for (e in lg_top) {
      eval(e, globalenv())
    }
    run_block <- as.list(lg_exprs[[grep(
      "^if \\(args\\[\\[1\\]\\] == \"run\"",
      lg_first_line
    )]][[3]])[-1]
    stmt_lhs <- vapply(
      run_block,
      \(e) {
        if (is.call(e) && identical(e[[1]], quote(`<-`))) {
          deparse(e[[2]])[[1]]
        } else {
          ""
        }
      },
      ""
    )
    methods_at <- grep("^fit_samples$", stmt_lhs)[[1]]
    sim_stmts <- run_block[setdiff(
      seq_len(methods_at - 1),
      which(stmt_lhs %in% c("task", "tmp", "tk"))
    )]
    simulate_long <- function(identity_effect) {
      env <- new.env(parent = globalenv())
      env$tk <- tk
      for (st in sim_stmts) {
        eval(st, env)
        if (identity_effect && identical(deparse(st[[2]])[[1]], "effect")) {
          env$effect <- lapply(env$effect, \(f) {
            list(shift = f$shift * 0, scale = f$scale^0)
          })
        }
      }
      list(samples = env$s, cells = env$x)
    }
    sim <- simulate_long(FALSE)
    oracle <- simulate_long(TRUE)
    stopifnot(
      "oracle has the same samples and populations" = identical(
        sim$cells$pop,
        oracle$cells$pop
      )
    )
    sim$truth <- as.matrix(oracle$cells[MARKERS])
    sim$samples$condition <- NA_character_
    sim$cells$cell_id <- seq_len(nrow(sim$cells))
    study_ids <- sim$samples$sample_id[!is.na(sim$samples$subject)]
    coldata_cols <- c("timepoint", "subject")
    rm(oracle)

    # No cyCombine_covar arm in the longitudinal designs (BV: cross-sectional only)
    run_body <- body(run_methods)
    covar_at <- which(vapply(
      run_body,
      \(e) is.call(e) && identical(e[[2]], quote(cyc_args)),
      NA
    ))
    stopifnot("one cyc_args statement" = length(covar_at) == 1)
    run_body[[covar_at]] <- quote(cyc_args <- list(cyCombine = list()))
    body(run_methods) <- run_body
  }

  # Methods ------------------------------------------------------------------
  corrected <- run_methods(sim, SOM_SEED)

  # Shared SOM ---------------------------------------------------------------
  # Trained on the uncorrected cells of all samples (bridges included), up to SHARED_SOM_TRAIN_CELLS per sample,
  # standardised; every arm's cells go to the nearest node.
  cells <- sim$cells
  set.seed(SOM_SEED)
  train_rows <- unlist(
    lapply(split(seq_len(nrow(cells)), cells$sample_id), \(i) {
      sample(i, min(SHARED_SOM_TRAIN_CELLS, length(i)))
    }),
    use.names = FALSE
  )
  som_x <- scale(corrected[[BASELINE]][train_rows, ])
  shared_soms <- lapply(
    setNames(
      SHARED_SOM_SIDES,
      sprintf("node_median_%dx%d", SHARED_SOM_SIDES, SHARED_SOM_SIDES)
    ),
    \(side) {
      set.seed(SOM_SEED)
      kohonen::som(
        som_x,
        grid = kohonen::somgrid(side, side),
        rlen = SHARED_SOM_RLEN,
        keep.data = TRUE # predict() needs the training data
      )
    }
  )
  assign_node <- function(m, som) {
    z <- scale(
      m,
      center = attr(som_x, "scaled:center"),
      scale = attr(som_x, "scaled:scale")
    )
    ok <- complete.cases(z) # cyCombine leaves a few cells NA
    node <- rep(NA_integer_, nrow(z))
    node[ok] <- predict(som, newdata = z[ok, , drop = FALSE])$unit.classif
    ifelse(is.na(node), NA_character_, sprintf("SOM%02d", node))
  }

  # Feature tables -----------------------------------------------------------
  # One matrix per arm: rows "<group>|<marker>", columns study samples in colData order.
  coldata <- sim$samples[
    match(study_ids, sim$samples$sample_id),
    c("sample_id", "batch", coldata_cols)
  ]
  rownames(coldata) <- coldata$sample_id
  group_medians <- function(m, group) {
    tibble(sample_id = cells$sample_id, group = group, as_tibble(m)) |>
      filter(sample_id %in% study_ids, !is.na(group)) |>
      summarise(
        across(all_of(MARKERS), \(v) median(v, na.rm = TRUE)),
        .by = c(sample_id, group)
      ) |>
      pivot_longer(all_of(MARKERS), names_to = "marker") |>
      mutate(feature = paste(group, marker, sep = "|")) |>
      pivot_wider(
        id_cols = feature,
        names_from = sample_id,
        values_from = value
      ) |>
      tibble::column_to_rownames("feature") |>
      as.matrix() |>
      _[, study_ids]
  }
  pop_freq <- cells |>
    filter(sample_id %in% study_ids) |>
    count(sample_id, pop) |>
    mutate(freq = n / sum(n), .by = sample_id) |>
    pivot_wider(
      id_cols = pop,
      names_from = sample_id,
      values_from = freq,
      values_fill = 0
    ) |>
    tibble::column_to_rownames("pop") |>
    as.matrix() |>
    _[, study_ids]
  features <- list(
    pop_median = lapply(corrected, group_medians, group = cells$pop),
    pop_freq = lapply(corrected, \(m) pop_freq)
  )
  for (nm in names(shared_soms)) {
    features[[nm]] <- lapply(corrected, \(m) {
      group_medians(m, assign_node(m, shared_soms[[nm]]))
    })
  }
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(task = tk, coldata = coldata, features = features), out_file)
  quit(save = "no")
}

# ==============================================================================
# evaluate (R 4.6)
# ==============================================================================

suppressPackageStartupMessages({
  library(BatchVaria)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})
run_dir <- cli[[2]]
# Second pass `evaluate <run_dir> batch_timepoint`: longitudinal designs only, ~ batch + timepoint (see LONG_FORMULA)
fallback <- length(cli) >= 3 && cli[[3]] == "batch_timepoint"
if (fallback) {
  LONG_FORMULA <- ~ batch + timepoint
}
out_dir <- if (fallback) file.path(run_dir, "batch_timepoint") else run_dir
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Shared-SOM grid -------------------------------------------------------------------
# Pad each arm to the union of features (a node can be empty in every sample of one arm), then keep features
# without an NA in any arm. Used for the grid choice and for every SummarizedExperiment.
pad_arms <- function(mats) {
  all_rows <- Reduce(union, map(mats, rownames))
  mats <- map(mats, \(m) {
    padded <- matrix(
      NA_real_,
      length(all_rows),
      ncol(m),
      dimnames = list(all_rows, colnames(m))
    )
    padded[rownames(m), ] <- m
    padded
  })
  list(mats = mats, keep = Reduce(`&`, map(mats, \(m) rowSums(is.na(m)) == 0)))
}
read_task <- \(task_id) {
  readRDS(file.path(run_dir, "features", sprintf("task_%02d.rds", task_id)))
}
grid_check <- map(seq_len(nrow(eval_tasks)), \(i) {
  feats <- read_task(i)$features
  map(feats[grepl("^node_median_", names(feats))], \(mats) {
    padded <- pad_arms(mats)
    tibble(
      n_features_union = length(padded$keep),
      n_features_complete = sum(padded$keep)
    )
  }) |>
    list_rbind(names_to = "grid") |>
    mutate(
      task_id = i,
      scenario = eval_tasks$scenario[[i]],
      seed = eval_tasks$seed[[i]]
    )
}) |>
  list_rbind()
write.csv(
  grid_check,
  file.path(out_dir, "node_grid_check.csv"),
  row.names = FALSE
)
NODE_GRID <- sprintf("node_median_%dx%d", NODE_GRID_SIDE, NODE_GRID_SIDE)
message("Shared SOM grid used for node_median: ", NODE_GRID)

# Build one SummarizedExperiment per task ----------------------------------------
build_se <- function(task_id, feature_set) {
  res <- read_task(task_id)
  padded <- pad_arms(res$features[[
    if (feature_set == "node_median") NODE_GRID else feature_set
  ]])
  cd <- res$coldata
  cd[] <- lapply(cd, \(v) if (is.character(v)) factor(v) else v)
  se <- SummarizedExperiment(
    assays = map(padded$mats, \(m) m[padded$keep, rownames(cd)]),
    colData = S4Vectors::DataFrame(cd)
  )
  list(se = se, n_total = length(padded$keep))
}

# Variance profile, tables and plot inputs for one SE -----------------------------
profile_se <- function(se, formula, with_basis = TRUE) {
  se <- profileVariance(se, formula, methods = "anova", assays = assayNames(se))
  vt <- varianceTable(
    se,
    assays = assayNames(se),
    method = "anova",
    baseline = BASELINE
  )
  list(
    se = se,
    percent = vt$percent,
    delta = vt$delta,
    assay_variance = assayVariance(se),
    retention = if (with_basis) basisRetention(se, reference = BASELINE), # no principal axes for one feature
    eval = evaluateCorrections(
      se,
      baseline = BASELINE,
      diagnostics = "variance",
      varianceMethod = "anova"
    )
  )
}

# Profile every task and feature set -----------------------------------------------
runs <- expand.grid(
  task_id = which(eval_tasks$kind == "long" | !fallback),
  feature_set = FEATURE_SETS,
  stringsAsFactors = FALSE
) |>
  as_tibble()
profiles <- pmap(runs, \(task_id, feature_set) {
  tk <- eval_tasks[task_id, ]
  built <- build_se(task_id, feature_set)
  se <- built$se
  n_total <- built$n_total
  pr <- profile_se(se, if (tk$kind == "cross") CROSS_FORMULA else LONG_FORMULA)
  pr$info <- tibble(
    scenario = tk$scenario,
    seed = tk$seed,
    feature_set = feature_set,
    n_features_total = n_total,
    n_features_kept = nrow(se)
  )
  pr
}) |>
  set_names(paste(runs$task_id, runs$feature_set))
tag <- \(tbl, key) {
  bind_cols(profiles[[key]]$info[c("scenario", "seed", "feature_set")], tbl)
}
all_keys <- names(profiles)

# Tables --------------------------------------------------------------------------
# percent and delta (varianceTable), total variance (assayVariance), retention (basisRetention): one long table each
feature_filter <- map(all_keys, \(k) profiles[[k]]$info) |> list_rbind()
variance_percent <- map(all_keys, \(k) {
  profiles[[k]]$percent |>
    pivot_longer(-component, names_to = "assay", values_to = "percent") |>
    tag(k)
}) |>
  list_rbind()
variance_delta <- map(all_keys, \(k) {
  profiles[[k]]$delta |>
    pivot_longer(
      -component,
      names_to = "assay",
      values_to = "delta_percent_points"
    ) |>
    tag(k)
}) |>
  list_rbind()
assay_variance <- map(all_keys, \(k) tag(profiles[[k]]$assay_variance, k)) |>
  list_rbind() |>
  mutate(
    total_variance_ratio = total_variance / total_variance[assay == BASELINE],
    .by = c(scenario, seed, feature_set)
  )
basis_retention <- map(all_keys, \(k) tag(profiles[[k]]$retention, k)) |>
  list_rbind()
write.csv(
  feature_filter,
  file.path(out_dir, "feature_filter.csv"),
  row.names = FALSE
)
write.csv(
  variance_percent,
  file.path(out_dir, "variance_percent.csv"),
  row.names = FALSE
)
write.csv(
  variance_delta,
  file.path(out_dir, "variance_delta.csv"),
  row.names = FALSE
)
write.csv(
  assay_variance,
  file.path(out_dir, "assay_variance.csv"),
  row.names = FALSE
)
write.csv(
  basis_retention,
  file.path(out_dir, "basis_retention.csv"),
  row.names = FALSE
)

# Summary: mean over seeds ---------------------------------------------------------
# group = condition (cross-sectional) or timepoint (longitudinal). *_abs_ratio = term percent x total variance,
# relative to the uncorrected assay, so it is not changed by the compositional rescaling (assayVariance help).
summary_seed <- variance_percent |>
  mutate(
    component = if_else(
      component %in% c("condition", "timepoint"),
      "group",
      component
    )
  ) |>
  pivot_wider(names_from = component, values_from = percent) |>
  left_join(
    assay_variance[c(
      "scenario",
      "seed",
      "feature_set",
      "assay",
      "total_variance",
      "total_variance_ratio"
    )],
    by = c("scenario", "seed", "feature_set", "assay"),
    relationship = "one-to-one"
  ) |>
  left_join(
    basis_retention[c("scenario", "seed", "feature_set", "assay", "retention")],
    by = c("scenario", "seed", "feature_set", "assay"),
    relationship = "one-to-one"
  ) |>
  mutate(
    batch_abs_ratio = batch *
      total_variance /
      (batch * total_variance)[assay == BASELINE],
    group_abs_ratio = group *
      total_variance /
      (group * total_variance)[assay == BASELINE],
    shared_negative = shared < 0,
    .by = c(scenario, seed, feature_set)
  )
summary_table <- summary_seed |>
  summarise(
    n_seeds = n(),
    batch_pct = mean(batch),
    group_pct = mean(group),
    subject_pct = mean(subject),
    shared_pct = mean(shared),
    shared_negative = mean(shared) < 0,
    shared_negative_seeds = sum(shared_negative),
    residual_pct = mean(residual),
    total_variance_ratio = mean(total_variance_ratio),
    batch_abs_ratio = mean(batch_abs_ratio),
    group_abs_ratio = mean(group_abs_ratio),
    basis_retention = mean(retention),
    .by = c(scenario, feature_set, assay)
  ) |>
  rename(method = assay) |>
  arrange(scenario, feature_set, method)
write.csv(summary_table, file.path(out_dir, "summary.csv"), row.names = FALSE)

# Effect feature alone ------------------------------------------------------------------
# The group effect sits in one pop_median feature (EFFECT_FEATURE, in every task: cross-sectional condition, and the
# longitudinal time effect), so the mean over all features dilutes it. The same variance profile on that feature only.
effect_feature <- map(which(runs$feature_set == "pop_median"), \(i) {
  tk <- eval_tasks[runs$task_id[[i]], ]
  se <- build_se(runs$task_id[[i]], "pop_median")$se[EFFECT_FEATURE, ]
  pr <- profile_se(
    se,
    if (tk$kind == "cross") CROSS_FORMULA else LONG_FORMULA,
    with_basis = FALSE
  )
  pr$percent |>
    pivot_longer(-component, names_to = "assay", values_to = "percent") |>
    mutate(scenario = tk$scenario, seed = tk$seed)
}) |>
  list_rbind()
effect_feature_summary <- effect_feature |>
  summarise(percent = mean(percent), .by = c(scenario, assay, component)) |>
  pivot_wider(names_from = component, values_from = percent) |>
  rename(method = assay) |>
  arrange(scenario, method)
write.csv(
  effect_feature,
  file.path(out_dir, "effect_feature_percent.csv"),
  row.names = FALSE
)
write.csv(
  effect_feature_summary,
  file.path(out_dir, "effect_feature_summary.csv"),
  row.names = FALSE
)

# Figures ---------------------------------------------------------------------------
# Delta: mean over seeds of evaluateCorrections()'s variance_delta, drawn by plotVarianceDelta().
# Composition: plotVarianceComposition() of seed 1.
for (sc in unique(eval_tasks$scenario[runs$task_id])) {
  for (fs in FEATURE_SETS) {
    keys <- all_keys[map_lgl(all_keys, \(k) {
      profiles[[k]]$info$scenario == sc && profiles[[k]]$info$feature_set == fs
    })]
    delta_mean <- map(keys, \(k) profiles[[k]]$eval$variance_delta) |>
      list_rbind() |>
      summarise(delta = mean(delta), .by = c(assay, method, term))
    p_delta <- plotVarianceDelta(list(variance_delta = delta_mean)) +
      ggtitle(sprintf(
        "%s, %s: variance change vs %s (mean of %d seeds)",
        sc,
        fs,
        BASELINE,
        length(keys)
      ))
    ggsave(
      file.path(fig_dir, sprintf("delta_%s_%s.png", sc, fs)),
      p_delta,
      width = 9,
      height = 5,
      dpi = 150,
      bg = "white"
    )
    first <- keys[map_lgl(keys, \(k) {
      profiles[[k]]$info$seed == min(EVAL_SEEDS)
    })]
    p_comp <- plotVarianceComposition(profiles[[first]]$se, method = "anova") +
      ggtitle(sprintf(
        "%s, %s: variance composition (seed %d)",
        sc,
        fs,
        min(EVAL_SEEDS)
      ))
    ggsave(
      file.path(fig_dir, sprintf("composition_%s_%s.png", sc, fs)),
      p_comp,
      width = 9,
      height = 5,
      dpi = 150,
      bg = "white"
    )
  }
}
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo_R46.txt")
)
