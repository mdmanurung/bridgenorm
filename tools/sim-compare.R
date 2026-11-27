# Purpose: Simulated comparison of bridgenorm with CytoNorm and cyCombine for the "Compare with CytoNorm and
#          cyCombine" article. The design and the batch-effect sizes are set to realistic values (numbers below);
#          no data are used. Four scenarios (batch effect global only, or
#          global plus cell-state-specific shift and scale) x (condition balanced within batch, or confounded with
#          batch), 10 simulation seeds x 3 SOM seeds each. A task is one scenario x simulation seed x SOM seed
#          (120 tasks); the data depend on the simulation seed only, the methods on the SOM seed.
# Run:     Rscript tools/sim-compare.R run <task> <tmp_dir>   task 1-120, e.g. as a job array
#          Rscript tools/sim-compare.R run <task> <tmp_dir> repro   task 1-12: simulation seed 1-3 with SOM seed
#                                                                   equal to it, the seeding of the first release
#          Rscript tools/sim-compare.R combine <tmp_dir> [out_dir]  merge the task results (default OUT below)
#          (package root; needs bridgenorm, kohonen, CytoNorm, flowCore, cyCombine)
# Inputs:  none (simulated)
# Outputs (OUT, default vignettes/compare-data/):
#   scores.csv        batch_effect, design, sim_seed, som_seed, method, heldout_emd (held-out bridge aliquot: worst
#                     batch pair, median over markers), effect (estimated condition effect), effect_true,
#                     rmse_truth (study cells, all markers, to the simulated values without batch effect)
#   densities.csv     batch_effect, design, method, batch, value: the effect marker in effect-population cells of
#                     the held-out bridge aliquot, simulation seed 1, first SOM seed, 300 cells for 4 batches
#   settings.csv      setting, value
#   (package versions are recorded in settings.csv; no sessionInfo, which would carry local paths)

suppressPackageStartupMessages({
  library(bridgenorm)
})

# Parameters -------------------------------------------------------------------
SIM_SEEDS <- 1:10 # simulated data sets per scenario (decision P2c)
SOM_SEEDS <- 473:475 # SOM seed of every method, crossed with the simulation seeds (decision P2c)
N_BATCHES <- 15
BATCHES <- sprintf("b%02d", seq_len(N_BATCHES)) # b01 = anchor, no batch effect
N_REF <- 3 # bridge aliquots per batch; ref1 and ref2 fit, ref3 held out
N_STUDY <- 12 # study samples (donors) per batch
N_CELLS <- 5000 # per sample
N_MARKERS <- 24
MARKERS <- sprintf("M%02d", seq_len(N_MARKERS))
POP_PROPS <- c(
  P1 = 0.25,
  P2 = 0.18,
  P3 = 0.15,
  P4 = 0.12,
  P5 = 0.10,
  P6 = 0.08,
  P7 = 0.07,
  P8 = 0.05
)
POPS <- names(POP_PROPS)
POP_SD <- 0.3 # within-population SD per marker
DONOR_SD <- 0.1 # donor offset per population and marker
# Batch-effect spreads: quantiles (25%, 50%, 75%) over markers of the SD of the global shift (sd_a), of the
# cell-state-specific shift on top (sd_c) and of the log within-state scale ratio (sd_s). Each marker draws its
# own SDs from a log-normal fitted to these quantiles.
CALIBRATION <- list(
  sd_a = c(q25 = 0.044, q50 = 0.084, q75 = 0.163),
  sd_c = c(q25 = 0.090, q50 = 0.100, q75 = 0.137),
  sd_s = c(q25 = 0.218, q50 = 0.236, q75 = 0.326) # includes node sampling noise: an upper estimate
)
EFFECT_POP <- "P3"
EFFECT_MARKER <- "M05"
EFFECT_SIZE <- 0.4 # cases minus controls
CASES <- list(
  balanced = rep(6, N_BATCHES),
  confounded = round(seq(1, 11, length.out = N_BATCHES))
) # cases per batch, of N_STUDY
GRID <- c(8, 8) # all three methods
# CytoNorm.train defaults (CytoNorm 2.0.10) except the grid: nCells 1e6, nClus 10, scale FALSE.
CYN_FLOWSOM <- list(
  nCells = 1e6,
  xdim = GRID[[1]],
  ydim = GRID[[2]],
  nClus = 10,
  scale = FALSE
)
OUT <- "vignettes/compare-data"

args <- commandArgs(trailingOnly = TRUE)
stopifnot("calibration filled" = !anyNA(CALIBRATION$sd_s))
scenarios <- expand.grid(
  som_seed = SOM_SEEDS,
  sim_seed = SIM_SEEDS,
  design = c("balanced", "confounded"),
  batch_effect = c("global", "calibrated"),
  stringsAsFactors = FALSE
)
# Reproduction of the first release, which used one seed for the simulation and the methods
repro <- unique(transform(
  scenarios[scenarios$sim_seed %in% 1:3, ],
  som_seed = sim_seed
))
if (length(args) >= 4 && args[[4]] == "repro") {
  scenarios <- repro
}

# Simulation -------------------------------------------------------------------
# Population means are drawn per seed. Each study donor has its own population mix and small marker offsets;
# cases add EFFECT_SIZE to EFFECT_MARKER in EFFECT_POP. Bridge aliquots come from one reference donor. Batch
# effect: "global" = one shift per batch and marker; "calibrated" = that shift plus a shift and a log-scale change
# per batch, population and marker. The anchor has no batch effect.

draw_sd <- function(q, n) {
  exp(rnorm(n, log(q[["q50"]]), (log(q[["q75"]]) - log(q[["q25"]])) / 1.349))
}
simulate <- function(seed, batch_effect, design) {
  set.seed(seed)
  n_pop <- length(POPS)
  means <- matrix(
    runif(n_pop * N_MARKERS, 0.5, 4),
    n_pop,
    N_MARKERS,
    dimnames = list(POPS, MARKERS)
  )
  sd_a <- draw_sd(CALIBRATION$sd_a, N_MARKERS)
  sd_c <- draw_sd(CALIBRATION$sd_c, N_MARKERS)
  sd_s <- draw_sd(CALIBRATION$sd_s, N_MARKERS)
  effect <- lapply(setNames(BATCHES, BATCHES), function(b) {
    shift <- matrix(rnorm(N_MARKERS, 0, sd_a), n_pop, N_MARKERS, byrow = TRUE)
    log_scale <- matrix(0, n_pop, N_MARKERS)
    if (batch_effect == "calibrated") {
      shift <- shift +
        matrix(
          rnorm(n_pop * N_MARKERS, 0, rep(sd_c, each = n_pop)),
          n_pop,
          N_MARKERS
        )
      log_scale <- matrix(
        rnorm(n_pop * N_MARKERS, 0, rep(sd_s, each = n_pop)),
        n_pop,
        N_MARKERS
      )
    }
    if (b == BATCHES[[1]]) {
      shift[] <- 0
      log_scale[] <- 0
    }
    list(shift = shift, scale = exp(log_scale))
  })
  ref_offset <- matrix(rnorm(n_pop * N_MARKERS, 0, DONOR_SD), n_pop, N_MARKERS)
  samples <- list()
  cells <- list()
  truths <- list() # simulated values without batch effect, same row order as cells
  for (bi in seq_along(BATCHES)) {
    b <- BATCHES[[bi]]
    n_case <- CASES[[design]][[bi]]
    ids <- c(
      sprintf("%s_ref%d", b, seq_len(N_REF)),
      sprintf("%s_s%02d", b, seq_len(N_STUDY))
    )
    samples[[b]] <- data.frame(
      sample_id = ids,
      batch = b,
      reference_id = c(rep("R", N_REF), rep(NA, N_STUDY)),
      aliquot = c(sprintf("ref%d", seq_len(N_REF)), rep(NA, N_STUDY)),
      condition = c(
        rep(NA, N_REF),
        rep(c("case", "control"), c(n_case, N_STUDY - n_case))
      )
    )
    for (i in seq_along(ids)) {
      is_ref <- i <= N_REF
      props <- if (is_ref) {
        POP_PROPS
      } else {
        p <- POP_PROPS * exp(rnorm(n_pop, 0, 0.2))
        p / sum(p)
      }
      offset <- if (is_ref) {
        ref_offset
      } else {
        matrix(rnorm(n_pop * N_MARKERS, 0, DONOR_SD), n_pop, N_MARKERS)
      }
      k <- sample.int(n_pop, N_CELLS, replace = TRUE, prob = props)
      centre <- means[k, ] + offset[k, ]
      truth <- centre +
        matrix(rnorm(N_CELLS * N_MARKERS, 0, POP_SD), N_CELLS, N_MARKERS)
      if (identical(samples[[b]]$condition[[i]], "case")) {
        hit <- k == match(EFFECT_POP, POPS)
        truth[hit, EFFECT_MARKER] <- truth[hit, EFFECT_MARKER] + EFFECT_SIZE
      }
      observed <- centre +
        effect[[b]]$scale[k, ] * (truth - centre) +
        effect[[b]]$shift[k, ]
      colnames(observed) <- MARKERS
      cells[[length(cells) + 1]] <- data.frame(
        sample_id = ids[[i]],
        batch = b,
        pop = POPS[k],
        observed
      )
      truths[[length(truths) + 1]] <- truth
    }
  }
  cells <- do.call(rbind, cells)
  cells$cell_id <- seq_len(nrow(cells))
  list(
    samples = do.call(rbind, samples),
    cells = cells,
    truth = do.call(rbind, truths)
  )
}

# Methods ----------------------------------------------------------------------
# Every method returns a cell x marker matrix in input row order. bridgenorm and CytoNorm fit on ref1 and ref2
# (ref3 held out); cyCombine has no apply step and corrects all cells in sample, as designed.

run_methods <- function(sim, som_seed) {
  s <- sim$samples
  x <- sim$cells
  fit_samples <- s[c("sample_id", "batch", "reference_id", "aliquot")]
  held_ids <- s$sample_id[s$aliquot %in% "ref3"]
  fit_samples$reference_id[fit_samples$sample_id %in% held_ids] <- NA
  fit_ids <- s$sample_id[s$aliquot %in% c("ref1", "ref2")]
  is_fit_ref <- x$sample_id %in% fit_ids
  population <- rep("all", nrow(x))
  transform <- list(name = "asinh", cofactor = 5)
  all_x <- x[c("sample_id", MARKERS)]
  out <- list(uncorrected = as.matrix(x[MARKERS]))

  paired <- bn_fit_paired(
    x[is_fit_ref, c("sample_id", MARKERS)],
    fit_samples,
    population = population[is_fit_ref],
    markers = MARKERS,
    transform = transform,
    anchor = BATCHES[[1]]
  )
  out$paired <- as.matrix(bn_apply(paired, all_x, fit_samples, population)$x[
    MARKERS
  ])

  model <- bn_fit(
    x = x[is_fit_ref, c("sample_id", MARKERS)],
    samples = fit_samples,
    population = population[is_fit_ref],
    markers = MARKERS,
    transform = transform,
    anchor = BATCHES[[1]],
    grid = GRID,
    mapping_markers = MARKERS,
    seed = som_seed
  )
  out$bridgenorm <- as.matrix(bn_apply(model, all_x, fit_samples, population)$x[
    MARKERS
  ])
  # shrink_k = 0: the article's shrinkage sentence reads this arm.
  model_k0 <- bn_fit(
    x = x[is_fit_ref, c("sample_id", MARKERS)],
    samples = fit_samples,
    population = population[is_fit_ref],
    markers = MARKERS,
    transform = transform,
    anchor = BATCHES[[1]],
    grid = GRID,
    mapping_markers = MARKERS,
    seed = som_seed,
    shrink_k = 0
  )
  # bridgenorm_v100: the 1.0.0 defaults (shrink_target "quantile", shrink_k 300), to show the 1.1.0 change (OP22).
  out$bridgenorm_v100 <- as.matrix(bn_apply(
    bn_fit(
      x = x[is_fit_ref, c("sample_id", MARKERS)],
      samples = fit_samples,
      population = population[is_fit_ref],
      markers = MARKERS,
      transform = transform,
      anchor = BATCHES[[1]],
      grid = GRID,
      mapping_markers = MARKERS,
      seed = som_seed,
      shrink_k = 300,
      shrink_target = "quantile"
    ),
    all_x,
    fit_samples,
    population
  )$x[MARKERS])
  out$bridgenorm_k0 <- as.matrix(bn_apply(
    model_k0,
    all_x,
    fit_samples,
    population
  )$x[MARKERS])

  to_fs <- function(ids) {
    flowCore::flowSet(lapply(setNames(ids, ids), function(id) {
      flowCore::flowFrame(as.matrix(x[x$sample_id == id, MARKERS]))
    }))
  }
  # CytoNorm: the package default goal ("mean" of the batches) is the main arm; CytoNorm_anchor sets goal to the
  # anchor batch label, so every cluster maps to the anchor's quantiles, like bridgenorm. A fresh folder per
  # run: CytoNorm.train() reuses a FlowSOM model it finds in outputDir (recompute = FALSE).
  ids <- unique(x$sample_id)
  cyn_goal <- c(CytoNorm = "mean", CytoNorm_anchor = BATCHES[[1]])
  for (arm in names(cyn_goal)) {
    out_dir <- tempfile("cytonorm_")
    cyn <- CytoNorm::CytoNorm.train(
      files = to_fs(fit_ids),
      labels = s$batch[match(fit_ids, s$sample_id)],
      channels = MARKERS,
      transformList = NULL,
      outputDir = out_dir,
      FlowSOM.params = CYN_FLOWSOM,
      normParams = list(nQ = 99, goal = cyn_goal[[arm]]),
      seed = som_seed
    )
    normed <- CytoNorm::CytoNorm.normalize(
      cyn,
      files = to_fs(ids),
      labels = s$batch[match(ids, s$sample_id)],
      transformList = NULL,
      transformList.reverse = NULL,
      outputDir = out_dir,
      write = FALSE
    )
    cyn_x <- out$uncorrected
    for (k in seq_along(ids)) {
      rows <- x$sample_id == ids[[k]]
      stopifnot("CytoNorm keeps cell counts" = nrow(normed[[k]]) == sum(rows))
      cyn_x[rows, ] <- flowCore::exprs(normed[[k]])[, MARKERS]
    }
    out[[arm]] <- cyn_x
  }

  # cyCombine (SC3): its default (no covariate); with the true condition as covariate, which gives it labels the
  # bridge-based methods never see (bridge aliquots take the placeholder level "ref"). Its bridge options
  # (ref.batch with anchor) were tried and dropped: ComBat stops when a SOM node has no anchor-batch cells.
  cond <- s$condition[match(x$sample_id, s$sample_id)]
  cyc_input <- data.frame(
    x[MARKERS],
    batch = x$batch,
    sample = x$sample_id,
    id = x$cell_id,
    cond = ifelse(is.na(cond), "ref", cond)
  )
  cyc_args <- list(
    cyCombine = list(),
    cyCombine_covar = list(covar = "cond")
  )
  for (arm in names(cyc_args)) {
    cyc <- do.call(
      cyCombine::batch_correct,
      c(
        list(
          cyc_input,
          markers = MARKERS,
          xdim = GRID[[1]],
          ydim = GRID[[2]],
          seed = som_seed
        ),
        cyc_args[[arm]]
      )
    )
    out[[arm]] <- as.matrix(cyc[match(x$cell_id, cyc$id), MARKERS])
  }
  # oracle: the simulated values without batch effect, the floor of every metric
  out$oracle <- sim$truth
  out
}

# Metrics ----------------------------------------------------------------------
# Held-out EMD: 1D Wasserstein distance (mean absolute quantile difference, 99 probabilities) of the ref3
# aliquot between every batch pair, per marker; worst pair, then median over markers. Effect: per-sample median
# of EFFECT_MARKER in EFFECT_POP cells (true labels), mean of cases minus mean of controls, no batch term.
# rmse_truth: root mean squared difference to the simulated values without batch effect, over all study cells
# and markers (0 for the oracle; it also counts a shift of the whole scale, such as CytoNorm's mean goal).

w1 <- function(a, b) {
  p <- seq(0.01, 0.99, by = 0.01)
  mean(abs(
    quantile(a, p, names = FALSE, na.rm = TRUE) -
      quantile(b, p, names = FALSE, na.rm = TRUE)
  ))
}
score <- function(sim, corrected) {
  if (all(is.na(corrected))) {
    return(c(heldout_emd = NA_real_, effect = NA_real_, rmse_truth = NA_real_))
  }
  s <- sim$samples
  x <- sim$cells
  held <- x$sample_id %in% s$sample_id[s$aliquot %in% "ref3"]
  pairs <- combn(BATCHES, 2)
  emd <- sapply(MARKERS, function(m) {
    by_batch <- split(corrected[held, m], x$batch[held])
    max(apply(pairs, 2, function(pr) {
      w1(by_batch[[pr[[1]]]], by_batch[[pr[[2]]]])
    }))
  })
  study <- s[!is.na(s$condition), ]
  in_pop <- x$pop == EFFECT_POP
  # na.rm: cyCombine returns NA for a few cells per run (3 of 1,125,000 in one checked run)
  med <- sapply(study$sample_id, function(id) {
    median(corrected[in_pop & x$sample_id == id, EFFECT_MARKER], na.rm = TRUE)
  })
  effect <- mean(med[study$condition == "case"]) -
    mean(med[study$condition == "control"])
  study_rows <- x$sample_id %in% study$sample_id
  rmse <- sqrt(mean(
    (corrected[study_rows, ] - sim$truth[study_rows, ])^2,
    na.rm = TRUE
  ))
  c(heldout_emd = median(emd), effect = effect, rmse_truth = rmse)
}

# Run one task -----------------------------------------------------------------

if (args[[1]] == "run") {
  task <- as.integer(args[[2]])
  tmp <- args[[3]]
  sc <- scenarios[task, ]
  sim <- simulate(sc$sim_seed, sc$batch_effect, sc$design)
  corrected <- run_methods(sim, sc$som_seed)
  scores <- do.call(
    rbind,
    lapply(names(corrected), function(m) {
      message(m, ": ", sum(is.na(corrected[[m]])), " NA values")
      v <- score(sim, corrected[[m]])
      data.frame(
        batch_effect = sc$batch_effect,
        design = sc$design,
        sim_seed = sc$sim_seed,
        som_seed = sc$som_seed,
        method = m,
        heldout_emd = v[["heldout_emd"]],
        effect = v[["effect"]],
        effect_true = EFFECT_SIZE,
        rmse_truth = v[["rmse_truth"]]
      )
    })
  )
  densities <- NULL
  if (sc$sim_seed == 1 && sc$som_seed == SOM_SEEDS[[1]]) {
    x <- sim$cells
    shown <- BATCHES[c(1, 5, 10, 15)]
    rows <- which(
      x$sample_id %in%
        sim$samples$sample_id[sim$samples$aliquot %in% "ref3"] &
        x$pop == EFFECT_POP &
        x$batch %in% shown
    )
    rows <- unlist(lapply(split(rows, x$batch[rows]), function(i) {
      i[seq_len(min(300, length(i)))]
    }))
    densities <- do.call(
      rbind,
      lapply(names(corrected), function(m) {
        data.frame(
          batch_effect = sc$batch_effect,
          design = sc$design,
          method = m,
          batch = x$batch[rows],
          value = round(corrected[[m]][rows, EFFECT_MARKER], 3)
        )
      })
    )
  }
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  saveRDS(
    list(scores = scores, densities = densities),
    file.path(tmp, sprintf("task_%02d.rds", task))
  )
  quit(save = "no")
}

# Combine ----------------------------------------------------------------------

if (args[[1]] == "combine") {
  tmp <- args[[2]]
  if (length(args) >= 3) {
    OUT <- args[[3]]
  }
  parts <- lapply(
    sprintf("task_%02d.rds", seq_len(nrow(scenarios))),
    function(f) readRDS(file.path(tmp, f))
  )
  stopifnot("output folder is new" = !dir.exists(OUT))
  dir.create(OUT, recursive = TRUE)
  write.csv(
    do.call(rbind, lapply(parts, `[[`, "scores")),
    file.path(OUT, "scores.csv"),
    row.names = FALSE
  )
  write.csv(
    do.call(rbind, lapply(parts, `[[`, "densities")),
    file.path(OUT, "densities.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      setting = c(
        "simulation seeds",
        "SOM seeds",
        "batches",
        "bridge aliquots per batch",
        "study samples per batch",
        "cells per sample",
        "markers",
        "populations",
        "effect",
        "cases per batch (balanced)",
        "cases per batch (confounded)",
        "global shift SD per marker (median)",
        "state-specific shift SD per marker (median)",
        "log scale SD per marker (median)",
        "SOM grid (all methods)",
        "bridgenorm",
        "CytoNorm",
        "cyCombine"
      ),
      value = c(
        paste(SIM_SEEDS, collapse = ", "),
        paste(SOM_SEEDS, collapse = ", "),
        N_BATCHES,
        sprintf("%d (2 fit, 1 held out)", N_REF),
        N_STUDY,
        N_CELLS,
        N_MARKERS,
        length(POPS),
        sprintf(
          "+%.1f on %s in population %s",
          EFFECT_SIZE,
          EFFECT_MARKER,
          EFFECT_POP
        ),
        paste(CASES$balanced, collapse = ", "),
        paste(CASES$confounded, collapse = ", "),
        CALIBRATION$sd_a[["q50"]],
        CALIBRATION$sd_c[["q50"]],
        CALIBRATION$sd_s[["q50"]],
        paste(GRID, collapse = "x"),
        sprintf(
          "bridgenorm %s, defaults, all markers mapped",
          packageVersion("bridgenorm")
        ),
        sprintf(
          "CytoNorm %s, FlowSOM %s; defaults except the grid; CytoNorm_anchor sets normParams goal to %s",
          packageVersion("CytoNorm"),
          paste(
            names(CYN_FLOWSOM),
            vapply(CYN_FLOWSOM, format, ""),
            sep = " = ",
            collapse = ", "
          ),
          BATCHES[[1]]
        ),
        sprintf(
          paste(
            "cyCombine %s, defaults (8x8); run twice: no covariate, and the true condition as covariate"
          ),
          packageVersion("cyCombine")
        )
      )
    ),
    file.path(OUT, "settings.csv"),
    row.names = FALSE
  )
}
