# Purpose: Simulated comparison of bridgenorm with CytoNorm and cyCombine for the "Compare with CytoNorm and
#          cyCombine" article. The design and the batch-effect sizes are set to realistic values (numbers below);
#          no data are used. Four scenarios (batch effect global only, or
#          global plus cell-state-specific shift and scale) x (condition balanced within batch, or confounded with
#          batch), 3 simulation seeds each.
# Run:     Rscript tools/sim-compare.R run <task> <tmp_dir>   one scenario x seed (task 1-12), e.g. as a job array
#          Rscript tools/sim-compare.R combine <tmp_dir>      merge the task results into vignettes/compare-data
#          (package root; needs bridgenorm, kohonen, CytoNorm, flowCore, cyCombine)
# Inputs:  none (simulated)
# Outputs (vignettes/compare-data/):
#   scores.csv        batch_effect, design, seed, method, heldout_emd (held-out bridge aliquot: worst batch pair,
#                     median over markers), effect (estimated condition effect), effect_true
#   densities.csv     batch_effect, design, method, batch, value: the effect marker in effect-population cells of
#                     the held-out bridge aliquot, seed 1, 300 cells for 4 of the batches
#   settings.csv      setting, value
#   (package versions are recorded in settings.csv; no sessionInfo, which would carry local paths)

suppressPackageStartupMessages({
  library(bridgenorm)
})

# Parameters -------------------------------------------------------------------
SEEDS <- 1:3
N_BATCHES <- 15
BATCHES <- sprintf("b%02d", seq_len(N_BATCHES)) # b01 = anchor, no batch effect
N_REF <- 3 # bridge aliquots per batch; ref1 and ref2 fit, ref3 held out
N_STUDY <- 12 # study samples (donors) per batch
N_CELLS <- 5000 # per sample
N_MARKERS <- 24
MARKERS <- sprintf("M%02d", seq_len(N_MARKERS))
POP_PROPS <- c(P1 = 0.25, P2 = 0.18, P3 = 0.15, P4 = 0.12, P5 = 0.10, P6 = 0.08, P7 = 0.07, P8 = 0.05)
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
CYN_FLOWSOM <- list(nCells = 1e6, xdim = GRID[[1]], ydim = GRID[[2]], nClus = 10, scale = FALSE)
OUT <- "vignettes/compare-data"

args <- commandArgs(trailingOnly = TRUE)
stopifnot("calibration filled" = !anyNA(CALIBRATION$sd_s))
scenarios <- expand.grid(
  seed = SEEDS,
  design = c("balanced", "confounded"),
  batch_effect = c("global", "calibrated"),
  stringsAsFactors = FALSE
)

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
  means <- matrix(runif(n_pop * N_MARKERS, 0.5, 4), n_pop, N_MARKERS, dimnames = list(POPS, MARKERS))
  sd_a <- draw_sd(CALIBRATION$sd_a, N_MARKERS)
  sd_c <- draw_sd(CALIBRATION$sd_c, N_MARKERS)
  sd_s <- draw_sd(CALIBRATION$sd_s, N_MARKERS)
  effect <- lapply(setNames(BATCHES, BATCHES), function(b) {
    shift <- matrix(rnorm(N_MARKERS, 0, sd_a), n_pop, N_MARKERS, byrow = TRUE)
    log_scale <- matrix(0, n_pop, N_MARKERS)
    if (batch_effect == "calibrated") {
      shift <- shift + matrix(rnorm(n_pop * N_MARKERS, 0, rep(sd_c, each = n_pop)), n_pop, N_MARKERS)
      log_scale <- matrix(rnorm(n_pop * N_MARKERS, 0, rep(sd_s, each = n_pop)), n_pop, N_MARKERS)
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
  for (bi in seq_along(BATCHES)) {
    b <- BATCHES[[bi]]
    n_case <- CASES[[design]][[bi]]
    ids <- c(sprintf("%s_ref%d", b, seq_len(N_REF)), sprintf("%s_s%02d", b, seq_len(N_STUDY)))
    samples[[b]] <- data.frame(
      sample_id = ids,
      batch = b,
      reference_id = c(rep("R", N_REF), rep(NA, N_STUDY)),
      aliquot = c(sprintf("ref%d", seq_len(N_REF)), rep(NA, N_STUDY)),
      condition = c(rep(NA, N_REF), rep(c("case", "control"), c(n_case, N_STUDY - n_case)))
    )
    for (i in seq_along(ids)) {
      is_ref <- i <= N_REF
      props <- if (is_ref) {
        POP_PROPS
      } else {
        p <- POP_PROPS * exp(rnorm(n_pop, 0, 0.2))
        p / sum(p)
      }
      offset <- if (is_ref) ref_offset else matrix(rnorm(n_pop * N_MARKERS, 0, DONOR_SD), n_pop, N_MARKERS)
      k <- sample.int(n_pop, N_CELLS, replace = TRUE, prob = props)
      centre <- means[k, ] + offset[k, ]
      truth <- centre + matrix(rnorm(N_CELLS * N_MARKERS, 0, POP_SD), N_CELLS, N_MARKERS)
      if (identical(samples[[b]]$condition[[i]], "case")) {
        hit <- k == match(EFFECT_POP, POPS)
        truth[hit, EFFECT_MARKER] <- truth[hit, EFFECT_MARKER] + EFFECT_SIZE
      }
      observed <- centre + effect[[b]]$scale[k, ] * (truth - centre) + effect[[b]]$shift[k, ]
      colnames(observed) <- MARKERS
      cells[[length(cells) + 1]] <- data.frame(sample_id = ids[[i]], batch = b, pop = POPS[k], observed)
    }
  }
  cells <- do.call(rbind, cells)
  cells$cell_id <- seq_len(nrow(cells))
  list(samples = do.call(rbind, samples), cells = cells)
}

# Methods ----------------------------------------------------------------------
# Every method returns a cell x marker matrix in input row order. bridgenorm and CytoNorm fit on ref1 and ref2
# (ref3 held out); cyCombine has no apply step and corrects all cells in sample, as designed.

run_methods <- function(sim, seed) {
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
  out$paired <- as.matrix(bn_apply(paired, all_x, fit_samples, population)$x[MARKERS])

  model <- bn_fit(
    x = x[is_fit_ref, c("sample_id", MARKERS)],
    samples = fit_samples,
    population = population[is_fit_ref],
    markers = MARKERS,
    transform = transform,
    anchor = BATCHES[[1]],
    grid = GRID,
    mapping_markers = MARKERS,
    seed = seed
  )
  out$bridgenorm <- as.matrix(bn_apply(model, all_x, fit_samples, population)$x[MARKERS])
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
    seed = seed,
    shrink_k = 0
  )
  out$bridgenorm_k0 <- as.matrix(bn_apply(model_k0, all_x, fit_samples, population)$x[MARKERS])

  to_fs <- function(ids) {
    flowCore::flowSet(lapply(setNames(ids, ids), function(id) {
      flowCore::flowFrame(as.matrix(x[x$sample_id == id, MARKERS]))
    }))
  }
  # A fresh folder per run: CytoNorm.train() reuses a FlowSOM model it finds in outputDir (recompute = FALSE).
  out_dir <- tempfile("cytonorm_")
  cyn <- CytoNorm::CytoNorm.train(
    files = to_fs(fit_ids),
    labels = s$batch[match(fit_ids, s$sample_id)],
    channels = MARKERS,
    transformList = NULL,
    outputDir = out_dir,
    FlowSOM.params = CYN_FLOWSOM,
    seed = seed
  )
  ids <- unique(x$sample_id)
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
  out$CytoNorm <- cyn_x

  # cyCombine twice (SC3): its default (no covariate), and with the true condition as covariate, which gives it
  # labels the bridge-based methods never see; bridge aliquots take the placeholder level "ref".
  cond <- s$condition[match(x$sample_id, s$sample_id)]
  cyc_input <- data.frame(
    x[MARKERS],
    batch = x$batch,
    sample = x$sample_id,
    id = x$cell_id,
    cond = ifelse(is.na(cond), "ref", cond)
  )
  for (arm in c("cyCombine", "cyCombine_covar")) {
    cyc <- cyCombine::batch_correct(
      cyc_input,
      markers = MARKERS,
      xdim = GRID[[1]],
      ydim = GRID[[2]],
      seed = seed,
      covar = if (arm == "cyCombine_covar") "cond" else NULL
    )
    cyc <- cyc[match(x$cell_id, cyc$id), ]
    out[[arm]] <- as.matrix(cyc[MARKERS])
  }
  out
}

# Metrics ----------------------------------------------------------------------
# Held-out EMD: 1D Wasserstein distance (mean absolute quantile difference, 99 probabilities) of the ref3
# aliquot between every batch pair, per marker; worst pair, then median over markers. Effect: per-sample median
# of EFFECT_MARKER in EFFECT_POP cells (true labels), mean of cases minus mean of controls, no batch term.

w1 <- function(a, b) {
  p <- seq(0.01, 0.99, by = 0.01)
  mean(abs(quantile(a, p, names = FALSE) - quantile(b, p, names = FALSE)))
}
score <- function(sim, corrected) {
  s <- sim$samples
  x <- sim$cells
  held <- x$sample_id %in% s$sample_id[s$aliquot %in% "ref3"]
  pairs <- combn(BATCHES, 2)
  emd <- sapply(MARKERS, function(m) {
    by_batch <- split(corrected[held, m], x$batch[held])
    max(apply(pairs, 2, function(pr) w1(by_batch[[pr[[1]]]], by_batch[[pr[[2]]]])))
  })
  study <- s[!is.na(s$condition), ]
  in_pop <- x$pop == EFFECT_POP
  # na.rm: cyCombine returns NA for a few cells per run (3 of 1,125,000 in one checked run)
  med <- sapply(study$sample_id, function(id) {
    median(corrected[in_pop & x$sample_id == id, EFFECT_MARKER], na.rm = TRUE)
  })
  effect <- mean(med[study$condition == "case"]) - mean(med[study$condition == "control"])
  c(heldout_emd = median(emd), effect = effect)
}

# Run one task -----------------------------------------------------------------

if (args[[1]] == "run") {
  task <- as.integer(args[[2]])
  tmp <- args[[3]]
  sc <- scenarios[task, ]
  sim <- simulate(sc$seed, sc$batch_effect, sc$design)
  corrected <- run_methods(sim, sc$seed)
  scores <- do.call(rbind, lapply(names(corrected), function(m) {
    v <- score(sim, corrected[[m]])
    data.frame(
      batch_effect = sc$batch_effect, design = sc$design, seed = sc$seed, method = m,
      heldout_emd = v[["heldout_emd"]], effect = v[["effect"]], effect_true = EFFECT_SIZE
    )
  }))
  densities <- NULL
  if (sc$seed == 1) {
    x <- sim$cells
    shown <- BATCHES[c(1, 5, 10, 15)]
    rows <- which(x$sample_id %in% sim$samples$sample_id[sim$samples$aliquot %in% "ref3"] &
      x$pop == EFFECT_POP & x$batch %in% shown)
    rows <- unlist(lapply(split(rows, x$batch[rows]), function(i) i[seq_len(min(300, length(i)))]))
    densities <- do.call(rbind, lapply(names(corrected), function(m) {
      data.frame(
        batch_effect = sc$batch_effect, design = sc$design, method = m,
        batch = x$batch[rows], value = round(corrected[[m]][rows, EFFECT_MARKER], 3)
      )
    }))
  }
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(scores = scores, densities = densities), file.path(tmp, sprintf("task_%02d.rds", task)))
  quit(save = "no")
}

# Combine ----------------------------------------------------------------------

if (args[[1]] == "combine") {
  tmp <- args[[2]]
  parts <- lapply(sprintf("task_%02d.rds", seq_len(nrow(scenarios))), function(f) readRDS(file.path(tmp, f)))
  stopifnot("output folder is new" = !dir.exists(OUT))
  dir.create(OUT, recursive = TRUE)
  write.csv(do.call(rbind, lapply(parts, `[[`, "scores")), file.path(OUT, "scores.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(parts, `[[`, "densities")), file.path(OUT, "densities.csv"), row.names = FALSE)
  write.csv(data.frame(
    setting = c(
      "seeds", "batches", "bridge aliquots per batch", "study samples per batch", "cells per sample",
      "markers", "populations", "effect", "cases per batch (balanced)", "cases per batch (confounded)",
      "global shift SD per marker (median)", "state-specific shift SD per marker (median)",
      "log scale SD per marker (median)", "SOM grid (all methods)", "bridgenorm", "CytoNorm", "cyCombine"
    ),
    value = c(
      paste(SEEDS, collapse = ", "), N_BATCHES, sprintf("%d (2 fit, 1 held out)", N_REF), N_STUDY, N_CELLS,
      N_MARKERS, length(POPS), sprintf("+%.1f on %s in population %s", EFFECT_SIZE, EFFECT_MARKER, EFFECT_POP),
      paste(CASES$balanced, collapse = ", "), paste(CASES$confounded, collapse = ", "),
      CALIBRATION$sd_a[["q50"]], CALIBRATION$sd_c[["q50"]], CALIBRATION$sd_s[["q50"]],
      paste(GRID, collapse = "x"),
      sprintf("bridgenorm %s, defaults, all markers mapped", packageVersion("bridgenorm")),
      sprintf(
        "CytoNorm %s, FlowSOM %s; defaults except the grid",
        packageVersion("CytoNorm"),
        paste(names(CYN_FLOWSOM), vapply(CYN_FLOWSOM, format, ""), sep = " = ", collapse = ", ")
      ),
      sprintf(
        "cyCombine %s, defaults (8x8); run twice: no covariate, and the true condition as covariate",
        packageVersion("cyCombine")
      )
    )
  ), file.path(OUT, "settings.csv"), row.names = FALSE)
}
