# Purpose: Simulated comparison of bridgenorm with CytoNorm and cyCombine for the "Compare with CytoNorm and
#          cyCombine" article. Four scenarios (batch effect global or state-specific x condition balanced within
#          batch or confounded with batch), 3 simulation seeds each. Writes small result tables that the article
#          reads; the article does not rerun the methods.
# Run:     Rscript tools/sim-compare.R   (package root; needs bridgenorm, kohonen, CytoNorm, flowCore, cyCombine)
# Inputs:  none (simulated)
# Outputs (vignettes/compare-data/):
#   scores.csv        scenario, seed, method, heldout_emd (worst batch pair, median over markers, held-out bridge
#                     aliquot), effect (estimated condition effect), effect_true
#   densities.csv     scenario, method, batch, value: marker M4 in population C cells of the held-out bridge
#                     aliquot, seed 1, 400 cells per batch
#   settings.csv      setting, value: the simulation and method settings used
#   (package versions are recorded in settings.csv; no sessionInfo, which would carry local paths)

suppressPackageStartupMessages({
  library(bridgenorm)
})

# Parameters -------------------------------------------------------------------
SEEDS <- 1:3
N_CELLS <- 1500 # per sample
N_STUDY <- 8 # study samples (donors) per batch
BATCHES <- c("b1", "b2", "b3", "b4") # b1 = anchor, no batch effect
MARKERS <- paste0("M", 1:6)
POPS <- c("A", "B", "C", "D")
POP_MEANS <- rbind(
  A = c(1, 3, 1, 1, 2, 1),
  B = c(3, 1, 1, 2, 1, 1),
  C = c(1, 1, 3, 1.5, 1, 2),
  D = c(3, 3, 2, 1, 1, 3)
) # rows = populations, columns = M1..M6 (asinh-like scale)
POP_SD <- 0.25
POP_PROPS <- c(A = 0.35, B = 0.25, C = 0.25, D = 0.15)
EFFECT_POP <- "C"
EFFECT_MARKER <- "M4"
EFFECT_SIZE <- 0.6 # cases minus controls, M4 in population C
CASES_PER_BATCH <- list(balanced = c(4, 4, 4, 4), confounded = c(1, 2, 6, 7))
GRID <- c(5, 5) # matched across methods (user choice)
CYN_NCLUS <- 10 # CytoNorm metaclusters on the 5x5 map; its default 30 exceeds 25 nodes
OUT <- "vignettes/compare-data"

stopifnot("output folder is new" = !dir.exists(OUT))
dir.create(OUT, recursive = TRUE)

# Simulation -------------------------------------------------------------------
# Every sample: cells from 4 populations. Study donors vary in population mix and carry small marker offsets;
# cases add EFFECT_SIZE to EFFECT_MARKER in EFFECT_POP. Bridge aliquots ref1 (fit) and ref2 (held out) come from
# one reference donor. Batch effect: "global" = one shift per batch and marker; "state" = shift and scale per
# batch, population and marker. b1 has no batch effect.

simulate <- function(seed, batch_effect, design) {
  set.seed(seed)
  effect <- lapply(setNames(BATCHES, BATCHES), function(b) {
    if (b == "b1") {
      return(list(shift = matrix(0, 4, 6, dimnames = list(POPS, MARKERS)), scale = matrix(1, 4, 6)))
    }
    if (batch_effect == "global") {
      s <- runif(6, -0.6, 0.6)
      list(shift = matrix(s, 4, 6, byrow = TRUE, dimnames = list(POPS, MARKERS)), scale = matrix(1, 4, 6))
    } else {
      list(
        shift = matrix(runif(24, -0.6, 0.6), 4, 6, dimnames = list(POPS, MARKERS)),
        scale = matrix(runif(24, 0.7, 1.4), 4, 6)
      )
    }
  })
  ref_offset <- matrix(rnorm(24, 0, 0.1), 4, 6)
  samples <- list()
  cells <- list()
  for (bi in seq_along(BATCHES)) {
    b <- BATCHES[[bi]]
    n_case <- CASES_PER_BATCH[[design]][[bi]]
    ids <- c(paste0(b, "_ref", 1:2), sprintf("%s_s%02d", b, seq_len(N_STUDY)))
    samples[[b]] <- data.frame(
      sample_id = ids,
      batch = b,
      reference_id = c("R", "R", rep(NA, N_STUDY)),
      aliquot = c("ref1", "ref2", rep(NA, N_STUDY)),
      condition = c(NA, NA, rep(c("case", "control"), c(n_case, N_STUDY - n_case)))
    )
    for (i in seq_along(ids)) {
      is_ref <- i <= 2
      props <- if (is_ref) POP_PROPS else {
        p <- POP_PROPS * exp(rnorm(4, 0, 0.2))
        p / sum(p)
      }
      offset <- if (is_ref) ref_offset else matrix(rnorm(24, 0, 0.1), 4, 6)
      pop <- sample(POPS, N_CELLS, replace = TRUE, prob = props)
      truth <- POP_MEANS[pop, ] + offset[match(pop, POPS), ] +
        matrix(rnorm(N_CELLS * 6, 0, POP_SD), N_CELLS, 6)
      if (identical(samples[[b]]$condition[[i]], "case")) {
        truth[pop == EFFECT_POP, EFFECT_MARKER == MARKERS] <-
          truth[pop == EFFECT_POP, EFFECT_MARKER == MARKERS] + EFFECT_SIZE
      }
      k <- match(pop, POPS)
      centre <- POP_MEANS[pop, ]
      observed <- centre + effect[[b]]$scale[k, ] * (truth - centre) + effect[[b]]$shift[k, ]
      colnames(observed) <- MARKERS
      cells[[length(cells) + 1]] <- data.frame(sample_id = ids[[i]], batch = b, pop = pop, observed)
    }
  }
  cells <- do.call(rbind, cells)
  cells$cell_id <- seq_len(nrow(cells))
  list(samples = do.call(rbind, samples), cells = cells)
}

# Methods ----------------------------------------------------------------------
# Every method returns a cell x marker matrix in input row order. bridgenorm and CytoNorm fit on ref1 aliquots
# only (ref2 is held out); cyCombine has no apply step and corrects all cells in sample (its design).

run_methods <- function(sim, seed) {
  s <- sim$samples
  x <- sim$cells
  fit_samples <- s
  fit_samples$reference_id[fit_samples$aliquot %in% "ref2"] <- NA # ref2 held out
  is_fit_ref <- x$sample_id %in% s$sample_id[s$aliquot %in% "ref1"]
  population <- rep("all", nrow(x))
  transform <- list(name = "asinh", cofactor = 5)
  out <- list(uncorrected = as.matrix(x[MARKERS]))

  paired <- bn_fit_paired(
    x[is_fit_ref, c("sample_id", MARKERS)],
    fit_samples[1:4],
    population = population[is_fit_ref],
    markers = MARKERS,
    transform = transform,
    anchor = "b1"
  )
  out$paired <- as.matrix(bn_apply(paired, x[c("sample_id", MARKERS)], fit_samples[1:4], population)$x[MARKERS])

  model <- bn_fit(
    x = x[is_fit_ref, c("sample_id", MARKERS)],
    samples = fit_samples[1:4],
    population = population[is_fit_ref],
    markers = MARKERS,
    transform = transform,
    anchor = "b1",
    grid = GRID,
    mapping_markers = MARKERS,
    seed = seed
  )
  out$bridgenorm <- as.matrix(bn_apply(model, x[c("sample_id", MARKERS)], fit_samples[1:4], population)$x[MARKERS])

  to_fs <- function(ids) {
    flowCore::flowSet(lapply(setNames(ids, ids), function(id) {
      flowCore::flowFrame(as.matrix(x[x$sample_id == id, MARKERS]))
    }))
  }
  ref_ids <- s$sample_id[s$aliquot %in% "ref1"]
  # A fresh folder per run: CytoNorm.train() reuses a FlowSOM model it finds in outputDir (recompute = FALSE).
  out_dir <- tempfile("cytonorm_")
  cyn <- CytoNorm::CytoNorm.train(
    files = to_fs(ref_ids),
    labels = s$batch[match(ref_ids, s$sample_id)],
    channels = MARKERS,
    transformList = NULL,
    outputDir = out_dir,
    FlowSOM.params = list(nCells = 1e6, xdim = GRID[[1]], ydim = GRID[[2]], nClus = CYN_NCLUS, scale = FALSE),
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

  cyc <- cyCombine::batch_correct(
    data.frame(x[MARKERS], batch = x$batch, sample = x$sample_id, id = x$cell_id),
    markers = MARKERS,
    xdim = GRID[[1]],
    ydim = GRID[[2]],
    seed = seed
  )
  cyc <- cyc[match(x$cell_id, cyc$id), ]
  out$cyCombine <- as.matrix(cyc[MARKERS])
  out
}

# Metrics ----------------------------------------------------------------------
# Held-out EMD: 1D Wasserstein distance (mean absolute quantile difference, 99 probabilities) of the ref2 aliquot
# between every batch pair, per marker; worst pair, then median over markers. Effect: per-sample median of
# EFFECT_MARKER in EFFECT_POP cells (true labels), mean of cases minus mean of controls, no batch term.

w1 <- function(a, b) {
  p <- seq(0.01, 0.99, by = 0.01)
  mean(abs(quantile(a, p, names = FALSE) - quantile(b, p, names = FALSE)))
}
score <- function(sim, corrected) {
  s <- sim$samples
  x <- sim$cells
  held <- x$sample_id %in% s$sample_id[s$aliquot %in% "ref2"]
  pairs <- combn(BATCHES, 2)
  emd <- sapply(MARKERS, function(m) {
    max(apply(pairs, 2, function(pr) {
      w1(corrected[held & x$batch == pr[[1]], m], corrected[held & x$batch == pr[[2]], m])
    }))
  })
  study <- s[!is.na(s$condition), ]
  in_pop <- x$pop == EFFECT_POP
  med <- sapply(study$sample_id, function(id) median(corrected[in_pop & x$sample_id == id, EFFECT_MARKER]))
  effect <- mean(med[study$condition == "case"]) - mean(med[study$condition == "control"])
  c(heldout_emd = median(emd), effect = effect)
}

# Run --------------------------------------------------------------------------

scenarios <- expand.grid(batch_effect = c("global", "state"), design = c("balanced", "confounded"),
  stringsAsFactors = FALSE)
scores <- list()
densities <- list()
for (r in seq_len(nrow(scenarios))) {
  sc <- scenarios[r, ]
  label <- paste(sc$batch_effect, sc$design, sep = "_")
  for (seed in SEEDS) {
    message("Scenario ", label, ", seed ", seed)
    sim <- simulate(seed, sc$batch_effect, sc$design)
    corrected <- run_methods(sim, seed)
    for (m in names(corrected)) {
      sc_m <- score(sim, corrected[[m]])
      scores[[length(scores) + 1]] <- data.frame(
        batch_effect = sc$batch_effect, design = sc$design, seed = seed, method = m,
        heldout_emd = sc_m[["heldout_emd"]], effect = sc_m[["effect"]], effect_true = EFFECT_SIZE
      )
      if (seed == 1) {
        x <- sim$cells
        rows <- which(x$sample_id %in% sim$samples$sample_id[sim$samples$aliquot %in% "ref2"] &
          x$pop == EFFECT_POP)
        rows <- unlist(lapply(split(rows, x$batch[rows]), function(i) i[seq_len(min(400, length(i)))]))
        densities[[length(densities) + 1]] <- data.frame(
          batch_effect = sc$batch_effect, design = sc$design, method = m,
          batch = x$batch[rows], value = round(corrected[[m]][rows, EFFECT_MARKER], 3)
        )
      }
    }
  }
}
write.csv(do.call(rbind, scores), file.path(OUT, "scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, densities), file.path(OUT, "densities.csv"), row.names = FALSE)
write.csv(data.frame(
  setting = c("seeds", "cells per sample", "study samples per batch", "batches", "markers", "populations",
    "effect", "cases per batch (balanced)", "cases per batch (confounded)", "SOM grid (all methods)",
    "bridgenorm", "CytoNorm", "cyCombine"),
  value = c(paste(SEEDS, collapse = ", "), N_CELLS, N_STUDY, length(BATCHES), length(MARKERS), length(POPS),
    sprintf("+%.1f on %s in population %s", EFFECT_SIZE, EFFECT_MARKER, EFFECT_POP),
    paste(CASES_PER_BATCH$balanced, collapse = ", "), paste(CASES_PER_BATCH$confounded, collapse = ", "),
    paste(GRID, collapse = "x"),
    sprintf("bridgenorm %s, defaults, all markers mapped", packageVersion("bridgenorm")),
    sprintf("CytoNorm %s, defaults except xdim/ydim 5, nClus %d", packageVersion("CytoNorm"), CYN_NCLUS),
    sprintf("cyCombine %s, defaults except xdim/ydim 5", packageVersion("cyCombine")))
), file.path(OUT, "settings.csv"), row.names = FALSE)
