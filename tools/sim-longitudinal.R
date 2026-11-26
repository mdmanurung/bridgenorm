# Purpose: Longitudinal simulation for the "Compare with CytoNorm and cyCombine" article. The calibrated simulator of
#          tools/sim-compare.R with a staggered clinical design: each batch holds t1 of a new 4-subject cohort, t2 of
#          the previous cohort and t3 of the one before, so time and batch overlap within every subject. Known time
#          effect on EFFECT_MARKER in EFFECT_POP (t2 +0.2, t3 +0.4), none on the other markers. Two batch effects
#          (global; global plus cell-state-specific) x 5 seeds.
# Run:     Rscript tools/sim-longitudinal.R run <task> <tmp_dir>   one batch effect x seed (task 1-10), e.g. as a job array
#          Rscript tools/sim-longitudinal.R combine <tmp_dir>      merge into vignettes/compare-data/longitudinal.csv
#          (package root; needs bridgenorm, kohonen, CytoNorm, flowCore)
# Inputs:  tools/sim-compare.R (design constants, calibration, batch-effect model)
# Outputs (vignettes/compare-data/):
#   longitudinal.csv  batch_effect, seed, method, heldout_emd, ratio_t2, ratio_t3 (estimated / true time effect),
#                     null_drift (median |mean t3 - t1 change| over no-change markers), noise_sd (median SD across
#                     subjects of the t3 - t1 change, no-change markers), t_t2, t_t3 (paired t of the time effect)

suppressPackageStartupMessages({
  library(bridgenorm)
})

# Parameters -------------------------------------------------------------------
LG_SEEDS <- 1:5
N_COHORTS <- 13 # cohorts enter in batches 1-13; t3 of the last lands in batch 15
COHORT_SIZE <- 4
TIME_EFFECT <- c(t1 = 0, t2 = 0.2, t3 = 0.4)
WITHIN_SD <- 0.03 # visit-to-visit noise per subject, population and marker (small against DONOR_SD 0.1)
CYN_GRIDS <- list(CytoNorm_default = c(15, 15), CytoNorm_8x8 = c(8, 8)) # CytoNorm default grid and the shared 8x8

# Design constants and calibration from tools/sim-compare.R (everything except its command-line blocks).
for (e in parse("tools/sim-compare.R")) {
  txt <- paste(deparse(e), collapse = "")
  if (grepl("^args <- commandArgs", txt) || grepl("^if \\(args", txt)) next
  eval(e, globalenv())
}
LG_OUT <- "vignettes/compare-data/longitudinal.csv"
tasks <- expand.grid(seed = LG_SEEDS, batch_effect = c("global", "calibrated"), stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# Run one task -----------------------------------------------------------------

if (args[[1]] == "run") {
  task <- as.integer(args[[2]])
  tmp <- args[[3]]
  tk <- tasks[task, ]

  # Simulation -------------------------------------------------------------------
  # Batch effects exactly as simulate() in tools/sim-compare.R. Each subject keeps one population mix and marker
  # offset across visits, plus WITHIN_SD noise per visit.

  set.seed(tk$seed)
  n_pop <- length(POPS)
  means <- matrix(runif(n_pop * N_MARKERS, 0.5, 4), n_pop, N_MARKERS, dimnames = list(POPS, MARKERS))
  sd_a <- draw_sd(CALIBRATION$sd_a, N_MARKERS)
  sd_c <- draw_sd(CALIBRATION$sd_c, N_MARKERS)
  sd_s <- draw_sd(CALIBRATION$sd_s, N_MARKERS)
  effect <- lapply(setNames(BATCHES, BATCHES), function(b) {
    shift <- matrix(rnorm(N_MARKERS, 0, sd_a), n_pop, N_MARKERS, byrow = TRUE)
    log_scale <- matrix(0, n_pop, N_MARKERS)
    if (tk$batch_effect == "calibrated") {
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
  subjects <- lapply(seq_len(N_COHORTS * COHORT_SIZE), function(j) {
    p <- POP_PROPS * exp(rnorm(n_pop, 0, 0.2))
    list(props = p / sum(p), offset = matrix(rnorm(n_pop * N_MARKERS, 0, DONOR_SD), n_pop, N_MARKERS))
  })

  draw_cells <- function(props, offset, b, shift_time) {
    k <- sample.int(n_pop, N_CELLS, replace = TRUE, prob = props)
    centre <- means[k, ] + offset[k, ]
    truth <- centre + matrix(rnorm(N_CELLS * N_MARKERS, 0, POP_SD), N_CELLS, N_MARKERS)
    hit <- k == match(EFFECT_POP, POPS)
    truth[hit, EFFECT_MARKER] <- truth[hit, EFFECT_MARKER] + shift_time
    observed <- centre + effect[[b]]$scale[k, ] * (truth - centre) + effect[[b]]$shift[k, ]
    colnames(observed) <- MARKERS
    data.frame(pop = POPS[k], observed)
  }
  samples <- list()
  cells <- list()
  for (bi in seq_along(BATCHES)) {
    b <- BATCHES[[bi]]
    for (r in seq_len(N_REF)) {
      id <- sprintf("%s_ref%d", b, r)
      samples[[id]] <- data.frame(sample_id = id, batch = b, reference_id = "R", aliquot = sprintf("ref%d", r),
        subject = NA, timepoint = NA)
      cells[[id]] <- data.frame(sample_id = id, batch = b, draw_cells(POP_PROPS, ref_offset, b, 0))
    }
    for (visit in 1:3) {
      cohort <- bi - visit + 1
      if (cohort < 1 || cohort > N_COHORTS) next
      for (j in (cohort - 1) * COHORT_SIZE + seq_len(COHORT_SIZE)) {
        id <- sprintf("%s_subj%02d_t%d", b, j, visit)
        sj <- subjects[[j]]
        visit_offset <- sj$offset + matrix(rnorm(n_pop * N_MARKERS, 0, WITHIN_SD), n_pop, N_MARKERS)
        samples[[id]] <- data.frame(sample_id = id, batch = b, reference_id = NA, aliquot = NA,
          subject = sprintf("subj%02d", j), timepoint = sprintf("t%d", visit))
        cells[[id]] <- data.frame(sample_id = id, batch = b, draw_cells(sj$props, visit_offset, b, TIME_EFFECT[[visit]]))
      }
    }
  }
  s <- do.call(rbind, samples)
  x <- do.call(rbind, cells)
  rownames(s) <- NULL
  rownames(x) <- NULL

  # Methods ----------------------------------------------------------------------
  # As tools/sim-compare.R: fit on ref1 and ref2, ref3 held out, all cells corrected. bridgenorm_k0 (shrink_k = 0)
  # is a diagnostic arm for the text, not the figure.

  fit_samples <- s[c("sample_id", "batch", "reference_id", "aliquot")]
  fit_samples$reference_id[s$aliquot %in% "ref3"] <- NA
  fit_ids <- s$sample_id[s$aliquot %in% c("ref1", "ref2")]
  is_fit_ref <- x$sample_id %in% fit_ids
  population <- rep("all", nrow(x))
  transform <- list(name = "asinh", cofactor = 5)
  all_x <- x[c("sample_id", MARKERS)]
  out <- list(uncorrected = as.matrix(x[MARKERS]))
  paired <- bn_fit_paired(x[is_fit_ref, c("sample_id", MARKERS)], fit_samples, population = population[is_fit_ref],
    markers = MARKERS, transform = transform, anchor = BATCHES[[1]])
  out$paired <- as.matrix(bn_apply(paired, all_x, fit_samples, population)$x[MARKERS])
  model <- bn_fit(x = x[is_fit_ref, c("sample_id", MARKERS)], samples = fit_samples, population = population[is_fit_ref],
    markers = MARKERS, transform = transform, anchor = BATCHES[[1]], grid = GRID, mapping_markers = MARKERS,
    seed = tk$seed)
  out$bridgenorm <- as.matrix(bn_apply(model, all_x, fit_samples, population)$x[MARKERS])

  out$bridgenorm_k0 <- as.matrix(bn_apply(
    bn_fit(x = x[is_fit_ref, c("sample_id", MARKERS)], samples = fit_samples, population = population[is_fit_ref],
      markers = MARKERS, transform = transform, anchor = BATCHES[[1]], grid = GRID, mapping_markers = MARKERS,
      seed = tk$seed, shrink_k = 0),
    all_x, fit_samples, population
  )$x[MARKERS])

  to_fs <- function(ids) {
    flowCore::flowSet(lapply(setNames(ids, ids), function(id) {
      flowCore::flowFrame(as.matrix(x[x$sample_id == id, MARKERS]))
    }))
  }
  ids <- unique(x$sample_id)
  for (arm in names(CYN_GRIDS)) {
    out_dir <- tempfile("cytonorm_")
    cyn <- CytoNorm::CytoNorm.train(
      files = to_fs(fit_ids),
      labels = s$batch[match(fit_ids, s$sample_id)],
      channels = MARKERS,
      transformList = NULL,
      outputDir = out_dir,
      FlowSOM.params = list(nCells = 1e6, xdim = CYN_GRIDS[[arm]][[1]], ydim = CYN_GRIDS[[arm]][[2]], nClus = 10,
        scale = FALSE),
      seed = tk$seed
    )
    normed <- CytoNorm::CytoNorm.normalize(cyn, files = to_fs(ids), labels = s$batch[match(ids, s$sample_id)],
      transformList = NULL, transformList.reverse = NULL, outputDir = out_dir, write = FALSE)
    cyn_x <- out$uncorrected
    for (k in seq_along(ids)) {
      rows <- x$sample_id == ids[[k]]
      stopifnot("CytoNorm keeps cell counts" = nrow(normed[[k]]) == sum(rows))
      cyn_x[rows, ] <- flowCore::exprs(normed[[k]])[, MARKERS]
    }
    out[[arm]] <- cyn_x
  }

  # Scores -----------------------------------------------------------------------
  # Held-out EMD as score() in tools/sim-compare.R. Per study sample, the median of every marker in EFFECT_POP cells;
  # per subject, t2 - t1 and t3 - t1; means over subjects. Ratio = mean change of EFFECT_MARKER / truth; null drift =
  # median over the other markers of |mean t3 - t1 change| (truth 0). Precision: noise_sd = median over the other
  # markers of the SD across subjects of the t3 - t1 change; t_t2, t_t3 = paired t of the effect marker's change.

  held <- x$sample_id %in% s$sample_id[s$aliquot %in% "ref3"]
  pairs <- combn(BATCHES, 2)
  study <- s[!is.na(s$subject), ]
  in_pop <- x$pop == EFFECT_POP
  scores <- do.call(rbind, lapply(names(out), function(m) {
    cx <- out[[m]]
    emd <- sapply(MARKERS, function(mk) {
      by_batch <- split(cx[held, mk], x$batch[held])
      max(apply(pairs, 2, function(pr) w1(by_batch[[pr[[1]]]], by_batch[[pr[[2]]]])))
    })
    med <- t(sapply(study$sample_id, function(id) {
      apply(cx[in_pop & x$sample_id == id, , drop = FALSE], 2, median, na.rm = TRUE)
    }))
    change <- function(tp) {
      a <- med[study$sample_id[study$timepoint == tp], , drop = FALSE]
      b0 <- med[study$sample_id[study$timepoint == "t1"], , drop = FALSE]
      sa <- study$subject[study$timepoint == tp]
      sb <- study$subject[study$timepoint == "t1"]
      a - b0[match(sa, sb), , drop = FALSE]
    }
    c2 <- change("t2")
    c3 <- change("t3")
    d2 <- colMeans(c2)
    d3 <- colMeans(c3)
    null <- setdiff(MARKERS, EFFECT_MARKER)
    paired_t <- function(v) mean(v) / (sd(v) / sqrt(length(v)))
    data.frame(
      batch_effect = tk$batch_effect, seed = tk$seed, method = m, heldout_emd = median(emd),
      ratio_t2 = d2[[EFFECT_MARKER]] / TIME_EFFECT[["t2"]], ratio_t3 = d3[[EFFECT_MARKER]] / TIME_EFFECT[["t3"]],
      null_drift = median(abs(d3[null])),
      noise_sd = median(apply(c3[, null, drop = FALSE], 2, sd)),
      t_t2 = paired_t(c2[, EFFECT_MARKER]), t_t3 = paired_t(c3[, EFFECT_MARKER])
    )
  }))
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  write.csv(scores, file.path(tmp, sprintf("longitudinal_%02d.csv", task)), row.names = FALSE)
  quit(save = "no")
}

# Combine ----------------------------------------------------------------------

if (args[[1]] == "combine") {
  tmp <- args[[2]]
  parts <- lapply(sprintf("longitudinal_%02d.csv", seq_len(nrow(tasks))), function(f) read.csv(file.path(tmp, f)))
  stopifnot("output file is new" = !file.exists(LG_OUT))
  write.csv(do.call(rbind, parts), LG_OUT, row.names = FALSE)
}
