# Simulations shared by the tests (ported from SpectraWeaveR test-corrector_heb*.R).

# Reference-only simulation: 3 donors measured in 6 batches (3 studies x 2 runs).
heb_sim <- function(seed = 5) {
  set.seed(seed)
  batches <- paste0("s", rep(1:3, each = 2), "_r", rep(1:2, 3))
  study_eff <- c(s1 = 0, s2 = 0.6, s3 = -0.4)
  run_eff <- stats::rnorm(6, 0, 0.3)
  names(run_eff) <- batches
  donors <- c("D1", "D2", "D3")
  donor_eff <- c(D1 = 0, D2 = 1, D3 = -1)
  rows <- list()
  for (b in batches) {
    for (d in donors) {
      for (a in 1:2) {
        mu <- donor_eff[[d]] + study_eff[[substr(b, 1, 2)]] + run_eff[[b]]
        rows[[length(rows) + 1]] <- data.frame(
          sample_id = paste(d, b, a, sep = "_"),
          batch = b,
          donor = d,
          aliquot = paste(d, b, a, sep = "_"),
          reference_id = d,
          mu = mu
        )
      }
    }
  }
  meta <- do.call(rbind, rows)
  cells <- do.call(
    rbind,
    lapply(seq_len(nrow(meta)), function(i) {
      data.frame(
        sample_id = meta$sample_id[i],
        m1 = stats::rnorm(40, meta$mu[i], 0.4)
      )
    })
  )
  list(
    x = tibble::as_tibble(cells),
    meta = meta,
    hierarchy = data.frame(batch = batches, study = substr(batches, 1, 2))
  )
}

# heb_fit() on heb_sim() with plain settings; any setting can be overridden.
fit_sim <- function(s, x = s$x, meta = s$meta, ...) {
  args <- utils::modifyList(
    list(
      x = x,
      samples = meta,
      population = rep("T", nrow(x)),
      markers = "m1",
      transform = list(name = "asinh", cofactor = 5),
      anchor = "s1_r1",
      hierarchy = s$hierarchy,
      min_cells = 10,
      unbridged = "stop",
      variance_floor = 1e-8,
      pool = "none",
      anchor_effect = "random",
      variance_prior = "none",
      residual_model = "estimated",
      scale_model = "none",
      moderation = "none",
      moderation_centre = "fitted",
      moderation_variance = "pev",
      batch_shrinkage = "blup",
      batch_shrinkage_centre = "study",
      batch_shrinkage_group = "study",
      batch_shrinkage_tau2 = "profile"
    ),
    list(...)
  )
  do.call(heb_fit, args)
}

# Keep only the samples where `keep` is TRUE.
subset_sim <- function(s, keep) {
  s$meta <- s$meta[keep, ]
  s$x <- s$x[s$x$sample_id %in% s$meta$sample_id, ]
  s
}

# 3 donors, 6 batches (3 studies x 2 runs), populations P and Q, 2 aliquots each.
# `few` = c(batch, population) cuts that cell group to 3 cells to test borrowing.
hebs_sim <- function(
  seed = 11,
  n_cells = 40,
  few = NULL,
  drop_anchor_pop = NULL
) {
  set.seed(seed)
  batches <- paste0("s", rep(1:3, each = 2), "_r", rep(1:2, 3))
  study_eff <- c(s1 = 0, s2 = 0.6, s3 = -0.4)
  run_eff <- stats::setNames(stats::rnorm(6, 0, 0.3), batches)
  dev_eff <- matrix(
    stats::rnorm(12, 0, 0.2),
    6,
    2,
    dimnames = list(batches, c("P", "Q"))
  )
  donor_eff <- outer(c(D1 = 0, D2 = 1, D3 = -1), c(P = 0, Q = 2), `+`)
  meta <- list()
  cells <- list()
  for (b in batches) {
    for (d in c("D1", "D2", "D3")) {
      for (a in 1:2) {
        for (k in c("P", "Q")) {
          sid <- paste(d, b, a, sep = "_")
          n <- if (!is.null(few) && b == few[[1]] && k == few[[2]]) {
            3
          } else {
            n_cells
          }
          if (
            !is.null(drop_anchor_pop) && b == "s1_r1" && k == drop_anchor_pop
          ) {
            n <- 3
          }
          mu <- donor_eff[d, k] +
            study_eff[[substr(b, 1, 2)]] +
            run_eff[[b]] +
            dev_eff[b, k]
          cells[[length(cells) + 1]] <- data.frame(
            sample_id = sid,
            pop = k,
            m1 = stats::rnorm(n, mu, 0.4)
          )
          if (k == "P") {
            meta[[length(meta) + 1]] <- data.frame(
              sample_id = sid,
              batch = b,
              donor = d,
              aliquot = sid,
              reference_id = d
            )
          }
        }
      }
    }
  }
  meta <- do.call(rbind, meta)
  cells <- do.call(rbind, cells)
  list(
    x = tibble::tibble(sample_id = cells$sample_id, m1 = cells$m1),
    pop = cells$pop,
    meta = meta,
    hierarchy = data.frame(batch = batches, study = substr(batches, 1, 2))
  )
}

fit_shared_sim <- function(s) {
  heb_fit_shared(
    s$x,
    s$meta,
    population = s$pop,
    markers = "m1",
    transform = list(name = "asinh", cofactor = 5),
    anchor = "s1_r1",
    hierarchy = s$hierarchy,
    min_cells = 10,
    variance_floor = 1e-8
  )
}

TF <- list(name = "asinh", cofactor = 5)

# Toy data: anchor b1 and batches b2, b3; two reference donors per batch plus one study sample in b2 and b3.
# Cells come from two well-separated clusters on the mapping marker CD3 (A near 1, B near 4).
# Marker M carries a batch shift that differs by cluster: b2 A +0.5, B -0.3; b3 A +0.2, B +0.2.
# small_b3: the b3 references have 30 cells per cluster (same cell mix, fewer cells), so with min_cells 40
# no b3 node summary survives while the b3 global summary (60 cells) does.
nodes_sim <- function(seed = 11, small_b3 = FALSE) {
  set.seed(seed)
  shift <- list(
    b1 = c(A = 0, B = 0),
    b2 = c(A = 0.5, B = -0.3),
    b3 = c(A = 0.2, B = 0.2)
  )
  base_m <- c(A = 2, B = 1)
  donor_m <- c(D1 = 0, D2 = 0.1, S = 0)
  centre_cd3 <- c(A = 1, B = 4)
  samples <- rbind(
    expand.grid(
      donor = c("D1", "D2"),
      batch = names(shift),
      stringsAsFactors = FALSE
    ),
    data.frame(donor = "S", batch = c("b2", "b3"))
  )
  samples$sample_id <- paste(samples$donor, samples$batch, sep = "_")
  samples$reference_id <- ifelse(
    samples$donor == "S",
    NA_character_,
    samples$donor
  )
  cells <- do.call(
    rbind,
    lapply(seq_len(nrow(samples)), function(i) {
      b <- samples$batch[i]
      d <- samples$donor[i]
      n <- if (small_b3 && b == "b3" && d != "S") 30 else 150
      do.call(
        rbind,
        lapply(c("A", "B"), function(k) {
          data.frame(
            sample_id = samples$sample_id[i],
            cluster = k,
            CD3 = stats::rnorm(n, centre_cd3[[k]], 0.2),
            M = stats::rnorm(
              n,
              base_m[[k]] + donor_m[[d]] + shift[[b]][[k]],
              0.1
            )
          )
        })
      )
    })
  )
  list(x = cells, samples = samples[, c("sample_id", "batch", "reference_id")])
}
