# Simulations shared by the tests (ported from SpectraWeaveR test-corrector_heb*.R).

# Reference-only simulation: 3 donors measured in 6 batches (3 studies x 2 runs).
bn_sim <- function(seed = 5) {
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

# bn_fit_paired() on bn_sim(); any setting can be overridden ('hierarchy' is ignored).
fit_sim <- function(s, x = s$x, meta = s$meta, ...) {
  args <- utils::modifyList(
    list(
      x = x,
      samples = meta,
      population = rep("T", nrow(x)),
      markers = "m1",
      transform = list(name = "asinh", cofactor = 5),
      anchor = "s1_r1",
      min_cells = 10,
      unbridged = "stop"
    ),
    list(...)
  )
  args$hierarchy <- NULL
  do.call(bn_fit_paired, args)
}

# Keep only the samples where `keep` is TRUE.
subset_sim <- function(s, keep) {
  s$meta <- s$meta[keep, ]
  s$x <- s$x[s$x$sample_id %in% s$meta$sample_id, ]
  s
}




TF <- list(name = "asinh", cofactor = 5)

# Toy data: anchor b1 and batches b2, b3; two reference donors per batch plus one study sample in b2 and b3.
# Cells come from two well-separated clusters on the mapping marker Lin (A near 1, B near 4).
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
  centre_lin <- c(A = 1, B = 4)
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
            Lin = stats::rnorm(n, centre_lin[[k]], 0.2),
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
