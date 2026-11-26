# Simulation shared by the 1.0.0 guard tests and the saved-model fixture (tests/testthat/fixtures).
# Two cell states (A, B) on a state marker S and a corrected marker M, three batches, two bridge wells per batch.
guard_sim <- function(
  ref_ids = c("donor", "donor"),
  n_ref = c(300, 300),
  donor_off = c(0, 0),
  study_ref = NA,
  zero_b = FALSE,
  n_b3 = 300
) {
  set.seed(1)
  centre <- c(A = 2, B = 1)
  s_pos <- c(A = 1, B = 4)
  shift <- list(
    b1 = c(A = 0, B = 0),
    b2 = c(A = 0.4, B = -0.3),
    b3 = c(A = -0.2, B = 0.3)
  )
  scale <- c(b1 = 1, b2 = 0.6, b3 = 1.4)
  samples <- data.frame(
    sample_id = c(outer(
      c("ref1", "ref2", "s1", "s2"),
      c("b1", "b2", "b3"),
      paste,
      sep = "_"
    )),
    batch = rep(c("b1", "b2", "b3"), each = 4),
    reference_id = rep(c(ref_ids, study_ref, study_ref), 3),
    aliquot = rep(c("ref1", "ref2", NA, NA), 3)
  )
  cells <- do.call(
    rbind,
    lapply(seq_len(nrow(samples)), function(i) {
      b <- samples$batch[i]
      j <- (i - 1) %% 4 + 1
      n <- if (j <= 2 && b == "b2") {
        n_ref[j]
      } else if (j <= 2 && b == "b3") {
        n_b3
      } else {
        300
      }
      off <- if (j <= 2) donor_off[j] else mean(donor_off)
      do.call(
        rbind,
        lapply(c("A", "B"), function(k) {
          truth <- stats::rnorm(n, centre[[k]] + off, 0.3)
          m <- centre[[k]] +
            scale[[b]] * (truth - centre[[k]]) +
            shift[[b]][[k]]
          if (zero_b && k == "B") {
            m <- rep(0, n)
          }
          data.frame(
            sample_id = samples$sample_id[i],
            S = stats::rnorm(n, s_pos[[k]], 0.2),
            M = m
          )
        })
      )
    })
  )
  list(x = cells, samples = samples)
}

guard_fit <- function(d, min_cells = 50) {
  bn_fit(
    d$x,
    d$samples,
    rep("T", nrow(d$x)),
    markers = c("S", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1",
    grid = c(2, 1),
    mapping_markers = "S",
    seed = 1,
    train_cells = 300,
    rlen = 50,
    min_cells = min_cells
  )
}
