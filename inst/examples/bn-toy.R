# Two batches, each with two bridge wells of one reference donor and one
# study sample. Marker M is shifted by 0.5 in batch b2; Lin separates two
# cell states.
set.seed(1)
samples <- data.frame(
  sample_id = c("r1_b1", "r2_b1", "s_b1", "r1_b2", "r2_b2", "s_b2"),
  batch = rep(c("b1", "b2"), each = 3),
  reference_id = rep(c("donor", "donor", NA), 2),
  aliquot = rep(c("a1", "a2", NA), 2)
)
x <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
  shift <- if (samples$batch[i] == "b2") 0.5 else 0
  data.frame(
    sample_id = samples$sample_id[i],
    Lin = c(rnorm(200, 1, 0.2), rnorm(200, 4, 0.2)),
    M = c(rnorm(200, 2, 0.3), rnorm(200, 1, 0.3)) + shift
  )
}))
