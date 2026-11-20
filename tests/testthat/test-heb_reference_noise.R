# heb_reference_noise() (roadmap R7): report-only noise per batch, against hand computations.

noise_sim <- function(s, meta = s$meta) {
  x <- s$x[s$x$sample_id %in% meta$sample_id, ]
  heb_reference_noise(
    x,
    meta,
    population = rep("T", nrow(x)),
    markers = "m1",
    transform = list(name = "asinh", cofactor = 5),
    min_cells = 10
  )
}

test_that("noise columns equal hand computations", {
  s <- heb_sim()
  got <- noise_sim(s)
  meta <- s$meta
  meta$am <- tapply(s$x$m1, s$x$sample_id, mean)[meta$sample_id]
  meta$sd <- tapply(s$x$m1, s$x$sample_id, stats::sd)[meta$sample_id]
  dm <- tapply(meta$am, list(meta$donor, meta$batch), mean)
  for (b in unique(meta$batch)) {
    mb <- meta[meta$batch == b, ]
    # Pooled within-donor SD of aliquot means.
    dev <- mb$am - stats::ave(mb$am, mb$donor)
    rep_sd <- sqrt(sum(dev^2) / (nrow(mb) - length(unique(mb$donor))))
    # Leave-one-batch-out z per donor, averaged over donors.
    others <- dm[, colnames(dm) != b, drop = FALSE]
    z <- mean((dm[, b] - rowMeans(others)) / apply(others, 1, stats::sd))
    row <- got[got$batch == b, ]
    expect_equal(row$n_references, 3L)
    expect_equal(row$n_aliquots, 6L)
    expect_equal(row$replicate_sd, rep_sd, tolerance = 1e-12)
    expect_equal(row$cell_sd, mean(mb$sd), tolerance = 1e-12)
    expect_equal(row$loo_z, z, tolerance = 1e-12)
  }
})

test_that("one aliquot per batch gives no replicate SD", {
  s <- heb_sim()
  meta <- s$meta[s$meta$donor == "D1" & endsWith(s$meta$sample_id, "_1"), ]
  got <- noise_sim(s, meta)
  expect_true(all(is.na(got$replicate_sd)))
  expect_true(all(is.finite(got$loo_z)))
})
