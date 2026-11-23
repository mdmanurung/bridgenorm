# Failure paths found on the data set A design: one reference donor, one aliquot per batch.

test_that("a batch whose references all fall below min_cells counts as unbridged", {
  s <- bn_sim()
  thin <- s$x$sample_id %in% s$meta$sample_id[s$meta$batch == "s3_r2"]
  keep <- !thin |
    stats::ave(seq_along(thin), s$x$sample_id, FUN = seq_along) <= 5
  s$x <- s$x[keep, ]
  expect_error(fit_sim(s, unbridged = "stop"), "s3_r2")
  m <- fit_sim(s, unbridged = "abstain")
  expect_false("s3_r2" %in% m$effects$batch)
  expect_false(m$bridges$bridged[m$bridges$batch == "s3_r2"])
  expect_error(
    bn_fit_paired(
      s$x,
      s$meta,
      population = rep("T", nrow(s$x)),
      markers = "m1",
      transform = TF,
      anchor = "s1_r1",
      min_cells = 10,
      unbridged = "stop"
    ),
    "s3_r2"
  )
})

test_that("bn_reference_noise checks min_cells", {
  s <- bn_sim()
  expect_error(
    bn_reference_noise(
      s$x,
      s$meta,
      population = rep("T", nrow(s$x)),
      markers = "m1",
      transform = TF,
      min_cells = NA
    ),
    "min_cells"
  )
})
