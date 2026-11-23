# Failure paths found on the data set A design: one reference donor, one aliquot per batch.

# One donor, five batches in one study, one aliquot each; shifts in asinh units.
one_ref_sim <- function(shift, n = 5000, seed = 1) {
  set.seed(seed)
  b <- paste0("b", seq_along(shift))
  x <- do.call(
    rbind,
    lapply(seq_along(b), function(i) {
      data.frame(
        sample_id = paste0("REF", i),
        M1 = stats::rnorm(n, 2 + shift[[i]], 0.5)
      )
    })
  )
  list(
    x = x,
    meta = data.frame(
      sample_id = paste0("REF", seq_along(b)),
      batch = b,
      reference_id = "D"
    ),
    hierarchy = data.frame(batch = b, study = "s1")
  )
}

test_that("a singular mixed-model system stops with the population and marker named", {
  s <- one_ref_sim(c(0, 0.3, -0.2, 0.5, 0.1) * 0.01)
  expect_error(
    bn_fit(
      s$x,
      s$meta,
      population = rep("all", nrow(s$x)),
      markers = "M1",
      transform = TF,
      anchor = "b1",
      hierarchy = s$hierarchy,
      min_cells = 50,
      unbridged = "abstain",
      variance_floor = 1e-8,
      pool = "none",
      anchor_effect = "random",
      variance_prior = "boundary_avoiding",
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
    "population 'all', marker 'M1'"
  )
})

test_that("bn_fit_shared fits one donor and one population", {
  s <- one_ref_sim(c(0, 0.3, -0.2, 0.5, 0.1))
  m <- bn_fit_shared(
    s$x,
    s$meta,
    population = rep("all", nrow(s$x)),
    markers = "M1",
    transform = TF,
    anchor = "b1",
    hierarchy = s$hierarchy,
    min_cells = 50,
    variance_floor = 1e-8
  )
  expect_equal(nrow(m$effects), 5)
  expect_true(all(is.finite(m$effects$estimate)))
})

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
