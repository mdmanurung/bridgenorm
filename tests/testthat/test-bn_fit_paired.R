# bn_fit_paired(): effects against a hand computation, abstention and integration.

paired_sim_fit <- function(s, meta = s$meta, unbridged = "stop") {
  x <- s$x[s$x$sample_id %in% meta$sample_id, ]
  bn_fit_paired(
    x,
    meta,
    population = rep("T", nrow(x)),
    markers = "m1",
    transform = list(name = "asinh", cofactor = 5),
    anchor = "s1_r1",
    min_cells = 10,
    unbridged = unbridged
  )
}

test_that("paired effects equal the mean over donors of aliquot-mean differences to the anchor", {
  s <- bn_sim()
  eff <- paired_sim_fit(s)$effects
  # Hand computation: aliquot means, donor means per batch, difference to the anchor, mean over donors.
  am <- tapply(s$x$m1, s$x$sample_id, mean)
  meta <- s$meta
  meta$am <- am[meta$sample_id]
  dm <- tapply(meta$am, list(meta$donor, meta$batch), mean)
  diff <- dm - dm[, "s1_r1"]
  want <- colMeans(diff)
  want_se <- apply(diff, 2, stats::sd) / sqrt(nrow(diff))
  got <- eff[match(names(want), eff$batch), ]
  expect_equal(got$estimate, unname(want), tolerance = 1e-12)
  expect_equal(got$se_donor[-1], unname(want_se[-1]), tolerance = 1e-12)
  expect_equal(got$support, c("anchor", rep("paired_reference", 5)))
  expect_equal(paired_sim_fit(s)$method, "paired_shift")
})

test_that("a batch sharing no reference with the anchor is abstained or stops", {
  s <- bn_sim()
  keep <- !(s$meta$batch == "s1_r1" & s$meta$donor == "D3") &
    !(s$meta$batch == "s3_r2" & s$meta$donor != "D3")
  meta <- s$meta[keep, ]
  expect_false("s3_r2" %in% paired_sim_fit(s, meta, "abstain")$effects$batch)
  expect_error(paired_sim_fit(s, meta, "stop"))
})

test_that("integrating a paired model subtracts its effects", {
  s <- bn_sim()
  m <- paired_sim_fit(s)
  res <- bn_integrate(
    m,
    s$x,
    s$meta,
    rep("T", nrow(s$x)),
    transform = m$transform,
    unsupported = "error"
  )
  b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)]
  shift <- m$effects$estimate[match(b, m$effects$batch)]
  expect_equal(res$x$m1, s$x$m1 - shift, tolerance = 1e-10)
})
