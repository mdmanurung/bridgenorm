# Shrinkage options of bn_fit(): shrink_target, shrink_reweight, quantile_maps.

test_that("the defaults are shrink_target = \"value\" and shrink_k = 1000, and are recorded", {
  d <- guard_sim()
  base <- guard_fit(d)
  same <- guard_fit(
    d,
    shrink_k = 1000,
    shrink_target = "value",
    shrink_reweight = FALSE,
    quantile_maps = TRUE
  )
  expect_identical(base$settings$shrink_target, "value")
  expect_identical(base$settings$shrink_k, 1000)
  expect_null(base$settings$shrink_reweight)
  expect_null(base$settings$quantile_maps)
  expect_identical(same$digest, base$digest)
  expect_identical(same$quantiles, base$quantiles)
})

test_that("the options are recorded only when set, and bad values are refused", {
  d <- guard_sim()
  m <- guard_fit(
    d,
    shrink_target = "value",
    shrink_reweight = TRUE,
    quantile_maps = FALSE
  )
  expect_identical(m$settings$shrink_target, "value")
  expect_true(m$settings$shrink_reweight)
  expect_false(m$settings$quantile_maps)
  expect_error(guard_fit(d, shrink_target = "other"), "shrink_target")
  expect_error(guard_fit(d, shrink_reweight = NA), "shrink_reweight")
  expect_error(guard_fit(d, quantile_maps = "no"), "quantile_maps")
})

test_that("with shrink_k = 0 the shrinkage options change nothing", {
  d <- guard_sim()
  base <- guard_fit(d, shrink_k = 0)
  value <- guard_fit(d, shrink_k = 0, shrink_target = "value")
  rew <- guard_fit(d, shrink_k = 0, shrink_reweight = TRUE)
  expect_identical(value$quantiles, base$quantiles)
  expect_identical(rew$quantiles, base$quantiles)
})

test_that("the value target and reweighting give finite, monotone tables for anchor and virtual targets", {
  d <- guard_sim()
  for (target in c("anchor", "virtual")) {
    m <- guard_fit(
      d,
      target = target,
      shrink_k = 1e6,
      shrink_target = "value",
      shrink_reweight = TRUE
    )
    for (q in m$quantiles) {
      expect_true(all(is.finite(q)))
      expect_false(any(apply(q, 2, is.unsorted)))
    }
    out <- bn_apply(m, d$x, d$samples, rep("T", nrow(d$x)))$x
    expect_true(all(is.finite(out$M)))
  }
})

test_that("quantile_maps = FALSE keeps no tables and corrects with the node shifts", {
  d <- guard_sim()
  pop <- rep("T", nrow(d$x))
  m <- guard_fit(d, quantile_maps = FALSE)
  expect_length(m$quantiles, 0)
  out <- bn_apply(m, d$x, d$samples, pop)
  expect_false(any(out$actions$reason_code == "node_quantile"))
  expect_true(all(is.finite(out$x$M)))
  mv <- guard_fit(d, quantile_maps = FALSE, target = "virtual")
  expect_length(mv$quantiles, 0)
})

test_that("the value target of a pure scale map is that scale map", {
  probs <- seq(0.1, 0.9, length.out = 9)
  target <- cbind(M = qnorm(probs, 2, 0.5))
  batch <- 0.3 + 1.4 * target # the batch-wide map is x -> 0.3 + 1.4 x
  node_q <- cbind(M = qnorm(probs, 2.2, 0.1)) # a node inside the batch-wide range
  got <- bridgenorm:::.bn_value_target(node_q, target, batch)
  expect_equal(got[, "M"], 0.3 + 1.4 * node_q[, "M"], tolerance = 0.02)
  # the quantile target adds the offset at the same probability instead
  expect_gt(
    max(abs(got[, "M"] - (node_q[, "M"] + (batch - target)[, "M"]))),
    0.05
  )
})

test_that("weighted tables with equal weights match the unweighted quantiles", {
  set.seed(2)
  xm <- cbind(M = rnorm(2000))
  probs <- c(0.1, 0.5, 0.9)
  tab <- bridgenorm:::.bn_weighted_table(xm, rep(1, 2000), probs)
  expect_equal(
    tab[, "M"],
    unname(quantile(xm[, "M"], probs)),
    tolerance = 0.05
  )
  heavy <- bridgenorm:::.bn_weighted_table(xm, ifelse(xm[, "M"] > 0, 3, 1), 0.5)
  expect_gt(heavy[1, "M"], 0)
})

test_that("a weighted table with fewer than two weighted cells returns the fallback", {
  xm <- cbind(M = c(1, 2, 3))
  fb <- matrix(9, 2, 1, dimnames = list(NULL, "M"))
  expect_identical(
    bridgenorm:::.bn_weighted_table(xm, c(0, 0, NA), c(0.25, 0.75), fb),
    fb
  )
  expect_error(
    bridgenorm:::.bn_weighted_table(xm, c(0, 0, 1), c(0.25, 0.75)),
    "positive weight"
  )
})

test_that("with reweighting, batches that differ only in cell-state mix are left in place", {
  # No batch effect; the batches differ only in their share of state A (0.8, 0.5, 0.2). The batch-wide curves
  # and the target curve must be compared at one cell-state mix, or the shrinkage moves aligned cells.
  set.seed(3)
  batches <- c("b1", "b2", "b3")
  share_a <- c(b1 = 0.8, b2 = 0.5, b3 = 0.2)
  samples <- data.frame(
    sample_id = c(outer(c("ref1", "ref2", "s1"), batches, paste, sep = "_")),
    batch = rep(batches, each = 3),
    reference_id = rep(c("donor", "donor", NA), 3),
    aliquot = rep(c("ref1", "ref2", NA), 3)
  )
  x <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
    n_a <- round(600 * share_a[[samples$batch[i]]])
    data.frame(
      sample_id = samples$sample_id[i],
      S = c(stats::rnorm(n_a, 1, 0.2), stats::rnorm(600 - n_a, 4, 0.2)),
      M = c(stats::rnorm(n_a, 2, 0.3), stats::rnorm(600 - n_a, 0.5, 0.3))
    )
  }))
  pop <- rep("T", nrow(x))
  for (target in c("anchor", "virtual")) {
    m <- bn_fit(x, samples, pop, markers = c("S", "M"), transform = list(name = "asinh", cofactor = 5),
      anchor = "b1", grid = c(2, 1), mapping_markers = "S", seed = 1, train_cells = 300, rlen = 50,
      target = target, shrink_k = 1e6, shrink_target = "value", shrink_reweight = TRUE)
    out <- bn_apply(m, x, samples, pop)$x
    expect_lt(median(abs(out$M - x$M)), 0.05)
  }
})
