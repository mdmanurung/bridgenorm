# goal_batches: target = "virtual" restricted to a set of goal batches (CytoNorm's goal distribution).

paired_fit <- function(d, ...) {
  bn_fit_paired(
    d$x,
    d$samples,
    rep("T", nrow(d$x)),
    markers = c("S", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1",
    ...
  )
}

test_that("a fit without goal_batches keeps the 1.0.0 digest", {
  d <- guard_sim()
  fixture <- bn_load(test_path("fixtures", "model-1.0.0.rds"))
  expect_identical(guard_fit(d)$digest, fixture$digest)
})

test_that("goal_batches = every batch gives the same correction as the all-batch virtual target", {
  d <- guard_sim()
  pop <- rep("T", nrow(d$x))
  all <- guard_fit(d, target = "virtual")
  goal <- guard_fit(d, target = "virtual", goal_batches = c("b1", "b2", "b3"))
  expect_equal(
    bn_apply(goal, d$x, d$samples, pop)$x$M,
    bn_apply(all, d$x, d$samples, pop)$x$M,
    tolerance = 1e-12
  )
})

test_that("the paired goal target is the equal-weight mean of the goal batches' summaries", {
  d <- guard_sim()
  goal <- paired_fit(d, target = "virtual", goal_batches = c("b2", "b3"))
  # With one donor per batch, the b2 and b3 anchor targets are those batches' donor summaries.
  t2 <- bn_fit_paired(
    d$x,
    d$samples,
    rep("T", nrow(d$x)),
    markers = c("S", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b2"
  )$targets
  t3 <- bn_fit_paired(
    d$x,
    d$samples,
    rep("T", nrow(d$x)),
    markers = c("S", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b3"
  )$targets
  expect_equal(
    unname(goal$targets$target),
    (t2$target + t3$target) / 2,
    tolerance = 1e-12
  )
})

test_that("a batch outside the goal does not move the goal target, and the anchor may lie outside it", {
  d <- guard_sim()
  a <- paired_fit(d, target = "virtual", goal_batches = c("b2", "b3"))
  b1 <- d$x$sample_id %in% d$samples$sample_id[d$samples$batch == "b1"]
  d$x$M[b1] <- d$x$M[b1] + 1
  b <- paired_fit(d, target = "virtual", goal_batches = c("b2", "b3"))
  expect_identical(a$targets$target, b$targets$target)
  anchor_effect <- b$effects$estimate[
    b$effects$batch == "b1" & b$effects$marker == "M"
  ]
  expect_gt(abs(anchor_effect), 0.5)
  q <- guard_fit(d, target = "virtual", goal_batches = c("b2", "b3"))
  res <- bn_apply(q, d$x, d$samples, rep("T", nrow(d$x)))
  expect_true(all(is.finite(res$x$M)))
})

test_that("goal_batches is checked", {
  d <- guard_sim()
  expect_error(paired_fit(d, goal_batches = c("b2", "b3")), "needs target")
  expect_error(
    paired_fit(d, target = "virtual", goal_batches = "b2"),
    "two or more"
  )
  expect_error(
    paired_fit(d, target = "virtual", goal_batches = c("b2", "b2")),
    "two or more"
  )
  expect_error(
    paired_fit(d, target = "virtual", goal_batches = c("b2", NA)),
    "two or more"
  )
  expect_error(
    paired_fit(d, target = "virtual", goal_batches = c("b2", "b9")),
    "not in 'samples'.*b9"
  )
})

test_that("the stored goal does not depend on the order it is given in", {
  d <- guard_sim()
  a <- guard_fit(d, target = "virtual", goal_batches = c("b3", "b2"))
  b <- guard_fit(d, target = "virtual", goal_batches = c("b2", "b3"))
  expect_identical(a$settings$goal_batches, c("b2", "b3"))
  expect_identical(a$digest, b$digest)
})

test_that("a goal without a summary the anchor has stops instead of giving NaN", {
  d <- guard_sim()
  pop <- rep("T", nrow(d$x))
  # Anchor reference cells of state A get their own label, absent from the goal batches.
  in_b1_ref <- d$x$sample_id %in% c("ref1_b1", "ref2_b1")
  pop[in_b1_ref & d$x$S < 2.5] <- "Z"
  expect_error(
    bn_fit_paired(
      d$x,
      d$samples,
      pop,
      markers = c("S", "M"),
      transform = list(name = "asinh", cofactor = 5),
      anchor = "b1",
      target = "virtual",
      goal_batches = c("b2", "b3")
    ),
    "goal batches have no reference summary"
  )
})
