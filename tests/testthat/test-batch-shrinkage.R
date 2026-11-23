# Adaptive shrinkage across batches (roadmap R2; plans/2026-11-20_r2-adaptive-batch-shrinkage.md).

# Independent oracle for given tau2 and centre.
eb_given <- function(y, s2, tau2, mu) {
  w <- tau2 / (tau2 + s2)
  z2 <- (y - mu)^2 / (s2 + tau2)
  a <- w + (1 - w) * z2 / (1 + z2)
  list(standard = mu + w * (y - mu), adaptive = mu + a * (y - mu))
}

test_that(".bn_moderate with a given tau2 applies the formula", {
  y <- c(0.1, 0.5, -0.2)
  s2 <- c(0.01, 0.02, 0.04)
  r <- bridgenorm:::.bn_moderate(y, s2, "adaptive", "zero", tau2 = 0.05)
  expect_equal(r$estimate, eb_given(y, s2, 0.05, 0)$adaptive)
  expect_equal(r$tau2, 0.05)
})

test_that("well variance: pooled excess of replicate variance over cell variance", {
  summ <- tibble::tibble(
    reference_id = "D",
    aliquot = c("a1", "a2", "a3", "b1", "b2", "c1"),
    batch = c("b1", "b1", "b1", "b2", "b2", "b3"),
    population = "T",
    marker = "m1",
    value = c(1.0, 1.2, 1.4, 2.0, 2.4, 3.0),
    n_cells = 100L,
    var_mean = 0.001
  )
  # Replicate SS: b1 0.08 (df 2), b2 0.08 (df 1): pooled 0.16 / 3; minus mean cell var 0.001.
  got <- bridgenorm:::.bn_well_variance(summ)
  expect_equal(got$well_var, 0.16 / 3 - 0.001)
})

# One donor (D1), 2 aliquots per batch, one marker; batch s3_r1 gets a large true extra shift.
one_donor_sim <- function(extra = 1.5) {
  s <- bn_sim()
  s <- subset_sim(s, s$meta$donor == "D1")
  b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)]
  s$x$m1 <- s$x$m1 + ifelse(b == "s3_r1", extra, 0)
  s
}

fit_shrink <- function(
  s,
  method,
  centre = "study",
  group = "study",
  tau2 = "profile"
) {
  fit_sim(
    s,
    anchor_effect = "fixed",
    variance_prior = "boundary_avoiding",
    residual_model = "known_plus_extra",
    batch_shrinkage = method,
    batch_shrinkage_centre = centre,
    batch_shrinkage_group = group,
    batch_shrinkage_tau2 = tau2
  )$effects
}

test_that("batch_shrinkage is required and 'blup' leaves the 0.3.0 effects", {
  s <- one_donor_sim()
  expect_error(
    fit_sim(s, batch_shrinkage = NULL, batch_shrinkage_centre = NULL),
    "batch_shrinkage"
  )
  blup <- fit_shrink(s, "blup")
  expect_false("estimate_blup" %in% names(blup))
})

test_that("adaptive shrinkage across batches equals the hand computation on paired contrasts", {
  s <- one_donor_sim()
  eff <- fit_shrink(s, "adaptive", centre = "zero", group = "all")
  # Hand: paired contrasts and their variances from the cells.
  am <- tapply(s$x$m1, s$x$sample_id, mean)
  vm <- tapply(s$x$m1, s$x$sample_id, stats::var) /
    tapply(s$x$m1, s$x$sample_id, length)
  meta <- s$meta
  meta$am <- am[meta$sample_id]
  meta$vm <- vm[meta$sample_id]
  dm <- tapply(meta$am, meta$batch, mean)
  y <- dm - dm[["s1_r1"]]
  cell <- tapply(meta$vm, meta$batch, sum) / 4
  rep_ss <- sum(tapply(meta$am, meta$batch, function(v) sum((v - mean(v))^2)))
  well <- max(rep_ss / (nrow(meta) - length(dm)) - mean(meta$vm), 0)
  nb <- names(y) != "s1_r1"
  s2 <- cell[nb] + cell[["s1_r1"]] + well * (1 / 2 + 1 / 2)
  want <- bridgenorm:::.bn_moderate(
    unname(y[nb]),
    unname(s2),
    "adaptive",
    "zero"
  )$estimate
  got <- eff$estimate[match(names(y)[nb], eff$batch)]
  expect_equal(got, as.vector(want), tolerance = 1e-8)
  expect_equal(eff$estimate[eff$batch == "s1_r1"], 0)
  expect_true(all(is.finite(eff$estimate_blup)))
})

test_that("an outlying batch keeps more of its shift under adaptive than standard EB", {
  s <- one_donor_sim(extra = 3)
  std <- fit_shrink(s, "standard_eb", centre = "zero", group = "all")
  ad <- fit_shrink(s, "adaptive", centre = "zero", group = "all")
  expect_gte(
    ad$estimate[ad$batch == "s3_r1"],
    std$estimate[std$batch == "s3_r1"]
  )
})

test_that("REML tau2 uses the fitted batch variance", {
  s <- one_donor_sim()
  eff <- fit_shrink(
    s,
    "standard_eb",
    centre = "study",
    group = "study",
    tau2 = "reml"
  )
  expect_true(all(is.finite(eff$estimate)))
  # Studies s2 and s3 have two non-anchor batches each; s1 has one, so its fitted centre is its own value.
  shrunk <- eff$batch %in% c("s2_r1", "s2_r2", "s3_r1", "s3_r2")
  expect_true(all(is.finite(eff$batch_shrinkage_tau2[shrunk])))
})

# Review fixes for 0.5.1 (plans/2026-11-21_heb-0.5.1-fixes.md) -----------------

test_that("a study with one non-anchor batch is labelled, not silently left unshrunk", {
  s <- one_donor_sim()
  eff <- fit_shrink(s, "standard_eb", centre = "study", group = "study")
  # s1 holds the anchor s1_r1, so s1_r2 is alone in its group; s2 and s3 have two batches each.
  expect_equal(
    eff$batch_shrinkage_status[eff$batch == "s1_r2"],
    "single_batch_group"
  )
  expect_equal(eff$batch_shrinkage_weight[eff$batch == "s1_r2"], 1)
  expect_true(all(
    eff$batch_shrinkage_status[
      eff$batch %in% c("s2_r1", "s2_r2", "s3_r1", "s3_r2")
    ] ==
      "ok"
  ))
})

test_that("REML tau2 equals the lme4 batch variance", {
  skip_if_not_installed("lme4")
  s <- bn_sim()
  eff <- fit_sim(
    s,
    batch_shrinkage = "standard_eb",
    batch_shrinkage_tau2 = "reml"
  )$effects
  d <- s$meta
  d$value <- vapply(
    d$sample_id,
    function(id) mean(s$x$m1[s$x$sample_id == id]),
    numeric(1)
  )
  d <- merge(d, s$hierarchy, by = "batch")
  fit <- lme4::lmer(
    value ~ 0 + donor + (1 | study) + (1 | batch),
    data = d,
    REML = TRUE,
    control = lme4::lmerControl(
      check.conv.grad = lme4::.makeCC("ignore", tol = 1e-3)
    )
  )
  want <- as.data.frame(lme4::VarCorr(fit))
  want <- want$vcov[want$grp == "batch"]
  got <- unique(stats::na.omit(eff$batch_shrinkage_tau2))
  expect_equal(got, want, tolerance = 1e-3)
})

test_that("shrinkage variances use donors and unequal aliquot counts as written by hand", {
  # Two donors; batch s2_r1 keeps one aliquot per donor, every other batch two.
  s <- bn_sim()
  s <- subset_sim(
    s,
    s$meta$donor %in%
      c("D1", "D2") &
      !(s$meta$batch == "s2_r1" & grepl("_2$", s$meta$sample_id))
  )
  eff <- fit_sim(
    s,
    batch_shrinkage = "standard_eb",
    batch_shrinkage_centre = "zero",
    batch_shrinkage_group = "all",
    batch_shrinkage_tau2 = "reml"
  )$effects
  meta <- s$meta
  meta$am <- vapply(
    meta$sample_id,
    function(id) mean(s$x$m1[s$x$sample_id == id]),
    numeric(1)
  )
  meta$vm <- vapply(
    meta$sample_id,
    function(id) {
      v <- s$x$m1[s$x$sample_id == id]
      stats::var(v) / length(v)
    },
    numeric(1)
  )
  # Donor means per batch, then the paired contrast to the anchor averaged over the 2 donors.
  dm <- tapply(meta$am, list(meta$batch, meta$donor), mean)
  y <- rowMeans(dm - matrix(dm["s1_r1", ], nrow(dm), 2, byrow = TRUE))
  # Cell variance of a batch's donor mean: sum of aliquot variances / n_aliquots^2, averaged over donors / 2.
  cv <- tapply(meta$vm, list(meta$batch, meta$donor), function(v) {
    sum(v) / length(v)^2
  })
  cell <- rowSums(cv) / 4
  rep_ss <- sum(
    tapply(meta$am, list(meta$batch, meta$donor), function(v) {
      sum((v - mean(v))^2)
    }),
    na.rm = TRUE
  )
  rep_df <- nrow(meta) - length(dm)
  well <- max(rep_ss / rep_df - mean(meta$vm), 0)
  k <- c(table(meta$batch)) / 2
  nb <- names(y) != "s1_r1"
  s2 <- cell[nb] +
    cell[["s1_r1"]] +
    well * (1 / k[names(y)[nb]] + 1 / k[["s1_r1"]]) / 2
  tau2 <- unique(stats::na.omit(eff$batch_shrinkage_tau2))
  want <- eb_given(y[nb], s2, tau2, 0)$standard
  got <- eff$estimate[match(names(y)[nb], eff$batch)]
  expect_equal(unname(got), unname(want), tolerance = 1e-8)
})
