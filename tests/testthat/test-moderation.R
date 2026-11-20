# Moderation of batch effects across markers (plans/2026-11-20_adaptive-moderation.md).

# Independent oracle: ComBat-refQL fit_adaptive_eb (R/shrink-batch-effects.R, MIT), rewritten from its
# published formulas, extended with an optional zero centre.
oracle_eb <- function(y, s2, centre = "fitted") {
  mean_at <- function(t) {
    if (centre == "zero") 0 else sum(y / (s2 + t)) / sum(1 / (s2 + t))
  }
  profile <- function(t) {
    0.5 * sum(log(s2 + t) + (y - mean_at(t))^2 / (s2 + t))
  }
  upper <- max(stats::var(y) * 100, max(s2), 1e-8)
  if (centre == "zero") {
    upper <- max(upper, max(y^2))
  }
  # Global search by a dense log grid (5000 points from 1e-14), then a local refine.
  g <- c(0, 10^seq(-14, log10(upper), length.out = 5000))
  nll <- vapply(g, profile, numeric(1))
  k <- which.min(nll)
  t <- if (k == 1) {
    0
  } else {
    exp(
      stats::optimize(
        function(lt) profile(exp(lt)),
        log(c(max(g[k - 1], g[k] / 10), g[min(k + 1, length(g))])),
        tol = 1e-12
      )$minimum
    )
  }
  if (profile(0) <= profile(t) + sqrt(.Machine$double.eps)) {
    t <- 0
  }
  mu <- mean_at(t)
  w <- t / (t + s2)
  z2 <- (y - mu)^2 / (s2 + t)
  a <- w + (1 - w) * z2 / (1 + z2)
  list(
    standard = w * y + (1 - w) * mu,
    adaptive = a * y + (1 - a) * mu,
    tau2 = t,
    centre = mu
  )
}

test_that("hand-computed cases: fitted centre", {
  r <- heb:::.heb_moderate(c(0, 2), c(2, 2), "standard_eb", "fitted")
  expect_equal(r$tau2, 0)
  expect_equal(r$estimate, c(1, 1))
  r <- heb:::.heb_moderate(c(0, 2), c(2, 2), "adaptive", "fitted")
  expect_equal(r$estimate, c(2 / 3, 4 / 3))
})

test_that("hand-computed case: zero centre", {
  r <- heb:::.heb_moderate(c(-2, 2), c(1, 1), "adaptive", "zero")
  expect_equal(r$tau2, 3, tolerance = 1e-4)
  expect_equal(r$estimate, c(-1.75, 1.75), tolerance = 1e-4)
})

test_that("moderation matches the oracle on heterogeneous variances", {
  set.seed(3)
  y <- stats::rnorm(24, 0.3, 0.4)
  s2 <- stats::runif(24, 0.01, 0.2)
  for (centre in c("fitted", "zero")) {
    o <- oracle_eb(y, s2, centre)
    expect_equal(
      heb:::.heb_moderate(y, s2, "standard_eb", centre)$estimate,
      o$standard,
      tolerance = 1e-6
    )
    expect_equal(
      heb:::.heb_moderate(y, s2, "adaptive", centre)$estimate,
      o$adaptive,
      tolerance = 1e-6
    )
  }
})

test_that("limits: precise effects stay, outliers keep their value, one marker passes through", {
  y <- c(0.1, 0.12, 0.09, 0.11, 3)
  s2 <- rep(0.01, 5)
  r <- heb:::.heb_moderate(y, s2, "adaptive", "fitted")
  expect_gt(r$weight[5], 0.9)
  r1 <- heb:::.heb_moderate(0.5, 0.01, "adaptive", "fitted")
  expect_equal(r1$estimate, 0.5)
  expect_equal(r1$weight, 1)
  expect_equal(r1$status, "too_few_markers")
})

test_that("cell variance of a contrast sums aliquot var / k^2 over batch and anchor", {
  summ <- tibble::tibble(
    reference_id = "D",
    aliquot = c("a1", "a2", "a3"),
    batch = c("b2", "b2", "b1"),
    population = "T",
    marker = "m1",
    value = c(1, 1, 0),
    n_cells = 100L,
    var_mean = c(4, 8, 1)
  )
  v <- heb:::.heb_cell_variance(summ, "b1")
  expect_equal(v$cell_var[v$batch == "b2"], (4 + 8) / 4 + 1)
})

# Two markers: m2 shares m1's batch shifts plus noise.
two_marker_sim <- function() {
  s <- heb_sim()
  s$x$m2 <- s$x$m1 + stats::rnorm(nrow(s$x), 0, 0.2)
  s
}

test_that("moderation is a required argument", {
  s <- two_marker_sim()
  args <- list(
    x = s$x,
    samples = s$meta,
    population = rep("T", nrow(s$x)),
    markers = c("m1", "m2"),
    transform = TF,
    anchor = "s1_r1",
    hierarchy = s$hierarchy,
    min_cells = 10,
    unbridged = "stop",
    variance_floor = 1e-8,
    pool = "none",
    anchor_effect = "fixed",
    variance_prior = "none",
    residual_model = "known_plus_extra",
    scale_model = "none"
  )
  expect_error(do.call(heb_fit, args), "moderation")
})

test_that("a moderated fit keeps raw estimates, anchor zero, and applies the formula per batch", {
  s <- two_marker_sim()
  fit <- function(moderation, variance = "pev") {
    fit_sim(
      s,
      markers = c("m1", "m2"),
      anchor_effect = "fixed",
      residual_model = "known_plus_extra",
      moderation = moderation,
      moderation_centre = "fitted",
      moderation_variance = variance
    )
  }
  none <- fit("none")$effects
  ad <- fit("adaptive")$effects
  expect_equal(ad$estimate_raw, none$estimate)
  expect_true(all(ad$estimate[ad$batch == "s1_r1"] == 0))
  b <- ad[ad$batch == "s2_r1", ]
  want <- oracle_eb(b$estimate_raw, b$se_donor^2)$adaptive
  expect_equal(b$estimate, want, tolerance = 1e-6)
  expect_false("estimate_raw" %in% names(none))
  ce <- fit("adaptive", "cell")$effects
  expect_true(all(is.finite(ce$moderation_weight[ce$batch != "s1_r1"])))
})

# Review fixes (Codex xhigh review, 2026-11-20) -------------------------------

test_that("a marker with zero sampling variance is left as it is and does not crash", {
  r <- heb:::.heb_moderate(
    c(0.1, 0.2, 0.3),
    c(0, 0.01, 0.01),
    "adaptive",
    "fitted"
  )
  expect_equal(r$estimate[1], 0.1)
  expect_equal(r$weight[1], 1)
  expect_true(all(is.finite(r$estimate)))
})

test_that("tau2 reaches the global minimum of the profile likelihood (brute-force grid)", {
  brute <- function(y, s2) {
    prof <- function(t) 0.5 * sum(log(s2 + t) + y^2 / (s2 + t))
    g <- c(0, 10^seq(-12, 2, length.out = 20000))
    g[which.min(vapply(g, prof, numeric(1)))]
  }
  # Two basins: many precise near-zero effects plus a few large noisy ones.
  y <- c(rep(c(-0.002, 0.002), 10), 2, -2, 1.8)
  s2 <- c(rep(1e-6, 20), 0.5, 0.5, 0.5)
  want <- brute(y, s2)
  got <- heb:::.heb_moderate(y, s2, "adaptive", "zero")$tau2
  prof <- function(t) 0.5 * sum(log(s2 + t) + y^2 / (s2 + t))
  expect_lte(prof(got), prof(want) + 1e-6)
  # Rescaling the analytic case by 1e-4 keeps tau2 = 3 x 1e-8.
  r <- heb:::.heb_moderate(c(-2, 2) * 1e-4, c(1, 1) * 1e-8, "adaptive", "zero")
  expect_equal(r$tau2, 3e-8, tolerance = 1e-3)
})

test_that("the adaptive weight keeps an outlier more than standard EB does", {
  y <- c(0.1, 0.12, 0.09, 0.11, 0.1, 0.6)
  s2 <- rep(0.01, 6)
  std <- heb:::.heb_moderate(y, s2, "standard_eb", "fitted")
  ad <- heb:::.heb_moderate(y, s2, "adaptive", "fitted")
  expect_gt(ad$weight[6], std$weight[6])
  expect_gt(ad$estimate[6], std$estimate[6])
})

# Two populations; donors D1, D2 measured in P2 everywhere, but P1 of batch s2_r1 only in D3, not in the anchor.
pop_gate_sim <- function() {
  s <- heb_sim()
  pop <- ifelse(seq_len(nrow(s$x)) %% 2 == 0, "P1", "P2")
  donor <- s$meta$donor[match(s$x$sample_id, s$meta$sample_id)]
  batch <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)]
  drop <- pop == "P1" &
    ((batch == "s2_r1" & donor != "D3") | (batch == "s1_r1" & donor == "D3"))
  list(s = s, keep = !drop, pop = pop)
}

test_that("a batch-population with no reference shared with the anchor gets no effect", {
  g <- pop_gate_sim()
  x <- g$s$x[g$keep, ]
  m <- fit_sim(g$s, x = x, population = g$pop[g$keep], unbridged = "abstain")
  expect_false(any(m$effects$batch == "s2_r1" & m$effects$population == "P1"))
  expect_true(any(m$effects$batch == "s2_r1" & m$effects$population == "P2"))
})

test_that("the log_sd centre ignores aliquots below min_cells", {
  s <- heb_sim()
  extra <- data.frame(sample_id = "D1_s2_r1_x", m1 = c(100, 100, 100))
  meta <- rbind(
    s$meta,
    data.frame(
      sample_id = "D1_s2_r1_x",
      batch = "s2_r1",
      donor = "D1",
      aliquot = "D1_s2_r1_x",
      reference_id = "D1",
      mu = 0
    )
  )
  base <- fit_sim(s, scale_model = "log_sd")$effects
  more <- fit_sim(
    s,
    x = rbind(s$x, extra),
    meta = meta,
    scale_model = "log_sd"
  )$effects
  b <- base$batch == "s2_r1"
  expect_equal(more$scale_centre[more$batch == "s2_r1"], base$scale_centre[b])
})

test_that("heb_fit needs min_cells of at least 2", {
  expect_error(fit_sim(heb_sim(), min_cells = 1), "min_cells")
})

test_that("new-batch paired rows of a moderated model carry identity moderation", {
  s <- heb_sim()
  s$x$m2 <- s$x$m1 + stats::rnorm(nrow(s$x), 0, 0.2)
  old <- s$meta$batch != "s3_r2"
  fit_s <- s
  fit_s$hierarchy <- s$hierarchy[s$hierarchy$batch != "s3_r2", ]
  m <- fit_sim(
    fit_s,
    x = s$x[s$x$sample_id %in% s$meta$sample_id[old], ],
    meta = s$meta[old, ],
    markers = c("m1", "m2"),
    anchor_effect = "fixed",
    residual_model = "known_plus_extra",
    moderation = "adaptive"
  )
  res <- heb_integrate(
    m,
    s$x,
    s$meta,
    rep("T", nrow(s$x)),
    transform = m$transform,
    unsupported = "error"
  )
  new <- res$model$effects[res$model$effects$batch == "s3_r2", ]
  expect_equal(new$estimate_raw, new$estimate)
  expect_equal(new$moderation_weight, c(1, 1))
  expect_equal(new$moderation_status, c("new_batch_paired", "new_batch_paired"))
})
