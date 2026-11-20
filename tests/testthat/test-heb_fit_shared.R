# heb_fit_shared(): lme4 oracle, borrowing, abstention and invariance.

test_that("heb_fit_shared needs hierarchy and variance_floor", {
  s <- hebs_sim()
  expect_error(
    heb_fit_shared(
      s$x,
      s$meta,
      population = s$pop,
      markers = "m1",
      transform = TF,
      anchor = "s1_r1",
      hierarchy = NULL,
      min_cells = 10,
      variance_floor = 1e-8
    ),
    "hierarchy"
  )
  expect_error(
    heb_fit_shared(
      s$x,
      s$meta,
      population = s$pop,
      markers = "m1",
      transform = TF,
      anchor = "s1_r1",
      hierarchy = s$hierarchy,
      min_cells = 10
    ),
    "variance_floor"
  )
})

test_that("heb_fit_shared effects match an lme4 REML fit with a batch-by-population term", {
  skip_if_not_installed("lme4")
  s <- hebs_sim()
  eff <- fit_shared_sim(s)$effects
  d <- stats::aggregate(
    m1 ~ sample_id + pop,
    data = cbind(s$x, pop = s$pop),
    FUN = mean
  )
  d <- merge(d, s$meta, by = "sample_id")
  d <- merge(d, s$hierarchy, by = "batch")
  d$dp <- paste(d$donor, d$pop)
  d$bp <- paste(d$batch, d$pop)
  fit <- lme4::lmer(
    m1 ~ 0 + dp + (1 | study) + (1 | batch) + (1 | bp),
    data = d,
    REML = TRUE,
    control = lme4::lmerControl(
      check.conv.grad = lme4::.makeCC("ignore", tol = 1e-3)
    )
  )
  re <- lme4::ranef(fit)
  beta <- function(b, k) {
    re$study[substr(b, 1, 2), 1] + re$batch[b, 1] + re$bp[paste(b, k), 1]
  }
  for (k in c("P", "Q")) {
    want <- vapply(
      s$hierarchy$batch,
      function(b) beta(b, k) - beta("s1_r1", k),
      numeric(1)
    )
    e <- eff[eff$population == k, ]
    expect_equal(
      unname(e$estimate[match(s$hierarchy$batch, e$batch)]),
      unname(want),
      tolerance = 2e-3
    )
    expect_equal(e$estimate[e$batch == "s1_r1"], 0)
  }
})

test_that("a population dropped by min_cells in a batch borrows the shared effect", {
  s <- hebs_sim(few = c("s2_r2", "Q"))
  eff <- fit_shared_sim(s)$effects
  q <- eff[eff$population == "Q" & eff$batch == "s2_r2", ]
  p <- eff[eff$population == "P" & eff$batch == "s2_r2", ]
  expect_equal(nrow(q), 1)
  expect_equal(q$support, "shared_borrowed")
  expect_equal(q$n_donors, 0L)
  expect_true(is.finite(q$estimate) && is.finite(q$se_donor))
  expect_gt(q$se_donor, p$se_donor)
  # The per-population model has no effect there.
  v0 <- fit_sim(s, population = s$pop)$effects
  expect_equal(nrow(v0[v0$population == "Q" & v0$batch == "s2_r2", ]), 0)
})

test_that("a population with no anchor summary abstains", {
  s <- hebs_sim(drop_anchor_pop = "Q")
  model <- fit_shared_sim(s)
  expect_false("Q" %in% model$effects$population)
  res <- heb_integrate(model, s$x, s$meta, s$pop, TF, unsupported = "abstain")
  expect_true(all(
    res$actions$action[res$actions$population == "Q"] == "abstain"
  ))
  expect_true(all(res$actions$action[res$actions$population == "P"] == "apply"))
})

test_that("heb_fit_shared stops on a batch that shares no reference with the anchor", {
  s <- hebs_sim()
  s$meta$reference_id[s$meta$batch == "s3_r2"] <- NA
  s$meta$aliquot[s$meta$batch == "s3_r2"] <- NA
  expect_error(fit_shared_sim(s), "share no reference.*s3_r2")
})
