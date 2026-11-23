# bn_fit(): required settings, lme4 oracle, limiting cases, scale term, unbridged batches.

test_that("bn_fit names every missing required argument", {
  s <- bn_sim()
  expect_error(
    bn_fit(s$x, s$meta, population = rep("T", nrow(s$x)), markers = "m1"),
    "transform.*anchor.*hierarchy.*min_cells.*unbridged.*variance_floor.*pool"
  )
  expect_error(fit_sim(s, hierarchy = NULL), "hierarchy")
  expect_error(fit_sim(s, variance_floor = 0), "variance_floor")
  expect_error(fit_sim(s, pool = "some"), "pool")
  expect_error(fit_sim(s, variance_prior = "flat"), "variance_prior")
  expect_error(
    fit_sim(s, scale_model = "sd"),
    "scale_model = 'none' or 'log_sd'"
  )
  expect_error(fit_sim(s, unbridged = "maybe"), "unbridged")
})

test_that("pooling over populations changes nothing with one population and marker", {
  s <- bn_sim()
  expect_equal(
    fit_sim(s, pool = "population")$effects$estimate,
    fit_sim(s)$effects$estimate,
    tolerance = 1e-6
  )
})

test_that("bn_fit effects match an lme4 REML fit on the same summaries", {
  skip_if_not_installed("lme4")
  s <- bn_sim()
  eff <- fit_sim(s)$effects
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
  re <- lme4::ranef(fit)
  beta <- function(b) re$study[substr(b, 1, 2), 1] + re$batch[b, 1]
  want <- vapply(
    s$hierarchy$batch,
    function(b) beta(b) - beta("s1_r1"),
    numeric(1)
  )
  got <- eff$estimate[match(s$hierarchy$batch, eff$batch)]
  expect_equal(unname(got), unname(want), tolerance = 1e-3)
  expect_equal(eff$estimate[eff$batch == "s1_r1"], 0)
})

test_that("the boundary prior keeps one-donor SEs off zero", {
  s1 <- subset_sim(bn_sim(), bn_sim()$meta$donor == "D1")
  eff <- fit_sim(s1, variance_prior = "boundary_avoiding")$effects
  expect_true(all(eff$se_donor[eff$batch != "s1_r1"] > 1e-3))
})

test_that("known sampling variance sets the one-donor, one-aliquot SE", {
  s <- bn_sim()
  s <- subset_sim(s, s$meta$donor == "D1" & endsWith(s$meta$sample_id, "_1"))
  eff <- fit_sim(s, residual_model = "known_plus_extra")$effects
  se <- eff$se_donor[eff$batch != "s1_r1"]
  expect_true(all(se > 0.5 * sqrt(2) * 0.4 / sqrt(40) / 5))
  expect_true(all(is.finite(eff$estimate)))
})

test_that("bn_fit fits a design with one reference donor", {
  s1 <- subset_sim(bn_sim(), bn_sim()$meta$donor == "D1")
  eff <- fit_sim(s1)$effects
  expect_equal(nrow(eff), 6)
  expect_true(all(is.finite(eff$estimate)))
})

test_that("Satterthwaite df reduces to the residual df without shrinkage", {
  # Batch variances far above the noise: df = 2 / Var(log sigma^2) = df_resid = 6 - 3 - 1 = 2.
  d <- .bn_design(
    data.frame(
      value = c(1, 2.1, 3, 2, 3, 4.2),
      reference_id = rep(c("D1", "D2", "D3"), 2),
      batch = rep(c("a", "b"), each = 3),
      aliquot = 1:6
    ),
    data.frame(batch = c("a", "b"), study = c("s", "s")),
    "a",
    "random",
    "estimated"
  )
  est <- .bn_blup(d, c(1, 1, 1e-8), "a", cov_log = diag(c(0, 0, 2 / 2)))
  expect_equal(est$df[est$batch == "b"], 2, tolerance = 1e-3)
  expect_true(is.na(est$df[est$batch == "a"]))
  eff <- fit_sim(bn_sim())$effects
  expect_true(all(eff$df[eff$batch != "s1_r1"] > 0))
  expect_true(all(is.na(eff$df[eff$batch == "s1_r1"])))
})

test_that("anchor_effect = 'fixed' pins the anchor at zero and reduces to paired shift", {
  s <- bn_sim()
  fixed <- fit_sim(s, anchor_effect = "fixed")$effects
  expect_equal(fixed$estimate[fixed$batch == "s1_r1"], 0)
  expect_equal(fixed$se_donor[fixed$batch == "s1_r1"], 0)
  d <- .bn_design(
    data.frame(
      value = c(1, 2, 3, 2, 3, 4),
      reference_id = rep(c("D1", "D2", "D3"), 2),
      batch = rep(c("a", "b"), each = 3),
      aliquot = 1:6
    ),
    data.frame(batch = c("a", "b"), study = c("s", "s")),
    "a",
    "fixed",
    "estimated"
  )
  expect_equal(
    .bn_blup(d, c(1e3, 1e3, 1e-3), "a")$estimate,
    c(0, 1),
    tolerance = 1e-3
  )
})

test_that("unbridged = 'stop' refuses and 'abstain' leaves the unbridged batch raw", {
  s <- bn_sim()
  # Batch s3_r2 keeps its cells, but its samples are not references.
  meta2 <- s$meta
  no_ref <- meta2$batch == "s3_r2"
  meta2$reference_id[no_ref] <- NA
  meta2$aliquot[no_ref] <- NA
  expect_error(fit_sim(s, meta = meta2), "share no reference.*s3_r2")
  model <- fit_sim(
    s,
    meta = meta2,
    unbridged = "abstain",
    variance_prior = "boundary_avoiding",
    residual_model = "known_plus_extra"
  )
  expect_setequal(model$effects$batch, setdiff(s$hierarchy$batch, "s3_r2"))
  out <- bn_integrate(
    model,
    s$x,
    meta2,
    rep("T", nrow(s$x)),
    TF,
    unsupported = "abstain"
  )
  in_b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)] == "s3_r2"
  expect_equal(out$x$m1[in_b], s$x$m1[in_b])
  ab <- out$actions[out$actions$batch == "s3_r2", ]
  expect_true(all(ab$action == "abstain"))
  expect_true(all(ab$reason_code == "unbridged_at_fit"))
  expect_false(isTRUE(all.equal(out$x$m1[!in_b], s$x$m1[!in_b])))
})

test_that("a batch whose references are not in the anchor batch gets no effect under abstain", {
  s <- bn_sim()
  meta2 <- s$meta
  # s3_r2 keeps references, but under ids the anchor batch never measured.
  meta2$reference_id[meta2$batch == "s3_r2"] <- paste0(
    "X",
    meta2$reference_id[meta2$batch == "s3_r2"]
  )
  model <- fit_sim(s, meta = meta2, unbridged = "abstain")
  expect_false("s3_r2" %in% model$effects$batch)
  expect_false(model$bridges$bridged[model$bridges$batch == "s3_r2"])
})

# Scale term on log SD (plan F12) -----------------------------------------------------------------
scale_sim <- function(spread) {
  s <- bn_sim()
  mu <- s$meta$mu[match(s$x$sample_id, s$meta$sample_id)]
  b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)]
  f <- ifelse(b == "s2_r1", spread, 1)
  s$x$m1 <- mu + f * (s$x$m1 - mu)
  s
}

test_that("log_sd scale is about 1 without a spread effect and recovers a 1.5x spread", {
  none <- fit_sim(
    scale_sim(1),
    residual_model = "known_plus_extra",
    scale_model = "log_sd"
  )$effects
  expect_true(all(abs(log(none$scale)) < 0.1))
  sp <- fit_sim(
    scale_sim(1.5),
    residual_model = "known_plus_extra",
    scale_model = "log_sd"
  )$effects
  e <- sp[sp$batch == "s2_r1", ]
  expect_lt(abs(-log(e$scale) - log(1.5)), 3 * e$scale_se)
  expect_equal(sp$scale[sp$batch == "s1_r1"], 1)
})

test_that("log_sd correction evens the spread, and scale_model = 'none' is the plain shift", {
  s <- scale_sim(1.5)
  fit <- fit_sim(s, residual_model = "known_plus_extra", scale_model = "log_sd")
  out <- bn_integrate(
    fit,
    s$x,
    s$meta,
    rep("T", nrow(s$x)),
    TF,
    unsupported = "abstain"
  )
  b <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)]
  within_sd <- function(bb) {
    v <- out$x$m1[b == bb]
    stats::sd(v - ave(v, s$x$sample_id[b == bb]))
  }
  expect_lt(abs(within_sd("s2_r1") / within_sd("s1_r1") - 1), 0.15)
  plain <- fit_sim(s, residual_model = "known_plus_extra", scale_model = "none")
  expect_false("scale" %in% names(plain$effects))
  expect_equal(plain$effects$estimate, fit$effects$estimate)
})

# Review fixes for 0.5.1 (plans/2026-11-21_heb-0.5.1-fixes.md) -----------------

test_that("n_aliquots counts aliquots per donor, not shared labels", {
  s <- bn_sim()
  # Label aliquots "1" and "2" within each donor and batch: 3 donors x 2 = 6 aliquots per batch.
  s$meta$aliquot <- sub(".*_", "", s$meta$sample_id)
  eff <- fit_sim(s)$effects
  expect_true(all(eff$n_aliquots == 6))
  expect_true(all(eff$n_aliquots_anchor == 6))
})

test_that("a missing study label stops with a clear message", {
  s <- bn_sim()
  s$hierarchy$study[2] <- NA
  expect_error(fit_sim(s), "study")
})
