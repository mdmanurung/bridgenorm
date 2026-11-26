# Guards added for 1.0.0: tied quantile tables, blank ids, multi-donor shares, min_cells gate, raw batches;
# the model contract. guard_sim() and guard_fit() are in helper-guards.R.

test_that("a marker with one value in a node keeps the node shift instead of stopping", {
  d <- guard_sim(zero_b = TRUE)
  model <- guard_fit(d)
  res <- bn_apply(model, d$x, d$samples, rep("T", nrow(d$x)))
  expect_true(all(is.finite(res$x$M)))
  tied <- res$actions[res$actions$reason_code == "node_quantile_tied", ]
  expect_gt(nrow(tied), 0)
  expect_true(all(tied$marker == "M"))
})

test_that("a blank reference_id stops with the sample names", {
  d <- guard_sim(study_ref = "")
  expect_error(guard_fit(d), "reference_id.*blank.*s1_b1")
  d$samples$aliquot <- NULL
  expect_error(guard_fit(d), "reference_id.*blank")
  expect_error(bn_bridges(d$samples, anchor = "b1"), "blank")
})

test_that("a carriage return in an id stops", {
  d <- guard_sim()
  d$samples$batch[d$samples$batch == "b3"] <- "b\r3"
  expect_error(guard_fit(d), "carriage return")
})

test_that("unequal bridge donor shares between a batch and the anchor warn", {
  d <- guard_sim(
    ref_ids = c("D1", "D2"),
    n_ref = c(50, 550),
    donor_off = c(0, 0.5)
  )
  expect_warning(guard_fit(d), "cell shares differ.*b2")
  d <- guard_sim(ref_ids = c("D1", "D2"), donor_off = c(0, 0.5))
  expect_no_warning(guard_fit(d))
})

test_that("a batch whose aliquots all fall below min_cells gets no quantile tables", {
  # b3 wells hold 30 cells per state (60 per well), below min_cells 100; pooled nodes would pass.
  d <- guard_sim(n_b3 = 30)
  expect_warning(model <- guard_fit(d, min_cells = 100), "No quantile maps.*b3")
  expect_false(any(startsWith(names(model$quantiles), "b3\r")))
})

test_that("bn_apply() names batches it left raw", {
  d <- guard_sim()
  model <- guard_fit(d)
  extra <- d$x[d$x$sample_id == "s1_b3", ]
  extra$sample_id <- "s1_b4"
  samples <- rbind(
    d$samples,
    data.frame(
      sample_id = "s1_b4",
      batch = "b4",
      reference_id = NA,
      aliquot = NA
    )
  )
  expect_message(
    bn_apply(
      model,
      rbind(d$x, extra),
      samples,
      rep("T", nrow(d$x) + nrow(extra))
    ),
    "no correction applied to batch\\(es\\) b4"
  )
  expect_no_message(bn_apply(model, d$x, d$samples, rep("T", nrow(d$x))))
})

# 1.0.0 model contract ------------------------------------------------------------

test_that("every model shares the bn_model parent class and a fixed field list", {
  d <- guard_sim()
  q <- guard_fit(d)
  p <- bn_fit_paired(
    d$x,
    d$samples,
    rep("T", nrow(d$x)),
    markers = c("S", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1"
  )
  expect_s3_class(q, c("bn_node_quantile_model", "bn_model"), exact = TRUE)
  expect_s3_class(p, c("bn_paired_model", "bn_model"), exact = TRUE)
  expect_s3_class(q$pass1, c("bn_node_model", "bn_model"), exact = TRUE)
  expect_identical(
    names(q),
    c(
      "schema_version",
      "settings",
      "markers",
      "transform",
      "anchor",
      "pass1",
      "pass2",
      "quantiles",
      "software",
      "created",
      "parent_digest",
      "digest"
    )
  )
  expect_identical(
    names(q$pass1),
    c(
      "schema_version",
      "method",
      "settings",
      "markers",
      "transform",
      "anchor",
      "map",
      "global_paired",
      "local",
      "software",
      "created",
      "parent_digest",
      "digest"
    )
  )
  expect_identical(
    names(p),
    c(
      "schema_version",
      "method",
      "settings",
      "markers",
      "transform",
      "population_levels",
      "anchor",
      "min_cells",
      "known_batches",
      "bridges",
      "effects",
      "targets",
      "excluded",
      "aliquot_source",
      "software",
      "created",
      "parent_digest",
      "digest"
    )
  )
  expect_identical(q$schema_version, "bn_node_quantile_model/2")
  expect_identical(p$schema_version, "bn_paired_model/2")
})

test_that("two fits on the same data give the same digest", {
  d <- guard_sim()
  a <- guard_fit(d)
  Sys.sleep(1.1)
  b <- guard_fit(d)
  expect_false(identical(a$created, b$created))
  expect_identical(a$digest, b$digest)
})

test_that("the fit does not depend on the collation locale", {
  d <- guard_sim()
  # Mixed-case ids sort differently in C and in most other locales.
  d$samples$sample_id <- sub("^ref1", "Ref1", d$samples$sample_id)
  d$x$sample_id <- sub("^ref1", "Ref1", d$x$sample_id)
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old))
  Sys.setlocale("LC_COLLATE", "C")
  a <- guard_fit(d)
  ok <- Sys.setlocale("LC_COLLATE", "en_US.UTF-8")
  skip_if(identical(ok, ""), "en_US.UTF-8 locale not available")
  b <- guard_fit(d)
  expect_identical(a$digest, b$digest)
})

test_that("seed must be a whole number", {
  d <- guard_sim()
  expect_error(
    bn_fit(
      d$x,
      d$samples,
      rep("T", nrow(d$x)),
      markers = c("S", "M"),
      transform = list(name = "asinh", cofactor = 5),
      anchor = "b1",
      grid = c(2, 1),
      mapping_markers = "S",
      seed = 1.7
    ),
    "whole number"
  )
})

test_that("a model saved by 1.0.0 still loads, validates and applies", {
  model <- bn_load(test_path("fixtures", "model-1.0.0.rds"))
  d <- guard_sim()
  res <- bn_apply(model, d$x, d$samples, rep("T", nrow(d$x)))
  expect_equal(
    res$x$M,
    readRDS(test_path("fixtures", "model-1.0.0-M.rds")),
    tolerance = 1e-12
  )
})
