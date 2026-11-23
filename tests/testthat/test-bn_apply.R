# bn_apply(), the frozen model and save/load.

test_that("integrating again leaves earlier corrected values unchanged", {
  s <- bn_sim()
  model <- fit_sim(s)
  old <- bn_apply(
    model,
    s$x,
    s$meta,
    rep("T", nrow(s$x)),
    TF,
    unsupported = "abstain"
  )
  again <- bn_apply(
    old$model,
    s$x,
    s$meta,
    rep("T", nrow(s$x)),
    TF,
    unsupported = "abstain"
  )
  expect_identical(old$x, again$x)
})

test_that("a new bridged batch gets a paired effect and old effects stay frozen", {
  s <- bn_sim()
  in_new <- s$meta$batch[match(s$x$sample_id, s$meta$sample_id)] == "s3_r2"
  s_old <- subset_sim(s, s$meta$batch != "s3_r2")
  model <- fit_sim(s_old, hierarchy = s$hierarchy)
  res <- bn_apply(
    model,
    s$x,
    s$meta,
    rep("T", nrow(s$x)),
    TF,
    unsupported = "abstain"
  )
  new_eff <- res$model$effects[res$model$effects$batch == "s3_r2", ]
  expect_equal(new_eff$source, "integrate")
  expect_equal(new_eff$support, "paired_reference")
  expect_identical(res$model$parent_digest, model$digest)
  expect_equal(
    res$model$effects[res$model$effects$batch != "s3_r2", ],
    model$effects
  )
  # The new effect is the mean paired difference to the frozen anchor summaries.
  dm <- function(b) {
    ids <- s$meta$sample_id[s$meta$batch == b]
    v <- tapply(
      s$x$m1[s$x$sample_id %in% ids],
      s$meta$reference_id[match(
        s$x$sample_id[s$x$sample_id %in% ids],
        s$meta$sample_id
      )],
      mean
    )
    v
  }
  expect_equal(
    new_eff$estimate,
    mean(dm("s3_r2") - dm("s1_r1")),
    tolerance = 1e-10
  )
  expect_false(isTRUE(all.equal(res$x$m1[in_new], s$x$m1[in_new])))
})

test_that("integrate refuses a different transform and reports abstentions under 'error'", {
  s <- bn_sim()
  model <- fit_sim(s)
  expect_error(
    bn_apply(
      model,
      s$x,
      s$meta,
      rep("T", nrow(s$x)),
      list(name = "asinh", cofactor = 6000),
      unsupported = "abstain"
    ),
    "Transform mismatch"
  )
  expect_error(
    bn_apply(
      model,
      s$x,
      s$meta,
      rep("Z", nrow(s$x)),
      TF,
      unsupported = "error"
    ),
    "cannot be corrected"
  )
})

test_that("save and load keep the digest; a changed model fails validation", {
  s <- bn_sim()
  model <- fit_sim(s)
  path <- tempfile(fileext = ".rds")
  bn_save(model, path)
  expect_identical(bn_load(path)$digest, model$digest)
  bad <- model
  bad$effects$estimate[2] <- bad$effects$estimate[2] + 1
  expect_error(bn_save(bad, path), "Digest mismatch")
})

test_that("bn_bridges lists batches without a reference shared with the anchor", {
  samples <- data.frame(
    sample_id = paste0("s", 1:5),
    batch = c("b1", "b1", "b2", "b3", "b3"),
    reference_id = c("R1", "R2", "R2", "R3", NA)
  )
  br <- bn_bridges(samples, anchor = "b1")
  expect_equal(br$batch, c("b2", "b3"))
  expect_equal(br$bridged, c(TRUE, FALSE))
  expect_equal(br$shared_ids, c("R2", ""))
})
