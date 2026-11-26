# Guards added for 1.0.0: tied quantile tables, blank ids, multi-donor shares, min_cells gate, raw batches.

# Two cell states (A, B) on a state marker S and a corrected marker M, three batches, two bridge wells per batch.
guard_sim <- function(
  ref_ids = c("donor", "donor"),
  n_ref = c(300, 300),
  donor_off = c(0, 0),
  study_ref = NA,
  zero_b = FALSE,
  n_b3 = 300
) {
  set.seed(1)
  centre <- c(A = 2, B = 1)
  s_pos <- c(A = 1, B = 4)
  shift <- list(
    b1 = c(A = 0, B = 0),
    b2 = c(A = 0.4, B = -0.3),
    b3 = c(A = -0.2, B = 0.3)
  )
  scale <- c(b1 = 1, b2 = 0.6, b3 = 1.4)
  samples <- data.frame(
    sample_id = c(outer(
      c("ref1", "ref2", "s1", "s2"),
      c("b1", "b2", "b3"),
      paste,
      sep = "_"
    )),
    batch = rep(c("b1", "b2", "b3"), each = 4),
    reference_id = rep(c(ref_ids, study_ref, study_ref), 3),
    aliquot = rep(c("ref1", "ref2", NA, NA), 3)
  )
  cells <- do.call(
    rbind,
    lapply(seq_len(nrow(samples)), function(i) {
      b <- samples$batch[i]
      j <- (i - 1) %% 4 + 1
      n <- if (j <= 2 && b == "b2") {
        n_ref[j]
      } else if (j <= 2 && b == "b3") {
        n_b3
      } else {
        300
      }
      off <- if (j <= 2) donor_off[j] else mean(donor_off)
      do.call(
        rbind,
        lapply(c("A", "B"), function(k) {
          truth <- stats::rnorm(n, centre[[k]] + off, 0.3)
          m <- centre[[k]] +
            scale[[b]] * (truth - centre[[k]]) +
            shift[[b]][[k]]
          if (zero_b && k == "B") {
            m <- rep(0, n)
          }
          data.frame(
            sample_id = samples$sample_id[i],
            S = stats::rnorm(n, s_pos[[k]], 0.2),
            M = m
          )
        })
      )
    })
  )
  list(x = cells, samples = samples)
}

guard_fit <- function(d, min_cells = 50) {
  bn_fit(
    d$x,
    d$samples,
    rep("T", nrow(d$x)),
    markers = c("S", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1",
    grid = c(2, 1),
    mapping_markers = "S",
    seed = 1,
    train_cells = 300,
    rlen = 50,
    min_cells = min_cells
  )
}

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
