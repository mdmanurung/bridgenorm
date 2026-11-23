# bn_fit(): two-pass node correction with per-node quantile maps (tusoskill c016).

# Toy data: anchor b1 and batch b2, two reference donors per batch plus one study sample in b2.
# Two clusters on the mapping marker CD3 (A near 1, B near 4). In b2, marker M is a known monotone
# function of the anchor-scale value z: M = centre + 0.6 * (z - centre) + shift, a scale change the
# node shift cannot remove. The anchor-scale value z of every cell is kept as the expected value.
quantile_sim <- function(seed = 21, n = 400) {
  set.seed(seed)
  centre_m <- c(A = 2, B = 1)
  centre_cd3 <- c(A = 1, B = 4)
  shift <- c(A = 0.4, B = -0.3)
  samples <- data.frame(
    sample_id = c("D1_b1", "D2_b1", "D1_b2", "D2_b2", "S_b2"),
    batch = c("b1", "b1", "b2", "b2", "b2"),
    reference_id = c("D1", "D2", "D1", "D2", NA)
  )
  cells <- do.call(
    rbind,
    lapply(seq_len(nrow(samples)), function(i) {
      do.call(
        rbind,
        lapply(c("A", "B"), function(k) {
          z <- stats::rnorm(n, centre_m[[k]], 0.3)
          m <- if (samples$batch[i] == "b2") {
            centre_m[[k]] + 0.6 * (z - centre_m[[k]]) + shift[[k]]
          } else {
            z
          }
          data.frame(
            sample_id = samples$sample_id[i],
            CD3 = stats::rnorm(n, centre_cd3[[k]], 0.2),
            M = m,
            z = z
          )
        })
      )
    })
  )
  list(x = cells, samples = samples)
}

fit_quantile_sim <- function(s, ...) {
  args <- utils::modifyList(
    list(
      x = s$x[, c("sample_id", "CD3", "M")],
      samples = s$samples,
      population = rep("T", nrow(s$x)),
      markers = c("CD3", "M"),
      transform = list(name = "asinh", cofactor = 5),
      anchor = "b1",
      min_cells = 20,
      unbridged = "stop",
      method = "paired",
      eb_settings = NULL,
      grid = c(2, 1),
      mapping_markers = "CD3",
      train_cells = 200,
      rlen = 50,
      seed = 473,
      assignment = "nearest",
      fallback = "global",
      quantile_probs = seq(0.025, 0.975, length.out = 21),
      tail = "shift",
      shrink_k = 0
    ),
    list(...)
  )
  do.call(bn_fit, args)
}

integrate_sim <- function(model, s) {
  bn_apply(
    model,
    s$x[, c("sample_id", "CD3", "M")],
    s$samples,
    rep("T", nrow(s$x)),
    transform = list(name = "asinh", cofactor = 5),
    unsupported = "abstain"
  )
}

test_that("node quantile maps recover a known monotone batch transform that node shifts cannot", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  model <- fit_quantile_sim(s)
  res <- integrate_sim(model, s)
  study <- s$x$sample_id == "S_b2"
  # Expected value: each study cell's anchor-scale value z, known from the simulation.
  err_quantile <- mean(abs(res$x$M[study] - s$x$z[study]))
  shift_only <- fit_quantile_sim(s)$pass1
  err_shift <- mean(abs(integrate_sim(shift_only, s)$x$M[study] - s$x$z[study]))
  expect_lt(err_quantile, 0.03)
  expect_lt(err_quantile, err_shift / 3)
})

test_that("a value beyond the outer batch quantile moves by the outer quantile difference", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  model <- fit_quantile_sim(s)
  # One study cell placed far above every b2 reference value in cluster A.
  s$x$M[which(s$x$sample_id == "S_b2")[1]] <- 4
  res <- integrate_sim(model, s)
  pass1 <- integrate_sim(model$pass1, s)
  i <- which(s$x$sample_id == "S_b2")[1]
  q <- model$quantiles
  qb <- q[[paste("b2", res$nodes[i], sep = "\r")]]
  qa <- q[[paste("b1", res$nodes[i], sep = "\r")]]
  expect_equal(
    res$x$M[i] - pass1$x$M[i],
    unname(qa[nrow(qa), "M"] - qb[nrow(qb), "M"]),
    tolerance = 1e-12
  )
})

test_that("a (batch, node) below min_cells pooled reference cells keeps the pass-2 shift and says so", {
  skip_if_not_installed("kohonen")
  # b3 references: 2 wells x 30 cells per node = 60 pooled; b2: 2 x 150 = 300. min_cells 70 sits between.
  s <- nodes_sim(small_b3 = TRUE)
  args <- list(
    x = s$x[, c("sample_id", "CD3", "M")],
    samples = s$samples,
    population = rep("T", nrow(s$x)),
    markers = c("CD3", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1",
    min_cells = 70,
    unbridged = "abstain",
    method = "paired",
    eb_settings = NULL,
    grid = c(2, 1),
    mapping_markers = "CD3",
    train_cells = 200,
    rlen = 50,
    seed = 473,
    assignment = "nearest",
    fallback = "global",
    quantile_probs = seq(0.025, 0.975, length.out = 21),
    tail = "shift",
    shrink_k = 0
  )
  model <- do.call(bn_fit, args)
  res <- bn_apply(
    model,
    s$x[, c("sample_id", "CD3", "M")],
    s$samples,
    rep("T", nrow(s$x)),
    transform = list(name = "asinh", cofactor = 5),
    unsupported = "abstain"
  )
  p2 <- res$actions[res$actions$pass == 2, ]
  expect_true(all(p2$reason_code[p2$batch == "b2"] == "node_quantile"))
  expect_false(any(p2$reason_code[p2$batch == "b3"] == "node_quantile"))
  expect_false(any(grepl("^b3\r", names(model$quantiles))))
})

test_that("save and load keep the digest, and a changed quantile table fails", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  model <- fit_quantile_sim(s)
  path <- tempfile(fileext = ".rds")
  bn_save(model, path)
  expect_identical(bn_load(path), model)
  bad <- model
  bad$quantiles[[1]][1, 1] <- bad$quantiles[[1]][1, 1] + 1
  saveRDS(bad, path)
  expect_error(bn_load(path), "Digest mismatch")
})

test_that("data-specific settings are required and invalid settings are refused", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  # modifyList() drops an element set to NULL, so these calls omit the argument.
  for (arg in c("anchor", "grid", "mapping_markers", "seed", "transform")) {
    expect_error(do.call(fit_quantile_sim, c(list(s), stats::setNames(list(NULL), arg))), arg)
  }
  expect_error(fit_quantile_sim(s, tail = "linear"), "tail")
  expect_error(fit_quantile_sim(s, shrink_k = -1), "shrink_k")
  expect_error(fit_quantile_sim(s, quantile_probs = c(0.5, 0.2)), "quantile_probs")
  expect_error(fit_quantile_sim(s, eb_settings = list()), "eb_settings")
})

test_that("the defaults are the benchmarked settings", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  explicit <- fit_quantile_sim(s, min_cells = 50, unbridged = "abstain", train_cells = 500, rlen = 100, shrink_k = 300)
  defaults <- bn_fit(
    x = s$x[, c("sample_id", "CD3", "M")],
    samples = s$samples,
    population = rep("T", nrow(s$x)),
    markers = c("CD3", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1",
    grid = c(2, 1),
    mapping_markers = "CD3",
    seed = 473
  )
  expect_identical(defaults$quantiles, explicit$quantiles)
  expect_identical(defaults$settings, explicit$settings)
  expect_identical(defaults$settings$quantile_probs, seq(0.025, 0.975, length.out = 21))
})

# Expected tables from first principles: pass-1-corrected reference cells, their pass-2 nodes,
# stats::quantile() per batch and per (batch, node), the shrinkage formula, then each column sorted
# (rearrangement, heb 0.7.1).
expect_shrunk_tables <- function(s, raw, shrunk, k) {
  probs <- seq(0.025, 0.975, length.out = 21)
  ref <- s$x[!is.na(s$samples$reference_id[match(s$x$sample_id, s$samples$sample_id)]), c("sample_id", "CD3", "M")]
  integ <- function(m, cells) {
    bn_apply(m, cells, s$samples, rep("T", nrow(cells)), transform = list(name = "asinh", cofactor = 5), unsupported = "abstain")
  }
  x1 <- integ(raw$pass1, ref)$x
  nodes <- integ(raw$pass2, x1)$nodes
  batch <- s$samples$batch[match(x1$sample_id, s$samples$sample_id)]
  q_of <- function(i) sapply(c("CD3", "M"), function(m) stats::quantile(x1[[m]][i], probs, names = FALSE, type = 7))
  glob_off <- q_of(batch == "b2") - q_of(batch == "b1")
  for (node in unique(nodes)) {
    g <- paste("b2", node, sep = "\r")
    ga <- paste("b1", node, sep = "\r")
    skip_if(is.null(raw$quantiles[[g]]) || is.null(raw$quantiles[[ga]]))
    n <- sum(batch == "b2" & nodes == node)
    w <- n / (n + k)
    expected <- q_of(batch == "b1" & nodes == node) +
      w * (q_of(batch == "b2" & nodes == node) - q_of(batch == "b1" & nodes == node)) +
      (1 - w) * glob_off
    expect_equal(unname(shrunk$quantiles[[g]]), unname(apply(expected, 2, sort)), tolerance = 1e-12)
    expect_identical(shrunk$quantiles[[ga]], raw$quantiles[[ga]])
  }
}

test_that("shrink_k pulls each node offset curve toward the batch-wide curve by n / (n + shrink_k), then sorts", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  raw <- fit_quantile_sim(s, shrink_k = 0)
  # shrink_k 3000 makes some shrunk columns of this fixture non-monotone before sorting (min step -0.017).
  for (k in c(300, 3000)) {
    shrunk <- fit_quantile_sim(s, shrink_k = k)
    expect_shrunk_tables(s, raw, shrunk, k)
    steps <- unlist(lapply(shrunk$quantiles, function(q) apply(q, 2, function(v) min(diff(v)))))
    expect_true(all(steps >= 0))
    expect_identical(shrunk$settings$shrink_k, k)
  }
})

test_that("shrink_k = 0 leaves every quantile table as fitted", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  model <- fit_quantile_sim(s, shrink_k = 0)
  probs <- seq(0.025, 0.975, length.out = 21)
  ref <- s$x[!is.na(s$samples$reference_id[match(s$x$sample_id, s$samples$sample_id)]), c("sample_id", "CD3", "M")]
  x1 <- bn_apply(model$pass1, ref, s$samples, rep("T", nrow(ref)), transform = list(name = "asinh", cofactor = 5), unsupported = "abstain")$x
  nodes <- bn_apply(model$pass2, x1, s$samples, rep("T", nrow(x1)), transform = list(name = "asinh", cofactor = 5), unsupported = "abstain")$nodes
  batch <- s$samples$batch[match(x1$sample_id, s$samples$sample_id)]
  g <- names(model$quantiles)[[1]]
  i <- paste(batch, nodes, sep = "\r") == g
  expect_equal(
    unname(model$quantiles[[g]][, "M"]),
    stats::quantile(x1$M[i], probs, names = FALSE, type = 7),
    tolerance = 1e-12
  )
})

test_that("only batches with exactly the anchor's reference donors get quantile maps (RV2)", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  # b3: unbridged (donor D3 only). b4: partly bridged (D1 shared, D3 not). Cells copied from b2.
  copy_sample <- function(from, to, batch, ref) {
    list(
      x = transform(s$x[s$x$sample_id == from, ], sample_id = to),
      samples = data.frame(sample_id = to, batch = batch, reference_id = ref)
    )
  }
  extra <- list(
    copy_sample("D2_b2", "D3_b3", "b3", "D3"),
    copy_sample("S_b2", "S_b3", "b3", NA),
    copy_sample("D1_b2", "D1_b4", "b4", "D1"),
    copy_sample("D2_b2", "D3_b4", "b4", "D3"),
    copy_sample("S_b2", "S_b4", "b4", NA)
  )
  s$x <- rbind(s$x, do.call(rbind, lapply(extra, `[[`, "x")))
  s$samples <- rbind(s$samples, do.call(rbind, lapply(extra, `[[`, "samples")))
  model <- fit_quantile_sim(s, unbridged = "abstain", shrink_k = 300)
  expect_false(any(grepl("^b[34]\r", names(model$quantiles))))
  expect_true(any(grepl("^b2\r", names(model$quantiles))))
  res <- integrate_sim(model, s)
  # Unbridged b3 stays raw, as the README says.
  s3 <- s$x$sample_id == "S_b3"
  expect_equal(res$x$M[s3], s$x$M[s3], tolerance = 1e-12)
  # Partly bridged b4 takes pass 1 then the pass-2 node shift, with no quantile map.
  r1 <- integrate_sim(model$pass1, s)
  r2 <- bn_apply(model$pass2, r1$x, s$samples, rep("T", nrow(s$x)), transform = list(name = "asinh", cofactor = 5), unsupported = "abstain")
  s4 <- s$x$sample_id == "S_b4"
  expect_equal(res$x$M[s4], r2$x$M[s4], tolerance = 1e-12)
})
