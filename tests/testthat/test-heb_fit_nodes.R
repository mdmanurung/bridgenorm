# heb_fit_nodes(): node-level shifts recover planted per-cluster batch effects; fallback, order, seed, save/load.

fit_nodes_sim <- function(s, ...) {
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
      heb_settings = NULL,
      grid = c(2, 1),
      mapping_markers = "CD3",
      train_cells = 200,
      rlen = 50,
      seed = 473,
      assignment = "nearest",
      fallback = "global"
    ),
    list(...)
  )
  do.call(heb_fit_nodes, args)
}

integrate_nodes_sim <- function(model, s) {
  heb_integrate(
    model,
    s$x[, c("sample_id", "CD3", "M")],
    s$samples,
    population = rep("T", nrow(s$x)),
    transform = list(name = "asinh", cofactor = 5),
    unsupported = "abstain"
  )
}

# Mean corrected M per study sample and cluster; the planted base values are A = 2, B = 1.
study_means <- function(res, s) {
  study <- s$x$sample_id %in% s$samples$sample_id[is.na(s$samples$reference_id)]
  stats::aggregate(
    res$x$M[study],
    by = list(sample_id = s$x$sample_id[study], cluster = s$x$cluster[study]),
    FUN = mean
  )
}

test_that("paired node model recovers planted per-cluster shifts", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  res <- integrate_nodes_sim(fit_nodes_sim(s), s)
  m <- study_means(res, s)
  expected <- ifelse(m$cluster == "A", 2, 1)
  expect_equal(m$x, expected, tolerance = 0.03, ignore_attr = TRUE)
})

test_that("the global CD3+ model cannot recover cluster-specific shifts (control)", {
  s <- nodes_sim()
  global <- heb_fit_paired(
    s$x[, c("sample_id", "CD3", "M")],
    s$samples,
    population = rep("T", nrow(s$x)),
    markers = c("CD3", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1",
    min_cells = 20,
    unbridged = "stop"
  )
  res <- heb_integrate(
    global,
    s$x[, c("sample_id", "CD3", "M")],
    s$samples,
    population = rep("T", nrow(s$x)),
    transform = list(name = "asinh", cofactor = 5),
    unsupported = "abstain"
  )
  m <- study_means(res, s)
  b2 <- m[m$sample_id == "S_b2", ]
  # One shift of about +0.1 (mean of +0.5 and -0.3) leaves cluster A 0.4 high and B 0.4 low.
  expect_gt(max(abs(b2$x - ifelse(b2$cluster == "A", 2, 1))), 0.3)
})

test_that("heb node model recovers planted per-cluster shifts", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  settings <- list(
    hierarchy = data.frame(batch = c("b1", "b2", "b3"), study = "s1"),
    variance_floor = 1e-8,
    pool = "none",
    anchor_effect = "fixed",
    variance_prior = "boundary_avoiding",
    residual_model = "known_plus_extra",
    scale_model = "none",
    moderation = "none",
    moderation_centre = "fitted",
    moderation_variance = "pev",
    batch_shrinkage = "blup",
    batch_shrinkage_centre = "study",
    batch_shrinkage_group = "study",
    batch_shrinkage_tau2 = "profile"
  )
  res <- integrate_nodes_sim(
    fit_nodes_sim(s, method = "heb", heb_settings = settings),
    s
  )
  m <- study_means(res, s)
  expect_equal(
    m$x,
    ifelse(m$cluster == "A", 2, 1),
    tolerance = 0.05,
    ignore_attr = TRUE
  )
})

test_that("a node without an effect falls back to the global shift or abstains", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim(small_b3 = TRUE)
  b3 <- s$x$sample_id == "S_b3"

  # Global fallback: b3's shift is +0.2 in both clusters, so the global shift recovers A = 2 and B = 1.
  res <- integrate_nodes_sim(
    fit_nodes_sim(s, fallback = "global", min_cells = 40),
    s
  )
  m <- study_means(res, s)
  b3_means <- m[m$sample_id == "S_b3", ]
  expect_equal(
    b3_means$x,
    ifelse(b3_means$cluster == "A", 2, 1),
    tolerance = 0.03,
    ignore_attr = TRUE
  )
  expect_true(all(
    res$actions$reason_code[res$actions$batch == "b3"] == "node_fallback_global"
  ))

  # Abstain: raw values kept and flagged.
  res <- integrate_nodes_sim(
    fit_nodes_sim(s, fallback = "abstain", min_cells = 40),
    s
  )
  expect_identical(res$x$M[b3], s$x$M[b3])
  expect_true(any(res$actions$action == "abstain"))
})

test_that("results follow cells, not row order", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  model <- fit_nodes_sim(s)
  res <- integrate_nodes_sim(model, s)
  perm <- sample(nrow(s$x))
  shuffled <- s
  shuffled$x <- s$x[perm, ]
  res_shuffled <- integrate_nodes_sim(model, shuffled)
  expect_equal(res_shuffled$x$M, res$x$M[perm])
})

test_that("one seed gives the same map; the seed is required", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  expect_identical(fit_nodes_sim(s)$map, fit_nodes_sim(s)$map)
  args <- list(
    x = s$x[, c("sample_id", "CD3", "M")],
    samples = s$samples,
    population = rep("T", nrow(s$x)),
    markers = c("CD3", "M"),
    transform = list(name = "asinh", cofactor = 5),
    anchor = "b1",
    min_cells = 20,
    unbridged = "stop",
    method = "paired",
    heb_settings = NULL,
    grid = c(2, 1),
    mapping_markers = "CD3",
    train_cells = 200,
    rlen = 50,
    assignment = "nearest",
    fallback = "global"
  )
  expect_error(do.call(heb_fit_nodes, args), "seed")
})

test_that("a node model survives heb_save and heb_load", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  model <- fit_nodes_sim(s)
  path <- tempfile(fileext = ".rds")
  heb_save(model, path)
  loaded <- heb_load(path)
  expect_equal(
    integrate_nodes_sim(loaded, s)$x,
    integrate_nodes_sim(model, s)$x
  )
  tampered <- model
  tampered$map$nodes[1, 1] <- tampered$map$nodes[1, 1] + 1
  expect_error(heb_save(tampered, path), "Digest mismatch")
})

# Review fixes for 0.5.1 (plans/2026-11-21_heb-0.5.1-fixes.md) -----------------

test_that("mapping markers outside the corrected markers are refused", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  expect_error(
    fit_nodes_sim(s, markers = "M", mapping_markers = "CD3"),
    "mapping_markers"
  )
})

test_that("a node row reads 'apply' only when every cell in it got the global effect", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim(small_b3 = TRUE)
  # Populations P and Q alternate by row, so every node holds both.
  pop <- rep_len(c("P", "Q"), nrow(s$x))
  batch <- s$samples$batch[match(s$x$sample_id, s$samples$sample_id)]
  is_ref <- !is.na(s$samples$reference_id[match(
    s$x$sample_id,
    s$samples$sample_id
  )])
  # Keep 5 Q cells per b3 reference sample: no global (b3, Q) effect, while (b3, P) keeps 30.
  rank_in_group <- stats::ave(
    seq_along(pop),
    s$x$sample_id,
    pop,
    FUN = seq_along
  )
  keep <- !(batch == "b3" & is_ref & pop == "Q" & rank_in_group > 5)
  s$x <- s$x[keep, ]
  pop <- pop[keep]
  batch <- batch[keep]
  # min_cells 25: each b3 node has 20 cells per sample (15 P + 5 Q), so no node effect; global P has 30.
  model <- fit_nodes_sim(s, population = pop, min_cells = 25)
  res <- heb_integrate(
    model,
    s$x[, c("sample_id", "CD3", "M")],
    s$samples,
    population = pop,
    transform = list(name = "asinh", cofactor = 5),
    unsupported = "abstain"
  )
  b3_rows <- res$actions[
    res$actions$batch == "b3" & res$actions$marker == "M",
  ]
  expect_true(all(b3_rows$action == "abstain"))
  expect_true(all(b3_rows$reason_code == "node_fallback_partial"))
  # Every b3 Q cell keeps its raw M; no b3 P cell does.
  raw <- res$x$M == s$x$M
  expect_equal(
    sum(raw[batch == "b3" & pop == "Q"]),
    sum(batch == "b3" & pop == "Q")
  )
  expect_equal(sum(raw[batch == "b3" & pop == "P"]), 0)
  expect_error(
    heb_integrate(
      model,
      s$x[, c("sample_id", "CD3", "M")],
      s$samples,
      population = pop,
      transform = list(name = "asinh", cofactor = 5),
      unsupported = "error"
    ),
    "cannot be corrected"
  )
})

test_that("a node model keeps the effects of a batch it estimated at integrate time", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  old <- s$samples$batch != "b3"
  s_old <- list(
    x = s$x[s$x$sample_id %in% s$samples$sample_id[old], ],
    samples = s$samples[old, ]
  )
  model <- fit_nodes_sim(s_old)
  first <- integrate_nodes_sim(model, s)
  expect_identical(first$model$parent_digest, model$digest)
  # The b3 study sample alone, with the returned model: corrected to the planted base A = 2, B = 1.
  s_study <- list(x = s$x[s$x$sample_id == "S_b3", ], samples = s$samples)
  second <- integrate_nodes_sim(first$model, s_study)
  expect_true(all(second$actions$action == "apply"))
  m <- study_means(second, s_study)
  expect_equal(
    m$x,
    ifelse(m$cluster == "A", 2, 1),
    tolerance = 0.03,
    ignore_attr = TRUE
  )
})

test_that("heb_fit_nodes() leaves the caller's random number stream alone", {
  skip_if_not_installed("kohonen")
  s <- nodes_sim()
  set.seed(1)
  want <- stats::runif(1)
  set.seed(1)
  invisible(fit_nodes_sim(s))
  expect_identical(stats::runif(1), want)
})
