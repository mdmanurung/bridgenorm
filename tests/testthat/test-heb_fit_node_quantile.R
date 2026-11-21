# heb_fit_node_quantile(): two-pass node correction with per-node quantile maps (tusoskill c016).

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
      heb_settings = NULL,
      grid = c(2, 1),
      mapping_markers = "CD3",
      train_cells = 200,
      rlen = 50,
      seed = 473,
      assignment = "nearest",
      fallback = "global",
      quantile_probs = seq(0.025, 0.975, length.out = 21),
      tail = "shift"
    ),
    list(...)
  )
  do.call(heb_fit_node_quantile, args)
}

integrate_sim <- function(model, s) {
  heb_integrate(
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
    heb_settings = NULL,
    grid = c(2, 1),
    mapping_markers = "CD3",
    train_cells = 200,
    rlen = 50,
    seed = 473,
    assignment = "nearest",
    fallback = "global",
    quantile_probs = seq(0.025, 0.975, length.out = 21),
    tail = "shift"
  )
  model <- do.call(heb_fit_node_quantile, args)
  res <- heb_integrate(
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
  heb_save(model, path)
  expect_identical(heb_load(path), model)
  bad <- model
  bad$quantiles[[1]][1, 1] <- bad$quantiles[[1]][1, 1] + 1
  saveRDS(bad, path)
  expect_error(heb_load(path), "Digest mismatch")
})

test_that("every estimate-changing argument is required", {
  skip_if_not_installed("kohonen")
  s <- quantile_sim()
  expect_error(
    fit_quantile_sim(s, quantile_probs = NULL),
    "quantile_probs"
  )
  expect_error(fit_quantile_sim(s, tail = "linear"), "tail")
  expect_error(
    fit_quantile_sim(s, quantile_probs = c(0.5, 0.2)),
    "quantile_probs"
  )
})
