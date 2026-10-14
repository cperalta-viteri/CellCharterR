.fixed_sce <- function() {
  ids <- paste0("cell", 1:8)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, 1, 8, dimnames = list("g", ids))),
    reducedDims = S4Vectors::SimpleList(
      X_fixed = matrix(seq_len(16), nrow = 8, dimnames = list(ids, c("a", "b")))
    )
  )
  CellCharterR:::.cc_sce_record_stage(
    sce, "aggregation:X_fixed", list(output = "X_fixed")
  )
}

test_that("fixed-K clustering fits exactly one requested K and attaches labels", {
  seen <- new.env(parent = emptyenv())
  seen$calls <- integer()
  fit <- function(x, n_clusters, seed, ...) {
    seen$calls <- c(seen$calls, n_clusters)
    list(k = n_clusters, seed = seed)
  }
  predict <- function(model, newdata, ...) (seq_len(nrow(newdata)) - 1L) %% model$k

  result <- cellcharter_cluster_fixed(
    .fixed_sce(), representation = "X_fixed", n_clusters = 3L,
    domain_col = "domain", seed = 9L, fit_fun = fit, predict_fun = predict
  )

  expect_identical(seen$calls, 3L)
  expect_identical(levels(result$domain), c("0", "1", "2"))
  info <- S4Vectors::metadata(result)$CellCharterR$fixed_k$domain
  expect_identical(info$n_clusters, 3L)
  expect_false(info$resumed)
  expect_true(is.list(info$model))
})

test_that("fixed-K compatible checkpoint resumes without refitting", {
  checkpoint <- tempfile("fixed-checkpoint-")
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  fit <- function(x, n_clusters, seed, ...) {
    calls$n <- calls$n + 1L
    list(k = n_clusters, seed = seed)
  }
  predict <- function(model, newdata, ...) (seq_len(nrow(newdata)) - 1L) %% model$k
  args <- list(
    sce = .fixed_sce(), representation = "X_fixed", n_clusters = 2L,
    checkpoint_dir = checkpoint, fit_fun = fit, predict_fun = predict
  )
  first <- do.call(cellcharter_cluster_fixed, args)
  second <- do.call(cellcharter_cluster_fixed, args)

  expect_identical(calls$n, 1L)
  expect_identical(as.character(first$cellcharter_domain),
                   as.character(second$cellcharter_domain))
  expect_true(S4Vectors::metadata(second)$CellCharterR$fixed_k$cellcharter_domain$resumed)
  expect_true(file.exists(file.path(checkpoint, "fixed_k2.rds")))
})

test_that("fixed-K rejects stale stages and invalid labels", {
  stale <- .fixed_sce()[, 1:6]
  expect_error(
    cellcharter_cluster_fixed(stale, representation = "X_fixed", n_clusters = 2L),
    "stale"
  )
  expect_error(
    cellcharter_cluster_fixed(
      .fixed_sce(), representation = "X_fixed", n_clusters = 2L,
      fit_fun = function(x, n_clusters, seed, ...) list(k = n_clusters),
      predict_fun = function(model, newdata, ...) rep(2L, nrow(newdata))
    ),
    "zero-based label"
  )
})

test_that("fixed-K validates integers and does not reuse another implementation", {
  expect_error(
    cellcharter_cluster_fixed(.fixed_sce(), "X_fixed", 2.5),
    "n_clusters"
  )
  checkpoint <- tempfile("fixed-implementation-")
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  predictor <- function(model, newdata, ...) rep.int(model$label, nrow(newdata))
  fit_a <- function(x, n_clusters, seed, ...) {
    calls$n <- calls$n + 1L
    list(label = 0L)
  }
  fit_b <- function(x, n_clusters, seed, ...) {
    calls$n <- calls$n + 1L
    list(label = 1L)
  }
  cellcharter_cluster_fixed(
    .fixed_sce(), "X_fixed", 2L, checkpoint_dir = checkpoint,
    fit_fun = fit_a, predict_fun = predictor
  )
  result <- cellcharter_cluster_fixed(
    .fixed_sce(), "X_fixed", 2L, checkpoint_dir = checkpoint,
    fit_fun = fit_b, predict_fun = predictor
  )
  expect_identical(calls$n, 2L)
  expect_true(all(as.character(result$cellcharter_domain) == "1"))
})
