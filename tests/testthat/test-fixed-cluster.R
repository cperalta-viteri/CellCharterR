.fixed_spe <- function() {
  ids <- paste0("spot", 1:8)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(0, 1, 8, dimnames = list("g", ids))),
    reducedDims = S4Vectors::SimpleList(
      X_fixed = matrix(seq_len(16), nrow = 8, dimnames = list(ids, c("a", "b")))
    ),
    spatialCoords = structure(
      cbind(x = seq_len(8), y = rep(0, 8)),
      dimnames = list(ids, c("x", "y"))
    ),
    sample_id = rep("sample1", 8)
  )
  spe
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
    .fixed_spe(), representation = "X_fixed", n_clusters = 3L,
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
    spe = .fixed_spe(), representation = "X_fixed", n_clusters = 2L,
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

test_that("fixed-K accepts first-class subset representations and rejects invalid labels", {
  stale <- .fixed_spe()[, 1:6]
  expect_silent(
    subset_result <- cellcharter_cluster_fixed(
      stale, representation = "X_fixed", n_clusters = 2L,
      fit_fun = function(x, n_clusters, seed, ...) list(k = n_clusters),
      predict_fun = function(model, newdata, ...) {
        (seq_len(nrow(newdata)) - 1L) %% model$k
      }
    )
  )
  expect_length(subset_result$cellcharter_domain, ncol(stale))
  expect_error(
    cellcharter_cluster_fixed(
      .fixed_spe(), representation = "X_fixed", n_clusters = 2L,
      fit_fun = function(x, n_clusters, seed, ...) list(k = n_clusters),
      predict_fun = function(model, newdata, ...) rep(2L, nrow(newdata))
    ),
    "zero-based label"
  )
})

test_that("fixed-K checkpoint identity follows representation content only", {
  checkpoint <- tempfile("fixed-content-")
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  fit <- function(x, n_clusters, seed, ...) {
    calls$n <- calls$n + 1L
    list(k = n_clusters)
  }
  predict <- function(model, newdata, ...) {
    (seq_len(nrow(newdata)) - 1L) %% model$k
  }
  spe <- .fixed_spe()
  cellcharter_cluster_fixed(
    spe, "X_fixed", 2L, checkpoint_dir = checkpoint,
    fit_fun = fit, predict_fun = predict
  )

  # Unrelated historical provenance cannot veto an unchanged representation.
  metadata <- S4Vectors::metadata(spe)
  metadata$CellCharterR <- list(stages = list(
    "aggregation:X_fixed" = list(filepath = "/moved/elsewhere")
  ))
  S4Vectors::metadata(spe) <- metadata
  resumed <- cellcharter_cluster_fixed(
    spe, "X_fixed", 2L, checkpoint_dir = checkpoint,
    fit_fun = fit, predict_fun = predict
  )
  expect_true(S4Vectors::metadata(resumed)$CellCharterR$fixed_k$cellcharter_domain$resumed)
  expect_identical(calls$n, 1L)

  changed <- spe
  value <- SingleCellExperiment::reducedDim(changed, "X_fixed")
  value[1L, 1L] <- value[1L, 1L] + 1
  SingleCellExperiment::reducedDim(changed, "X_fixed") <- value
  cellcharter_cluster_fixed(
    changed, "X_fixed", 2L, checkpoint_dir = checkpoint,
    fit_fun = fit, predict_fun = predict
  )
  expect_identical(calls$n, 2L)
})

test_that("fixed-K validates integers and does not reuse another implementation", {
  expect_error(
    cellcharter_cluster_fixed(.fixed_spe(), "X_fixed", 2.5),
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
    .fixed_spe(), "X_fixed", 2L, checkpoint_dir = checkpoint,
    fit_fun = fit_a, predict_fun = predictor
  )
  result <- cellcharter_cluster_fixed(
    .fixed_spe(), "X_fixed", 2L, checkpoint_dir = checkpoint,
    fit_fun = fit_b, predict_fun = predictor
  )
  expect_identical(calls$n, 2L)
  expect_true(all(as.character(result$cellcharter_domain) == "1"))
})
