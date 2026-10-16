cpp_gmm_fixture <- function() {
  set.seed(104)
  rbind(
    cbind(stats::rnorm(36L, -3, 0.35), stats::rnorm(36L, 1, 0.25)),
    cbind(stats::rnorm(36L, 2, 0.30), stats::rnorm(36L, -2, 0.40)),
    cbind(stats::rnorm(36L, 5, 0.25), stats::rnorm(36L, 4, 0.30))
  )
}

test_that("compiled scoring matches the pinned PyCave K=11 fixture", {
  fixture_path <- system.file(
    "extdata", "python_gmm_k11_fixture.rds", package = "CellCharterR"
  )
  if (!nzchar(fixture_path)) {
    fixture_path <- test_path(
      "..", "..", "inst", "extdata", "python_gmm_k11_fixture.rds"
    )
  }
  fixture <- readRDS(fixture_path)
  predicted <- cpp_gmm_predict(fixture$model, fixture$x, batch_size = 73L)
  expect_identical(predicted, fixture$expected_labels)
  expect_identical(fixture$provenance$backend, "PyCave 3.2.1 / PyTorch 1.12.1")
})

test_that("compiled GMM fits are deterministic across outer worker counts", {
  x <- cpp_gmm_fixture()
  arguments <- list(
    x = x, n_clusters = c(2L, 3L, 4L), seeds = c(11L, 12L, 13L),
    max_iter = 40L, tolerance = 1e-8, covariance_regularization = 1e-6
  )

  serial <- do.call(cpp_gmm_fit_many, c(arguments, list(workers = 1L)))
  parallel_two <- do.call(cpp_gmm_fit_many, c(arguments, list(workers = 2L)))
  parallel_four <- do.call(cpp_gmm_fit_many, c(arguments, list(workers = 4L)))
  inner_four <- do.call(cpp_gmm_fit_many, c(
    arguments, list(workers = 1L, threads_per_fit = 4L)
  ))

  for (i in seq_along(serial)) {
    expect_identical(parallel_two[[i]]$labels, serial[[i]]$labels)
    expect_identical(parallel_four[[i]]$labels, serial[[i]]$labels)
    expect_equal(parallel_two[[i]]$nll, serial[[i]]$nll, tolerance = 1e-12)
    expect_equal(parallel_four[[i]]$nll, serial[[i]]$nll, tolerance = 1e-12)
    expect_equal(parallel_two[[i]]$means, serial[[i]]$means, tolerance = 1e-12)
    expect_equal(parallel_four[[i]]$covariances, serial[[i]]$covariances,
                 tolerance = 1e-12)
    expect_identical(inner_four[[i]]$labels, serial[[i]]$labels)
    expect_equal(inner_four[[i]]$nll, serial[[i]]$nll, tolerance = 1e-10)
    expect_equal(inner_four[[i]]$means, serial[[i]]$means, tolerance = 1e-10)
    expect_equal(inner_four[[i]]$covariances, serial[[i]]$covariances,
                 tolerance = 1e-10)
  }
})

test_that("compiled GMM stored labels, likelihood, and parameters are coherent", {
  x <- cpp_gmm_fixture()
  ridge <- 1e-5
  model <- cpp_gmm_fit(x, 3L, max_iter = 50L, tolerance = 1e-8,
                       covariance_regularization = ridge, seed = 27L)

  expect_identical(model$labels, cpp_gmm_predict(model, x, batch_size = 17L))
  expect_equal(model$nll, -cpp_gmm_log_likelihood(model, x, batch_size = 19L),
               tolerance = 1e-10)
  expect_true(all(is.finite(model$weights)))
  expect_equal(sum(model$weights), 1, tolerance = 1e-12)
  expect_true(all(model$weights > 0))

  for (covariance in model$covariances) {
    expect_true(all(is.finite(covariance)))
    expect_equal(covariance, t(covariance), tolerance = 1e-12)
    regularized <- covariance + diag(ridge, ncol(x))
    expect_silent(chol(regularized))
    expect_gt(min(eigen(regularized, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
})

test_that("compiled scoring validates non-finite and empty inputs", {
  x <- cpp_gmm_fixture()
  model <- cpp_gmm_fit(x, 3L, max_iter = 5L, seed = 2L)
  expect_error(cpp_gmm_predict_proba(model, matrix(NA_real_, 1L, ncol(x))),
               "finite numeric")
  expect_error(cpp_gmm_log_likelihood(model, matrix(numeric(), 0L, ncol(x))),
               "at least one observation")
  empty <- cpp_gmm_predict(model, matrix(numeric(), 0L, ncol(x)))
  expect_identical(empty, integer())
})

test_that("compiled GMM is stable under a large common translation", {
  x <- cpp_gmm_fixture()
  offset <- matrix(c(1e9, -1e9), nrow(x), ncol(x), byrow = TRUE)
  original <- cpp_gmm_fit(x, 3L, max_iter = 50L, tolerance = 1e-8, seed = 91L)
  translated <- cpp_gmm_fit(x + offset, 3L, max_iter = 50L,
                            tolerance = 1e-8, seed = 91L)

  expect_identical(translated$labels, original$labels)
  expect_equal(sweep(translated$means, 2L, c(1e9, -1e9), "-"),
               original$means, tolerance = 2e-6)
  expect_equal(translated$covariances, original$covariances, tolerance = 2e-7)
  expect_equal(translated$nll, original$nll, tolerance = 2e-7)
})

test_that("compiled K=20 GMM handles collapsed and exactly degenerate observations", {
  x <- rbind(matrix(c(0, 0), 12L, 2L, byrow = TRUE),
             matrix(c(5, 5), 12L, 2L, byrow = TRUE))
  ridge <- 1e-4
  model <- cpp_gmm_fit(x, 20L, max_iter = 30L, tolerance = 1e-8,
                       covariance_regularization = ridge, seed = 3L)

  expect_true(is.finite(model$nll))
  expect_true(all(is.finite(model$weights)))
  expect_equal(sum(model$weights), 1, tolerance = 1e-12)
  expect_length(model$labels, nrow(x))
  expect_true(all(model$labels %in% 0:19))
  expect_true(length(model$initialization_collapsed_components) > 0L)
  for (covariance in model$covariances) {
    expect_silent(chol(covariance + diag(ridge, ncol(x))))
  }
})

test_that("checkpoint-free C++ AutoK creates no directory and embeds usable models", {
  x <- cpp_gmm_fixture()
  checkpoint_dir <- tempfile("cpp-autok-none-")
  fit <- cellcharter_autok_fit(
    x, n_clusters = 2:3, checkpoint_dir = checkpoint_dir,
    max_runs = 2L, convergence_tol = 0, seed = 44L,
    backend = "cpp", workers = 2L, checkpoint_mode = "none",
    fit_args = list(max_iter = 30L, tolerance = 1e-6), progress = FALSE
  )

  expect_false(dir.exists(checkpoint_dir))
  expect_identical(fit$checkpoint_mode, "none")
  expect_true(is.list(fit$embedded_best_artifacts))
  for (k in fit$fit_k) {
    model <- cellcharter_autok_model(fit, k)
    expect_s3_class(model, "torch_gmm")
    expect_identical(model$labels, cpp_gmm_predict(model, x, batch_size = 23L))
  }
})

test_that("round-checkpoint C++ AutoK resume is computationally equivalent", {
  x <- cpp_gmm_fixture()
  checkpoint_dir <- tempfile("cpp-autok-round-")
  arguments <- list(
    x = x, n_clusters = 2:3, checkpoint_dir = checkpoint_dir,
    max_runs = 2L, convergence_tol = 0, seed = 57L,
    backend = "cpp", workers = 2L, checkpoint_mode = "round",
    fit_args = list(max_iter = 30L, tolerance = 1e-6), progress = FALSE
  )
  first <- do.call(cellcharter_autok_fit, arguments)
  resumed <- do.call(cellcharter_autok_fit, arguments)

  expect_equal(resumed$stability, first$stability, tolerance = 0)
  expect_identical(resumed$best_k, first$best_k)
  expect_equal(resumed$best_index$nll, first$best_index$nll, tolerance = 0)
  expect_true(all(resumed$fit_index$resumed))
  for (k in first$fit_k) {
    expect_identical(cellcharter_autok_model(resumed, k)$labels,
                     cellcharter_autok_model(first, k)$labels)
  }
})
