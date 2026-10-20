cpp_gmm_optimization_fixture <- function() {
  set.seed(812)
  group <- rep(0:3, each = 18L)
  noise <- matrix(stats::rnorm(length(group) * 4L, sd = 0.18), ncol = 4L)
  centers <- rbind(
    c(-3.0, 2.0, -1.0, 0.5),
    c(0.0, -2.5, 2.0, 1.5),
    c(3.0, 1.0, 4.0, -2.0),
    c(6.0, 4.0, 0.5, 3.0)
  )
  noise + centers[group + 1L, ]
}

expect_cpp_gmm_equivalent <- function(observed, expected, tolerance = 1e-12) {
  expect_identical(observed$labels, expected$labels)
  expect_equal(observed$nll, expected$nll, tolerance = tolerance)
  expect_equal(observed$weights, expected$weights, tolerance = tolerance)
  expect_equal(observed$means, expected$means, tolerance = tolerance)
  expect_equal(observed$covariances, expected$covariances,
               tolerance = tolerance)
}

test_that("largest-K-first C++ scheduling preserves caller result order", {
  x <- cpp_gmm_optimization_fixture()
  ks <- c(2L, 5L, 3L)
  seeds <- c(102L, 505L, 303L)

  models <- cpp_gmm_fit_many(
    x, n_clusters = ks, seeds = seeds, workers = 2L,
    threads_per_fit = 1L, max_iter = 6L, tolerance = 1e-7,
    covariance_regularization = 1e-5, progress = FALSE
  )

  expect_identical(names(models), paste0("k", ks))
  expect_identical(
    unname(vapply(models, `[[`, integer(1L), "n_clusters")), ks
  )
  expect_identical(
    unname(vapply(models, function(model) model$runtime$seed, integer(1L))),
    seeds
  )
  expect_identical(
    unname(vapply(models, `[[`, numeric(1L), "seed")), as.numeric(seeds)
  )
})

test_that("optimized C++ fits are equivalent across worker topologies", {
  x <- cpp_gmm_optimization_fixture()
  arguments <- list(
    x = x, n_clusters = c(3L, 4L), seeds = c(413L, 714L),
    max_iter = 10L, tolerance = 1e-8,
    covariance_regularization = 1e-5, progress = FALSE
  )

  serial <- do.call(
    cpp_gmm_fit_many,
    c(arguments, list(workers = 1L, threads_per_fit = 1L))
  )
  outer_parallel <- do.call(
    cpp_gmm_fit_many,
    c(arguments, list(workers = 2L, threads_per_fit = 1L))
  )
  inner_parallel <- do.call(
    cpp_gmm_fit_many,
    c(arguments, list(workers = 1L, threads_per_fit = 2L))
  )

  for (i in seq_along(serial)) {
    expect_cpp_gmm_equivalent(outer_parallel[[i]], serial[[i]])
    expect_cpp_gmm_equivalent(inner_parallel[[i]], serial[[i]])
  }
})

test_that("optimized C++ covariance updates are stable under translation", {
  x <- cpp_gmm_optimization_fixture()
  translation <- c(1e8, -2e8, 3e8, -4e8)
  translated_x <- sweep(x, 2L, translation, "+")
  arguments <- list(
    n_clusters = 4L, max_iter = 15L, tolerance = 1e-8,
    covariance_regularization = 1e-5, seed = 919L, num_threads = 2L
  )

  original <- do.call(cpp_gmm_fit, c(list(x = x), arguments))
  translated <- do.call(cpp_gmm_fit, c(list(x = translated_x), arguments))

  expect_identical(translated$labels, original$labels)
  expect_equal(
    sweep(translated$means, 2L, translation, "-"), original$means,
    tolerance = 1e-6
  )
  expect_equal(translated$covariances, original$covariances,
               tolerance = 1e-7)
  expect_equal(translated$nll, original$nll, tolerance = 1e-7)
})

test_that("optimized covariance avoids cancellation for distant tight clusters", {
  observations_per_cluster <- 24L
  offsets <- as.numeric(scale(
    seq_len(observations_per_cluster), center = TRUE, scale = TRUE
  ))
  expected_group <- rep(0:1, each = observations_per_cluster)
  x <- matrix(c(-1e8 + offsets, 1e8 + offsets), ncol = 1L)
  translation <- 3e8
  arguments <- list(
    n_clusters = 2L, max_iter = 10L, tolerance = 1e-8,
    covariance_regularization = 1e-6, seed = 731L, num_threads = 2L
  )

  original <- do.call(cpp_gmm_fit, c(list(x = x), arguments))
  translated <- do.call(cpp_gmm_fit, c(list(x = x + translation), arguments))
  original_variances <- vapply(
    original$covariances, function(covariance) covariance[1L, 1L], numeric(1L)
  )
  translated_variances <- vapply(
    translated$covariances,
    function(covariance) covariance[1L, 1L], numeric(1L)
  )

  # Compare the partition rather than component numbers, which are arbitrary.
  expect_identical(
    original$labels == original$labels[[1L]],
    expected_group == expected_group[[1L]]
  )
  expect_true(is.finite(original$nll))
  expect_true(all(original_variances > 0.8 & original_variances < 1.2))

  expect_identical(translated$labels, original$labels)
  expect_equal(translated$means - translation, original$means,
               tolerance = 1e-6)
  expect_equal(translated_variances, original_variances, tolerance = 1e-7)
  expect_equal(translated$nll, original$nll, tolerance = 1e-8)
})

test_that("optimized C++ fitting and prediction do not dispatch to Torch", {
  x <- cpp_gmm_optimization_fixture()
  torch_called <- function(...) stop("Torch path was dispatched")
  testthat::local_mocked_bindings(
    .cc_require_torch = torch_called,
    torch_gmm_device = torch_called,
    torch_gmm_fit = torch_called,
    torch_gmm_predict = torch_called,
    torch_gmm_predict_proba = torch_called,
    torch_gmm_log_likelihood = torch_called,
    .package = "CellCharterR"
  )

  model <- cpp_gmm_fit(
    x, n_clusters = 4L, max_iter = 8L, tolerance = 1e-7,
    covariance_regularization = 1e-5, seed = 122L, num_threads = 2L
  )

  expect_identical(model$runtime$device, "cpu")
  expect_identical(model$runtime$backend, "compiled-cpp-double")
  expect_identical(predict(model, x, batch_size = 13L), model$labels)
  expect_equal(
    cpp_gmm_log_likelihood(model, x, batch_size = 11L), -model$nll,
    tolerance = 1e-10
  )
})
