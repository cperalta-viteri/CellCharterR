torch_kmeans_fixture <- function() {
  # Deliberately well-separated, mildly asymmetric groups make the optimum
  # unambiguous while still exercising more than one Lloyd iteration.
  local_seed <- 2718L
  withr::with_seed(local_seed, {
    rbind(
      matrix(stats::rnorm(240L, sd = 0.30), ncol = 3L) +
        matrix(rep(c(-5, 0, 2), each = 80L), ncol = 3L),
      matrix(stats::rnorm(270L, sd = 0.35), ncol = 3L) +
        matrix(rep(c(1, 5, -2), each = 90L), ncol = 3L),
      matrix(stats::rnorm(300L, sd = 0.40), ncol = 3L) +
        matrix(rep(c(5, -4, 1), each = 100L), ncol = 3L)
    )
  })
}

torch_kmeans_objective <- function(x, centers) {
  sum(apply(.rc_kmeans_distances(x, centers), 1L, min))
}

test_that("Torch k-means initialization repeats and preserves the R RNG", {
  skip_if_not_installed("torch")
  x <- torch_kmeans_fixture()
  reader <- .rc_reader(x)

  set.seed(99L)
  expected_next <- stats::runif(5L)
  set.seed(99L)
  first <- .rc_torch_stream_kmeans(
    reader, 3L, 47L, 123L, torch_gmm_device("cpu"), max_iter = 30L
  )
  observed_next <- stats::runif(5L)
  repeated <- .rc_torch_stream_kmeans(
    reader, 3L, 47L, 123L, torch_gmm_device("cpu"), max_iter = 30L
  )

  expect_identical(observed_next, expected_next)
  expect_identical(repeated$means, first$means)
  expect_identical(repeated$iterations, first$iterations)
  expect_identical(repeated$converged, first$converged)
  expect_match(first$method, "[Tt]orch")
})

test_that("Torch CPU initialization is chunk-invariant and matches legacy quality", {
  skip_if_not_installed("torch")
  x <- torch_kmeans_fixture()
  reader <- .rc_reader(x)
  dev <- torch_gmm_device("cpu")

  small <- .rc_torch_stream_kmeans(
    reader, 3L, 31L, 123L, dev, max_iter = 30L
  )
  odd <- .rc_torch_stream_kmeans(
    reader, 3L, 73L, 123L, dev, max_iter = 30L
  )
  one <- .rc_torch_stream_kmeans(
    reader, 3L, nrow(x), 123L, dev, max_iter = 30L
  )
  legacy <- .rc_stream_kmeans(reader, 3L, 31L, 123L, max_iter = 30L)

  objective <- function(result) torch_kmeans_objective(x, result$means)
  # Streaming reductions may differ at float32 precision across chunk sizes,
  # but must converge to the same partition and essentially the same loss.
  labels <- function(result) {
    max.col(-.rc_kmeans_distances(x, result$means), ties.method = "first")
  }
  expect_equal(objective(odd), objective(small), tolerance = 1e-5)
  expect_equal(objective(one), objective(small), tolerance = 1e-5)
  expect_equal(
    cellcharter_fowlkes_mallows(labels(small), labels(odd)), 1,
    tolerance = 1e-7
  )
  expect_equal(
    cellcharter_fowlkes_mallows(labels(small), labels(one)), 1,
    tolerance = 1e-7
  )
  # Acceleration must not buy speed by materially weakening initialization.
  expect_lte(objective(small), objective(legacy) * 1.02 + 1e-6)
})

test_that("Torch CUDA initialization repeats and agrees with Torch CPU", {
  skip_if_not_installed("torch")
  if (identical(Sys.getenv("CELLCHARTERR_CPU_ONLY"), "true")) {
    skip("CUDA validation is outside this CPU-only test scope")
  }
  if (!isTRUE(torch::cuda_is_available())) {
    skip("CUDA is not available to R Torch")
  }
  x <- torch_kmeans_fixture()
  reader <- .rc_reader(x)
  run <- function(device, batch_size) {
    .rc_torch_stream_kmeans(
      reader, 3L, batch_size, 123L, torch_gmm_device(device), max_iter = 30L
    )
  }
  cpu <- run("cpu", 47L)
  cuda <- run("cuda", 47L)
  repeated_cuda <- run("cuda", 47L)

  expect_identical(repeated_cuda$means, cuda$means)
  expect_equal(
    torch_kmeans_objective(x, cuda$means),
    torch_kmeans_objective(x, cpu$means),
    tolerance = 1e-4
  )
  cpu_labels <- max.col(-.rc_kmeans_distances(x, cpu$means), ties.method = "first")
  cuda_labels <- max.col(-.rc_kmeans_distances(x, cuda$means), ties.method = "first")
  expect_equal(
    cellcharter_fowlkes_mallows(cpu_labels, cuda_labels), 1,
    tolerance = 1e-7
  )
})

test_that("torch_gmm_fit integrates the Torch initializer and auto batch", {
  skip_if_not_installed("torch")
  x <- torch_kmeans_fixture()[seq_len(90L), , drop = FALSE]
  fit <- torch_gmm_fit(
    x, n_clusters = 2L, batch_size = "auto", max_iter = 2L,
    device = "cpu", seed = 17L, num_threads = 1L,
    initialization_method = "torch"
  )
  expect_s3_class(fit, "torch_gmm")
  expect_match(fit$initialization$method, "Torch")
  expect_true(is.finite(fit$nll))
  expect_identical(fit$runtime$batch_size, 90L)
})
