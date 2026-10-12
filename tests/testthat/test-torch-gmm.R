.gmm_parity_fixture <- function() {
  set.seed(271828)
  z <- matrix(stats::rnorm(360L * 3L), ncol = 3L)
  transforms <- list(
    matrix(c(0.80, 0, 0, 0.35, 0.55, 0, -0.20, 0.15, 0.45), 3L, byrow = TRUE),
    matrix(c(0.65, 0, 0, -0.25, 0.75, 0, 0.15, -0.30, 0.55), 3L, byrow = TRUE),
    matrix(c(0.70, 0, 0, 0.10, 0.60, 0, 0.30, 0.20, 0.80), 3L, byrow = TRUE)
  )
  locations <- list(c(-5, -1, 1), c(0, 5, -2), c(5, -3, 2.5))
  do.call(rbind, lapply(seq_along(transforms), function(component) {
    rows <- seq.int((component - 1L) * 120L + 1L, component * 120L)
    sweep(z[rows, , drop = FALSE] %*% transforms[[component]], 2L,
          locations[[component]], FUN = "+")
  }))
}

.gmm_permutations <- function(values) {
  if (length(values) == 1L) return(matrix(values, nrow = 1L))
  do.call(rbind, lapply(values, function(first) {
    cbind(first, .gmm_permutations(values[values != first]))
  }))
}

.gmm_component_match <- function(reference, candidate) {
  stopifnot(reference$n_clusters == candidate$n_clusters)
  candidates <- .gmm_permutations(seq_len(reference$n_clusters))
  losses <- apply(candidates, 1L, function(order) {
    sum((reference$means - candidate$means[order, , drop = FALSE])^2)
  })
  as.integer(candidates[which.min(losses), ])
}

.gmm_fowlkes_mallows <- function(first, second) {
  contingency <- table(first, second)
  true_positive <- sum(contingency * (contingency - 1) / 2)
  first_positive <- sum(rowSums(contingency) * (rowSums(contingency) - 1) / 2)
  second_positive <- sum(colSums(contingency) * (colSums(contingency) - 1) / 2)
  true_positive / sqrt(first_positive * second_positive)
}

.gmm_parity_metrics <- function(reference, candidate, x,
                                reference_device = "cpu",
                                candidate_device = "cpu") {
  component_order <- .gmm_component_match(reference, candidate)
  reference_labels <- torch_gmm_predict(
    reference, x, batch_size = 29L, device = reference_device
  )
  candidate_labels <- torch_gmm_predict(
    candidate, x, batch_size = 29L, device = candidate_device
  )
  candidate_to_reference <- integer(reference$n_clusters)
  candidate_to_reference[component_order] <- seq_len(reference$n_clusters)
  matched_candidate_labels <- candidate_to_reference[candidate_labels + 1L] - 1L
  list(
    component_order = component_order,
    fmi = .gmm_fowlkes_mallows(reference_labels, candidate_labels),
    matched_label_agreement = mean(reference_labels == matched_candidate_labels),
    max_mean_error = max(abs(
      reference$means - candidate$means[component_order, , drop = FALSE]
    )),
    max_covariance_error = max(vapply(
      seq_along(component_order),
      function(component) max(abs(
        reference$covariances[[component]] -
          candidate$covariances[[component_order[[component]]]]
      )),
      numeric(1L)
    )),
    max_weight_error = max(abs(
      reference$weights - candidate$weights[component_order]
    )),
    nll_error = abs(reference$nll - candidate$nll)
  )
}

.expect_gmm_float32_parity <- function(reference, candidate, x,
                                       reference_device = "cpu",
                                       candidate_device = "cpu") {
  metrics <- .gmm_parity_metrics(
    reference, candidate, x, reference_device, candidate_device
  )
  # Float32 EM reductions change order across CUDA kernels and chunk layouts.
  # These absolute bounds are hundreds of float32 epsilons at the fixture's
  # order-one scale, yet <0.04% of its smallest within-component SD (about
  # 0.45), so they admit reduction noise without masking a scientific change.
  expect_gte(metrics$fmi, 0.9999)
  expect_gte(metrics$matched_label_agreement, 0.9999)
  expect_lt(metrics$max_mean_error, 2e-4)
  expect_lt(metrics$max_covariance_error, 2e-4)
  expect_lt(metrics$max_weight_error, 2e-5)
  expect_lt(metrics$nll_error, 2e-5)
  invisible(metrics)
}

test_that("Torch GMM separates a simple full-covariance mixture", {
  skip_if_not_installed("torch")
  set.seed(7)
  x <- rbind(MASS::mvrnorm(300, c(-2, 0), matrix(c(1, .6, .6, 1), 2)),
             MASS::mvrnorm(300, c(2, 0), matrix(c(1, -.4, -.4, 1), 2)))
  gmm <- torch_gmm_fit(x, 2, batch_size = 64, max_iter = 60, device = "cpu", seed = 7)
  expect_true(gmm$converged || gmm$iterations == 60)
  labels <- torch_gmm_predict(gmm, x, batch_size = 64, device = "cpu")
  expect_length(labels, nrow(x))
  expect_true(all(labels %in% 0:1))
  expect_equal(ncol(torch_gmm_predict_proba(gmm, x, 64, "cpu")), 2)
  mean_order <- order(gmm$means[, 1])
  expect_gt(gmm$covariances[[mean_order[[1L]]]][1, 2], 0.15)
  expect_lt(gmm$covariances[[mean_order[[2L]]]][1, 2], -0.05)
})

test_that("centered covariance is translation-stable across odd batches", {
  skip_if_not_installed("torch")
  i <- seq.int(-15, 15)
  x <- cbind(i, ((seq_along(i) - 1L) * 7L) %% 31L - 15L)
  shifted <- x + 2^20
  expected <- crossprod(sweep(shifted, 2, colMeans(shifted))) / nrow(shifted)

  fit <- torch_gmm_fit(
    shifted, 1, batch_size = 7, max_iter = 1, device = "cpu", seed = 23,
    num_threads = 1
  )
  unshifted_fit <- torch_gmm_fit(
    x, 1, batch_size = 31, max_iter = 1, device = "cpu", seed = 23,
    num_threads = 1
  )

  expect_equal(unname(drop(fit$means)), unname(colMeans(shifted)), tolerance = 1e-6)
  expect_lt(max(abs(fit$covariances[[1]] - expected)), 0.02)
  expect_lt(
    max(abs(fit$covariances[[1]] - unshifted_fit$covariances[[1]])), 0.02
  )
})

test_that("fractional effective memberships use one normalization", {
  skip_if_not_installed("torch")
  device <- torch::torch_device("cpu")
  nk_cpu <- c(0, 0.25, 3.75)
  nk <- .rc_tensor(nk_cpu, device)
  sx <- .rc_tensor(rbind(c(0, 0), c(0.5, -0.25), c(15, 7.5)), device)

  update <- .rc_mstep_means_weights(
    nk, sx, nk_cpu, collapsed = 1L, collapse_threshold = 1e-6,
    device = device
  )
  means <- torch::as_array(update$means$to(device = "cpu"))
  weights <- as.numeric(torch::as_array(update$weights$to(device = "cpu")))

  expect_equal(means[2:3, , drop = FALSE], rbind(c(2, -1), c(4, 2)))
  expect_equal(weights, c(1e-6, 0.25, 3.75) / (4 + 1e-6), tolerance = 1e-7)
})

test_that("hard initialization uses nearest-center counts, means, and centered covariance", {
  skip_if_not_installed("torch")
  x <- rbind(c(0, 0), c(2, 0), c(10, 0), c(14, 0))
  device <- torch::torch_device("cpu")
  eye <- .rc_eye(2L, device)
  initialized <- .rc_hard_initialize(
    .rc_reader(x), batch_size = 3L,
    initialized_means = rbind(c(0, 0), c(10, 0)),
    device = device, eye = eye
  )

  expect_equal(initialized$counts, c(2, 2))
  expect_equal(
    torch::as_array(initialized$means$to(device = "cpu")),
    rbind(c(1, 0), c(12, 0))
  )
  expect_equal(
    as.numeric(torch::as_array(initialized$weights$to(device = "cpu"))),
    c(0.5, 0.5)
  )
  expected <- list(diag(c(1, 0)), diag(c(4, 0)))
  observed <- lapply(initialized$covariances, function(value) {
    torch::as_array(value$to(device = "cpu"))
  })
  expect_equal(observed, expected)

  # PyCave's CovarianceAggregator adds one global ridge after dividing the
  # centered sums by the hard-assignment counts.
  chol <- .rc_cholesky(initialized$covariances, 0.25, eye)
  regularized <- lapply(chol, function(value) {
    lower <- torch::as_array(value$to(device = "cpu"))
    lower %*% t(lower)
  })
  expect_equal(regularized, lapply(expected, `+`, diag(0.25, 2)), tolerance = 1e-6)
})

test_that("one streamed EM update matches independent full-batch equations", {
  skip_if_not_installed("torch")
  x <- rbind(
    c(-4, 0), c(-3, 1), c(-2, -1), c(-3, 0.25),
    c(3, 0), c(4, 1), c(5, -1), c(4, -0.25)
  )
  ridge <- 1e-4
  fit <- torch_gmm_fit(
    x, 2L, batch_size = 3L, max_iter = 1L, tolerance = 1e-8,
    covariance_regularization = ridge, device = "cpu", seed = 41L,
    num_threads = 1L
  )

  log_scores <- function(data, means, covariances, weights) {
    vapply(seq_len(nrow(means)), function(component) {
      covariance <- covariances[[component]] + diag(ridge, ncol(data))
      centered <- sweep(data, 2L, means[component, ], FUN = "-")
      precision_distance <- rowSums((centered %*% solve(covariance)) * centered)
      log(weights[[component]]) - 0.5 * (
        ncol(data) * log(2 * pi) + determinant(covariance, logarithm = TRUE)$modulus +
          precision_distance
      )
    }, numeric(nrow(data)))
  }
  centroids <- fit$initialization$centers
  initial_assignment <- max.col(
    -.rc_kmeans_distances(x, centroids), ties.method = "first"
  )
  initial_weights <- tabulate(initial_assignment, nbins = 2L) / nrow(x)
  initial_means <- do.call(rbind, lapply(1:2, function(component) {
    colMeans(x[initial_assignment == component, , drop = FALSE])
  }))
  initial_covariances <- lapply(1:2, function(component) {
    member <- x[initial_assignment == component, , drop = FALSE]
    centered <- sweep(member, 2L, initial_means[component, ], FUN = "-")
    crossprod(centered) / nrow(member)
  })
  initial_scores <- log_scores(
    x, initial_means, initial_covariances, initial_weights
  )
  row_max <- apply(initial_scores, 1L, max)
  responsibilities <- exp(initial_scores - row_max)
  responsibilities <- responsibilities / rowSums(responsibilities)
  expected_nk <- colSums(responsibilities)
  expected_weights <- expected_nk / nrow(x)
  expected_means <- crossprod(responsibilities, x) / expected_nk
  expected_covariances <- lapply(1:2, function(component) {
    centered <- sweep(x, 2L, expected_means[component, ], FUN = "-")
    crossprod(centered * sqrt(responsibilities[, component])) /
      expected_nk[[component]]
  })
  expected_scores <- log_scores(
    x, expected_means, expected_covariances, expected_weights
  )
  expected_nll <- -mean(
    apply(expected_scores, 1L, function(value) {
      maximum <- max(value)
      maximum + log(sum(exp(value - maximum)))
    })
  )

  expect_equal(fit$weights, expected_weights, tolerance = 2e-5)
  expect_equal(fit$means, expected_means, tolerance = 2e-5)
  expect_equal(fit$covariances, expected_covariances, tolerance = 2e-5)
  expect_equal(fit$nll, expected_nll, tolerance = 2e-5)
  expect_identical(
    torch_gmm_predict(fit, x, batch_size = 3L, device = "cpu"),
    max.col(expected_scores, ties.method = "first") - 1L
  )

  full_chunk <- torch_gmm_fit(
    x, 2L, batch_size = nrow(x), max_iter = 1L, tolerance = 1e-8,
    covariance_regularization = ridge, device = "cpu", seed = 41L,
    num_threads = 1L
  )
  expect_equal(full_chunk$weights, fit$weights, tolerance = 2e-5)
  expect_equal(full_chunk$means, fit$means, tolerance = 2e-5)
  expect_identical(
    torch_gmm_predict(full_chunk, x, device = "cpu"),
    torch_gmm_predict(fit, x, device = "cpu")
  )
})

test_that("early stopping monitors pre-update NLL with min_delta and patience two", {
  skip_if_not_installed("torch")
  x <- rbind(c(-2, 0), c(-1, 1), c(0, -1), c(1, 1), c(2, -1))
  fit <- torch_gmm_fit(
    x, 1L, batch_size = 2L, max_iter = 10L, tolerance = 0.1,
    device = "cpu", seed = 3L, num_threads = 1L
  )

  expect_true(fit$converged)
  expect_identical(fit$iterations, 3L)
  expect_length(fit$nll_history, 3L)
  expect_equal(diff(fit$nll_history), c(0, 0), tolerance = 1e-6)
  expect_identical(fit$early_stopping$monitor, "mean_nll_before_m_step")
  expect_equal(fit$early_stopping$min_delta, 0.1)
  expect_identical(fit$early_stopping$patience, 2L)
  expect_identical(fit$early_stopping$wait_count, 2L)
})

test_that("Cholesky failure restarts the complete fit with one effective ridge", {
  skip_if_not_installed("torch")
  # The two identical feature columns make the float32 covariance rank one.
  # Ridges below float32 resolution are swallowed, forcing global restarts.
  x <- cbind(-3:3, -3:3)
  fit <- suppressMessages(torch_gmm_fit(
    x, 1L, batch_size = 3L, max_iter = 3L,
    covariance_regularization = 1e-12, device = "cpu", seed = 1L,
    num_threads = 1L
  ))

  expect_equal(fit$requested_covariance_regularization, 1e-12)
  expect_identical(fit$regularization_max_attempts, 8L)
  expect_lte(length(fit$regularization_attempts), fit$regularization_max_attempts)
  expect_gt(fit$covariance_regularization, 1e-12)
  expect_equal(
    fit$covariance_regularization,
    tail(fit$regularization_attempts, 1L)
  )
  expect_equal(
    fit$regularization_attempts[-length(fit$regularization_attempts)] * 10,
    fit$regularization_attempts[-1L],
    tolerance = 1e-15
  )
  expect_length(
    fit$regularization_restart_errors,
    length(fit$regularization_attempts) - 1L
  )
  expect_true(all(grepl("Cholesky failed for component 1", fit$regularization_restart_errors)))

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_identical(
    restored$covariance_regularization, fit$covariance_regularization
  )
  expect_identical(
    torch_gmm_predict(restored, x, batch_size = 2L, device = "cpu"),
    torch_gmm_predict(fit, x, batch_size = 4L, device = "cpu")
  )
  expect_equal(
    torch_gmm_predict_proba(restored, x, batch_size = 2L, device = "cpu"),
    torch_gmm_predict_proba(fit, x, batch_size = 4L, device = "cpu"),
    tolerance = 1e-7
  )
  expect_equal(
    fit$nll, -torch_gmm_log_likelihood(restored, x, 2L, "cpu"),
    tolerance = 1e-6
  )
})

test_that("Torch GMM is reproducible and serializable on CPU", {
  skip_if_not_installed("torch")
  set.seed(11)
  x <- rbind(
    MASS::mvrnorm(80, c(-3, 1), matrix(c(1, 0.7, 0.7, 1.5), 2)),
    MASS::mvrnorm(80, c(3, -1), matrix(c(1.3, -0.5, -0.5, 1), 2))
  )
  fit1 <- torch_gmm_fit(x, 2, batch_size = 37, max_iter = 30,
                        device = "cpu", seed = 11, num_threads = 1)
  fit2 <- torch_gmm_fit(x, 2, batch_size = 37, max_iter = 30,
                        device = "cpu", seed = 11, num_threads = 1)
  labels1 <- torch_gmm_predict(fit1, x, batch_size = 29, device = "cpu")
  expect_equal(labels1, torch_gmm_predict(fit2, x, 29, "cpu"))
  expect_equal(fit1$nll, fit2$nll, tolerance = 1e-6)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit1, path)
  expect_equal(
    torch_gmm_predict(readRDS(path), x, 31, "cpu"), labels1
  )
})

test_that("full-covariance fits are deterministic on CPU and chunk invariant", {
  skip_if_not_installed("torch")
  x <- .gmm_parity_fixture()
  fit <- function(batch_size) suppressMessages(torch_gmm_fit(
    x, 3L, batch_size = batch_size, max_iter = 35L, tolerance = 1e-5,
    covariance_regularization = 1e-5, device = "cpu", seed = 31415L,
    num_threads = 1L
  ))
  first <- fit(37L)
  repeated <- fit(37L)
  odd_boundary <- fit(41L)
  one_batch <- fit(nrow(x))

  expect_identical(first$runtime$device, "cpu")
  expect_identical(repeated$means, first$means)
  expect_identical(repeated$covariances, first$covariances)
  expect_identical(repeated$weights, first$weights)
  expect_identical(repeated$nll, first$nll)
  expect_identical(
    torch_gmm_predict(repeated, x, batch_size = 31L, device = "cpu"),
    torch_gmm_predict(first, x, batch_size = 43L, device = "cpu")
  )
  .expect_gmm_float32_parity(first, odd_boundary, x)
  .expect_gmm_float32_parity(first, one_batch, x)
})

test_that("full-covariance fits repeat on CUDA and agree with CPU", {
  skip_if_not_installed("torch")
  if (identical(Sys.getenv("CELLCHARTERR_CPU_ONLY"), "true")) {
    skip("CUDA validation is outside this CPU-only acceptance scope")
  }
  if (!isTRUE(torch::cuda_is_available())) {
    skip("CUDA is not available to R Torch; CPU reproducibility was exercised")
  }
  x <- .gmm_parity_fixture()
  fit <- function(device, batch_size) suppressMessages(torch_gmm_fit(
    x, 3L, batch_size = batch_size, max_iter = 35L, tolerance = 1e-5,
    covariance_regularization = 1e-5, device = device, seed = 31415L,
    num_threads = 1L
  ))
  cpu <- fit("cpu", 37L)
  cuda <- fit("cuda", 37L)
  repeated_cuda <- fit("cuda", 37L)
  odd_boundary_cuda <- fit("cuda", 41L)

  expect_identical(cuda$runtime$device, "cuda")
  expect_true(cuda$runtime$cuda_available)
  expect_true(is.character(cuda$runtime$gpu_name))
  expect_identical(repeated_cuda$means, cuda$means)
  expect_identical(repeated_cuda$covariances, cuda$covariances)
  expect_identical(repeated_cuda$weights, cuda$weights)
  expect_identical(repeated_cuda$nll, cuda$nll)
  expect_identical(
    torch_gmm_predict(repeated_cuda, x, batch_size = 31L, device = "cuda"),
    torch_gmm_predict(cuda, x, batch_size = 43L, device = "cuda")
  )
  .expect_gmm_float32_parity(
    cpu, cuda, x, reference_device = "cpu", candidate_device = "cuda"
  )
  .expect_gmm_float32_parity(
    cuda, odd_boundary_cuda, x,
    reference_device = "cuda", candidate_device = "cuda"
  )
})

test_that("Torch GMM auto-device smoke test returns finite likelihood", {
  skip_if_not_installed("torch")
  set.seed(19)
  x <- rbind(matrix(rnorm(60, -2), ncol = 3), matrix(rnorm(60, 2), ncol = 3))
  requested_device <- if (identical(Sys.getenv("CELLCHARTERR_CPU_ONLY"), "true")) {
    "cpu"
  } else {
    "auto"
  }
  fit <- torch_gmm_fit(x, 2, batch_size = 16, max_iter = 8,
                       device = requested_device, seed = 19, num_threads = 1)
  labels <- torch_gmm_predict(fit, x, batch_size = 13, device = requested_device)
  expect_true(is.finite(fit$nll))
  expect_true(is.finite(torch_gmm_log_likelihood(fit, x, 17, requested_device)))
  expect_length(labels, nrow(x))
  expect_true(all(labels %in% 0:1))
  expect_identical(fit$runtime$cuda_available, torch::cuda_is_available())
  if (identical(requested_device, "cpu")) {
    expect_identical(fit$runtime$device, "cpu")
  } else if (torch::cuda_is_available()) {
    expect_identical(fit$runtime$device, "cuda")
  }
})

test_that("R prediction matches the pinned PyCave K=11 fixture", {
  skip_if_not_installed("torch")
  path <- system.file(
    "extdata", "python_gmm_k11_fixture.rds", package = "CellCharterR",
    mustWork = TRUE
  )
  fixture <- readRDS(path)
  predicted <- torch_gmm_predict(
    fixture$model, fixture$x, batch_size = 73L, device = "cpu"
  )
  expect_identical(predicted, fixture$expected_labels)
  expect_identical(
    fixture$provenance$cellcharter_commit,
    "2a8fdde8213e32e5a00f1f9571bffd3dbb825a79"
  )
  expect_identical(fixture$provenance$backend, "PyCave 3.2.1 / PyTorch 1.12.1")
  expect_equal(fixture$provenance$expected_agreement, 1)
})

test_that("full R fit matches an independent Python EM fixture from fixed centroids", {
  skip_if_not_installed("torch")
  path <- system.file(
    "extdata", "python_gmm_fit_fixture.json", package = "CellCharterR",
    mustWork = TRUE
  )
  fixture <- jsonlite::read_json(path, simplifyVector = TRUE)
  x <- as.matrix(fixture$x)
  initial_centers <- as.matrix(fixture$initial_centers)
  testthat::local_mocked_bindings(
    .rc_stream_kmeans = function(reader, k, batch_size, seed,
                                 max_iter = 100L, tolerance = 1e-4) {
      list(
        means = initial_centers, iterations = 0L, converged = TRUE,
        method = "controlled Python fixture centroids", seed = seed
      )
    },
    .package = "CellCharterR"
  )
  fit <- torch_gmm_fit(
    x, n_clusters = 2L, batch_size = nrow(x),
    max_iter = fixture$parameters$max_iter,
    tolerance = fixture$parameters$tolerance,
    covariance_regularization = fixture$parameters$covariance_regularization,
    device = "cpu", seed = 17L, num_threads = 1L
  )

  expect_equal(fit$weights, fixture$weights, tolerance = 2e-5)
  expect_equal(fit$means, fixture$means, tolerance = 2e-5)
  for (component in seq_along(fit$covariances)) {
    expect_equal(
      fit$covariances[[component]], fixture$covariances[component, , ],
      tolerance = 2e-5
    )
  }
  expect_equal(fit$nll, fixture$nll, tolerance = 2e-5)
  expect_identical(
    torch_gmm_predict(fit, x, batch_size = 3L, device = "cpu"),
    as.integer(fixture$labels)
  )
  expect_identical(fit$iterations, as.integer(fixture$iterations))
})

test_that("Torch GMM rejects invalid fitting controls before EM", {
  skip_if_not_installed("torch")
  x <- matrix(rnorm(30), ncol = 3)
  expect_error(torch_gmm_fit(x, 0), "n_clusters")
  expect_error(torch_gmm_fit(x, 2, batch_size = 0), "batch_size")
  expect_error(torch_gmm_fit(x, 2, max_iter = 0), "max_iter")
  expect_error(torch_gmm_fit(x, 2, tolerance = 0), "tolerance")
  expect_error(
    torch_gmm_fit(x, 2, covariance_regularization = 0),
    "covariance_regularization"
  )
  x[1, 1] <- NA_real_
  expect_error(torch_gmm_fit(x, 2), "non-finite")
})

test_that("Torch GMM helpers preserve the caller's R RNG state", {
  skip_if_not_installed("torch")
  set.seed(8675309)
  before <- .Random.seed
  invisible(torch_gmm_sample_rows(matrix(seq_len(40), nrow = 20), n = 5, seed = 2))
  expect_identical(.Random.seed, before)

  x <- rbind(
    matrix(c(-2, -2, -2.1, -1.9, -1.9, -2.1), ncol = 2, byrow = TRUE),
    matrix(c(2, 2, 2.1, 1.9, 1.9, 2.1), ncol = 2, byrow = TRUE)
  )
  invisible(torch_gmm_fit(
    x, 2L, batch_size = 3L, max_iter = 2L,
    device = "cpu", seed = 3L, num_threads = 1L
  ))
  expect_identical(.Random.seed, before)
})

test_that("Torch GMM calls restore the caller's thread setting", {
  skip_if_not_installed("torch")
  before <- torch::torch_get_num_threads()
  requested <- if (before == 1L) 2L else 1L
  invisible(torch_gmm_runtime_info(
    device = "cpu", batch_size = 4L, seed = 1L, num_threads = requested
  ))
  expect_identical(torch::torch_get_num_threads(), before)

  x <- rbind(c(-2, -2), c(-2.1, -1.9), c(-1.9, -2.1),
             c(2, 2), c(2.1, 1.9), c(1.9, 2.1))
  invisible(torch_gmm_fit(
    x, 2L, batch_size = 3L, max_iter = 2L,
    device = "cpu", seed = 4L, num_threads = requested
  ))
  expect_identical(torch::torch_get_num_threads(), before)
})
