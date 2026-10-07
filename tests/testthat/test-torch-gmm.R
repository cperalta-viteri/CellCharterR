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

test_that("Torch GMM auto-device smoke test returns finite likelihood", {
  skip_if_not_installed("torch")
  set.seed(19)
  x <- rbind(matrix(rnorm(60, -2), ncol = 3), matrix(rnorm(60, 2), ncol = 3))
  fit <- torch_gmm_fit(x, 2, batch_size = 16, max_iter = 8,
                       device = "auto", seed = 19, num_threads = 1)
  labels <- torch_gmm_predict(fit, x, batch_size = 13, device = "auto")
  expect_true(is.finite(fit$nll))
  expect_true(is.finite(torch_gmm_log_likelihood(fit, x, 17, "auto")))
  expect_length(labels, nrow(x))
  expect_true(all(labels %in% 0:1))
  expect_identical(fit$runtime$cuda_available, torch::cuda_is_available())
  if (torch::cuda_is_available()) expect_identical(fit$runtime$device, "cuda")
})

test_that("R prediction matches the pinned PyCave K=11 fixture", {
  skip_if_not_installed("torch")
  path <- system.file(
    "extdata", "python_gmm_k11_fixture.rds", package = "rCellCharter",
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
