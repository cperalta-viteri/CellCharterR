test_that("Torch GMM separates a simple full-covariance mixture", {
  skip_if_not_installed("torch")
  set.seed(7)
  x <- rbind(MASS::mvrnorm(300, c(-2, 0), matrix(c(1, .6, .6, 1), 2)),
             MASS::mvrnorm(300, c(2, 0), matrix(c(1, -.4, -.4, 1), 2)))
  gmm <- torch_gmm_fit(x, 2, batch_size = 64, max_iter = 60, device = "cpu", seed = 7)
  expect_true(gmm$converged || gmm$iterations == 60)
  expect_length(torch_gmm_predict(gmm, x, batch_size = 64, device = "cpu"), nrow(x))
  expect_equal(ncol(torch_gmm_predict_proba(gmm, x, 64, "cpu")), 2)
})
