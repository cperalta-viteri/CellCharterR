test_that("RAM representation cache preserves values, dimensions, and row order", {
  skip_if_not_installed("torch")
  x <- outer(seq_len(7L), seq_len(5L), function(i, j) 1000 * i + j)
  expected <- torch::as_array(torch::torch_tensor(x, dtype = torch::torch_float32()))

  for (batch_size in c(1L, 2L, 3L, 7L, 10L)) {
    cache <- .rc_ram_cache(x, batch_size)
    expect_identical(dim(cache$tensor), c(7L, 5L))
    expect_equal(torch::as_array(cache$tensor), expected, tolerance = 0)
    observed <- .rc_get(.rc_reader(cache), c(7L, 2L, 5L))
    expect_equal(observed, expected[c(7L, 2L, 5L), , drop = FALSE], tolerance = 0)
  }
})

test_that("RAM representation cache rejects overflow and malformed objects", {
  skip_if_not_installed("torch")
  expect_error(.rc_ram_cache(matrix(1e40, nrow = 1L), 1L), "float32 range")

  malformed <- structure(
    list(format_version = 2L, validated = TRUE, n = 2L, p = 2L,
         tensor = torch::torch_zeros(c(1L, 2L))),
    class = "cellcharter_torch_cache"
  )
  expect_error(.rc_reader(malformed), "Malformed")
})

test_that("RAM cache and streamed readers give the same small K=1 fit", {
  skip_if_not_installed("torch")
  set.seed(19)
  x <- matrix(rnorm(240), ncol = 4L)
  cache <- .rc_ram_cache(x, 13L)
  args <- list(
    n_clusters = 1L, batch_size = 13L, max_iter = 2L, device = "cpu",
    seed = 23L, num_threads = 1L, initialization_method = "torch"
  )
  streamed <- do.call(torch_gmm_fit, c(list(x = x), args))
  cached <- do.call(torch_gmm_fit, c(list(x = cache), args))
  expect_equal(cached$nll_history, streamed$nll_history, tolerance = 1e-6)
  expect_equal(cached$nll, streamed$nll, tolerance = 1e-6)
  expect_equal(cached$means, streamed$means, tolerance = 1e-6)
  expect_equal(cached$covariances, streamed$covariances, tolerance = 1e-6)
  expect_identical(
    torch_gmm_predict(cached, cache, batch_size = 13L, device = "cpu"),
    torch_gmm_predict(streamed, x, batch_size = 13L, device = "cpu")
  )
})

test_that("RAM caching does not invalidate compatible AutoK checkpoints", {
  skip_if_not_installed("torch")
  checkpoint <- tempfile("autok-ram-cache-")
  x <- matrix(c(-2, -2, -1.8, -2.1, 2, 2, 2.1, 1.8), ncol = 2L, byrow = TRUE)
  common <- list(
    x = x, n_clusters = c(2L, 2L), checkpoint_dir = checkpoint,
    max_runs = 1L, seed = 31L,
    fit_args = list(batch_size = 4L, max_iter = 2L, device = "cpu",
                    num_threads = 1L, initialization_method = "torch"),
    predict_args = list(batch_size = 4L, device = "cpu"),
    data_id = "cache-resume-fixture"
  )
  first <- do.call(cellcharter_autok_fit, c(common, list(cache_representation = "ram")))
  resumed <- do.call(cellcharter_autok_fit, c(common, list(cache_representation = "none")))
  expect_identical(first$data_id, resumed$data_id)
  expect_true(all(resumed$fit_index$resumed))
})
