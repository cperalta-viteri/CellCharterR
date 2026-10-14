test_that("CPU autotuning has deterministic explicit behavior", {
  first <- torch_gmm_autotune_batch_size(
    30L, 12L, n_obs = 5000L, device = "cpu", cpu_batch_size = 2048L
  )
  second <- torch_gmm_autotune_batch_size(
    30L, 12L, n_obs = 5000L, device = "cpu", cpu_batch_size = 2048L
  )
  expect_identical(first, second)
  expect_identical(first$batch_size, 2048L)
  expect_identical(first$memory_source, "cpu_default")
})

test_that("CUDA estimate is capped, deterministic, and allocation-free", {
  skip_if_not_installed("torch")
  cuda <- tryCatch(isTRUE(torch::cuda_is_available()), error = function(e) FALSE)
  skip_if_not(cuda)
  args <- list(
    n_features = 30L, n_clusters = 12L, n_obs = 4000000L,
    device = "cuda", cap = 65536L, safety_fraction = 0.5,
    free_memory_bytes = 8 * 1024^3
  )
  first <- do.call(torch_gmm_autotune_batch_size, args)
  second <- do.call(torch_gmm_autotune_batch_size, args)
  expect_identical(first, second)
  expect_identical(first$batch_size, 65536L)
  expect_identical(first$memory_source, "argument")
  expect_gt(first$estimated_bytes_per_observation, 0)
})

test_that("pure CUDA memory estimate is conservative and bounded", {
  small <- .cc_gmm_batch_estimate(30L, 12L, 512 * 1024^2, 0.5, 100000L)
  large <- .cc_gmm_batch_estimate(30L, 12L, 1024 * 1024^2, 0.5, 100000L)
  capped <- .cc_gmm_batch_estimate(30L, 12L, 8 * 1024^3, 0.5, 4096L)
  expect_gte(large$batch_size, small$batch_size)
  expect_identical(capped$batch_size, 4096L)
  expect_gt(small$estimated_bytes_per_observation, 0)
})

test_that("autotuner validates dimensions and memory controls", {
  expect_error(torch_gmm_autotune_batch_size(0, 2, device = "cpu"), "n_features")
  expect_error(torch_gmm_autotune_batch_size(2, 0, device = "cpu"), "n_clusters")
  expect_error(torch_gmm_autotune_batch_size(2, 2, device = "cpu", cap = 1.5), "cap")
  expect_error(
    torch_gmm_autotune_batch_size(2, 2, device = "cpu", safety_fraction = 1.1),
    "safety_fraction"
  )
  expect_error(
    torch_gmm_autotune_batch_size(2, 2, device = "cpu", free_memory_bytes = -1),
    "free_memory_bytes"
  )
})
