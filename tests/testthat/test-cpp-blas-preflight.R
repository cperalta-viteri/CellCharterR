test_that("parallel compiled GMM fails before materialisation with unsafe BLAS", {
  withr::local_envvar(c(
    OPENBLAS_NUM_THREADS = "2",
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    GOTO_NUM_THREADS = "1"
  ))

  expect_error(
    cpp_gmm_fit_many(
      matrix(0, nrow = 12L, ncol = 2L),
      n_clusters = 2L, workers = 2L, threads_per_fit = 1L,
      max_iter = 1L, progress = FALSE
    ),
    "single-threaded BLAS"
  )
})

test_that("safe BLAS controls are recorded in the C++ execution contract", {
  controls <- c(
    OPENBLAS_NUM_THREADS = "1",
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    GOTO_NUM_THREADS = "1"
  )
  withr::local_envvar(controls)

  contract <- CellCharterR:::.autok_execution_contract(
    "cpp", list(), list(), FALSE, 4L
  )
  expect_identical(contract$schema_version, 2L)
  expect_identical(unlist(contract$blas$environment, use.names = TRUE), controls)
  expect_true(is.character(contract$blas$library))
  expect_length(contract$blas$library, 1L)
})

test_that("AutoK rejects unsafe parallel BLAS before fingerprinting", {
  withr::local_envvar(c(
    OPENBLAS_NUM_THREADS = "32",
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    GOTO_NUM_THREADS = "1"
  ))
  fingerprint_called <- FALSE
  testthat::local_mocked_bindings(
    .autok_representation_sha256 = function(...) {
      fingerprint_called <<- TRUE
      stop("fingerprint should not run")
    },
    .package = "CellCharterR"
  )

  expect_error(
    cellcharter_autok_fit(
      matrix(0, nrow = 12L, ncol = 2L),
      n_clusters = c(2L, 3L),
      checkpoint_dir = tempfile("autok-blas-preflight-"),
      max_runs = 1L,
      backend = "cpp",
      workers = 2L,
      threads_per_fit = 1L,
      checkpoint_mode = "none",
      progress = FALSE
    ),
    "single-threaded BLAS"
  )
  expect_false(fingerprint_called)
})
