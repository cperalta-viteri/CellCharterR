test_that("numeric AutoK input is unambiguously a pair of selectable bounds", {
  expect_identical(
    .autok_k_grid(c(2L, 20L)),
    list(fit_k = 1L:21L, requested_k = 2L:20L)
  )
  expect_identical(
    .autok_k_grid(as.list(1L:21L)),
    list(fit_k = 1L:21L, requested_k = 2L:20L)
  )
  expect_error(
    .autok_k_grid(2L:20L),
    "exactly two selectable bounds.*as.list",
    perl = TRUE
  )
})

test_that("round-mode FMI reads each prior round bundle once per pass", {
  checkpoint_dir <- tempfile("autok-round-io-")
  memory_dir <- tempfile("autok-memory-parity-")
  on.exit(unlink(c(checkpoint_dir, memory_dir), recursive = TRUE), add = TRUE)
  x <- matrix(seq_len(48), ncol = 3)

  fake_fit_many <- function(x, n_clusters, seeds, workers,
                            threads_per_fit, progress, ...) {
    stats::setNames(lapply(seq_along(n_clusters), function(i) {
      k <- n_clusters[[i]]
      run_seed <- seeds[[i]]
      list(
        labels = as.integer((seq_len(nrow(x)) - 1L + run_seed) %% k),
        nll = as.numeric(k + run_seed / 1000)
      )
    }), as.character(n_clusters))
  }

  reads <- character()
  original_read_round <- CellCharterR:::.autok_read_round_labels
  testthat::local_mocked_bindings(
    cpp_gmm_fit_many = fake_fit_many,
    .autok_read_round_labels = function(path, n_obs, k_values) {
      reads <<- c(reads, basename(path))
      original_read_round(path, n_obs, k_values)
    },
    .package = "CellCharterR"
  )

  round_fit <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
    max_runs = 3L, convergence_tol = 0, seed = 11L,
    backend = "cpp", checkpoint_mode = "round", progress = FALSE
  )
  # Run 2 compares with round 1. Run 3 compares with rounds 1 and 2. Each
  # bundle is therefore decompressed exactly once in each current-round pass.
  expect_identical(reads, c("run_001.rds", "run_001.rds", "run_002.rds"))

  memory_fit <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 3L), checkpoint_dir = memory_dir,
    max_runs = 3L, convergence_tol = 0, seed = 11L,
    backend = "cpp", checkpoint_mode = "none", progress = FALSE
  )
  expect_equal(round_fit$edge_stability, memory_fit$edge_stability, tolerance = 0)
  expect_equal(round_fit$stability, memory_fit$stability, tolerance = 0)
})

test_that("round label reader rejects incomplete bundles", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(
    format = "cellcharter_autok_round_v1",
    artifacts = list(`2` = list(labels = c(0L, 1L)))
  ), path)

  expect_error(
    .autok_read_round_labels(path, n_obs = 2L, k_values = c(2L, 3L)),
    "has no K = 3"
  )
})

test_that("AutoK records the effective compiled topology used by every fit", {
  checkpoint_dir <- tempfile("autok-effective-topology-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  x <- matrix(seq_len(48), ncol = 3)
  observed <- list()

  fake_fit_many <- function(x, n_clusters, seeds, workers,
                            threads_per_fit, progress, ...) {
    observed[[length(observed) + 1L]] <<- list(
      workers = workers,
      threads_per_fit = threads_per_fit
    )
    stats::setNames(lapply(seq_along(n_clusters), function(i) {
      k <- n_clusters[[i]]
      list(
        labels = as.integer((seq_len(nrow(x)) - 1L + seeds[[i]]) %% k),
        nll = as.numeric(k),
        selection_nll = as.numeric(k)
      )
    }), as.character(n_clusters))
  }

  testthat::local_mocked_bindings(
    cpp_gmm_fit_many = fake_fit_many,
    .cpp_gmm_require_single_thread_blas = function(...) invisible(list()),
    .package = "CellCharterR"
  )

  expect_warning(
    fit <- cellcharter_autok_fit(
      x, n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
      max_runs = 2L, convergence_tol = 0, seed = 11L,
      backend = "cpp", workers = 20L, threads_per_fit = 2L,
      checkpoint_mode = "round", progress = FALSE
    ),
    "Capping `threads_per_fit` at 1"
  )

  expect_length(observed, 2L)
  expect_true(all(vapply(observed, `[[`, integer(1), "workers") == 20L))
  expect_true(all(vapply(observed, `[[`, integer(1), "threads_per_fit") == 1L))
  expect_identical(fit$workers, 20L)
  expect_identical(fit$threads_per_fit, 1L)
  expect_identical(fit$requested_workers, 20L)
  expect_identical(fit$requested_threads_per_fit, 2L)
  expect_identical(fit$execution_contract$threads_per_fit, 1L)

  state <- readRDS(file.path(checkpoint_dir, "autok_state.rds"))
  expect_identical(state$workers, 20L)
  expect_identical(state$threads_per_fit, 1L)
  expect_identical(state$requested_workers, 20L)
  expect_identical(state$requested_threads_per_fit, 2L)
  expect_identical(state$execution_contract$threads_per_fit, 1L)

  bundle <- readRDS(CellCharterR:::.autok_round_path(checkpoint_dir, 1L))
  expect_true(all(vapply(
    bundle$artifacts,
    function(artifact) identical(artifact$threads_per_fit, 1L),
    logical(1)
  )))
  expect_true(all(vapply(
    bundle$artifacts,
    function(artifact) identical(artifact$requested_threads_per_fit, 2L),
    logical(1)
  )))
})
