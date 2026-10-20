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

test_that("compiled multi-fit progress is optional and numerically inert", {
  x <- cpp_gmm_fixture()
  arguments <- list(
    x = x, n_clusters = 2:4, seeds = 101:103, workers = 3L,
    threads_per_fit = 1L, max_iter = 12L, tolerance = 1e-6
  )

  expect_output(
    quiet <- do.call(cpp_gmm_fit_many, c(arguments, list(progress = FALSE))),
    NA
  )
  expect_output(
    visible <- do.call(cpp_gmm_fit_many, c(arguments, list(progress = TRUE))),
    "3/3 jobs finished",
    fixed = TRUE
  )
  for (i in seq_along(quiet)) {
    expect_identical(visible[[i]]$labels, quiet[[i]]$labels)
    expect_equal(visible[[i]]$nll, quiet[[i]]$nll, tolerance = 0)
    expect_equal(visible[[i]]$means, quiet[[i]]$means, tolerance = 0)
    expect_equal(visible[[i]]$covariances, quiet[[i]]$covariances, tolerance = 0)
  }
  expect_error(
    do.call(cpp_gmm_fit_many, c(arguments, list(progress = NA))),
    "progress"
  )
})

test_that("compiled progress path joins all workers after native fit failures", {
  x <- matrix(rep(c(-1e308, 1e308), 100L), nrow = 100L, ncol = 2L)
  expect_output(
    expect_error(
      cpp_gmm_fit_many(
        x, 2:4, seeds = 1:3, workers = 3L, threads_per_fit = 1L,
        max_iter = 2L, progress = TRUE
      ),
      "Native GMM job 1 failed.*Cholesky factorization failed"
    ),
    "3/3 jobs finished",
    fixed = TRUE
  )
})

test_that("compiled GMM stored labels, likelihood, and parameters are coherent", {
  x <- cpp_gmm_fixture()
  ridge <- 1e-5
  model <- cpp_gmm_fit(x, 3L, max_iter = 50L, tolerance = 1e-8,
                       covariance_regularization = ridge, seed = 27L)

  expect_s3_class(model, "cpp_gmm")
  expect_s3_class(model, "torch_gmm")
  expect_identical(model$labels, cpp_gmm_predict(model, x, batch_size = 17L))
  expect_identical(model$labels, predict(model, x, batch_size = 17L))
  probability <- cpp_gmm_predict_proba(model, x, batch_size = 13L)
  expect_identical(
    cpp_gmm_predict(model, x, batch_size = 11L),
    max.col(probability, ties.method = "first") - 1L
  )
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

test_that("compiled label prediction is genuinely batch bounded", {
  means <- cbind(seq(-7, 7), seq(-3.5, 3.5, length.out = 15L))
  model <- structure(list(
    means = means,
    covariances = rep(list(diag(2L)), nrow(means)),
    weights = rep(1 / nrow(means), nrow(means)),
    n_features = 2L,
    n_clusters = nrow(means),
    covariance_regularization = 1e-6
  ), class = c("cpp_gmm", "torch_gmm"))
  set.seed(920)
  x <- matrix(stats::rnorm(40002L), ncol = 2L)

  observed_batch_rows <- integer()
  original_scores <- CellCharterR:::.cpp_gmm_log_scores
  testthat::local_mocked_bindings(
    .cpp_gmm_log_scores = function(x, prepared) {
      observed_batch_rows <<- c(observed_batch_rows, nrow(x))
      original_scores(x, prepared)
    },
    cpp_gmm_predict_proba = function(...) {
      stop("label prediction attempted to construct the full posterior")
    },
    .package = "CellCharterR"
  )

  labels <- cpp_gmm_predict(model, x, batch_size = 37L)
  expect_length(labels, nrow(x))
  expect_true(all(labels %in% 0:14))
  expect_lte(max(observed_batch_rows), 37L)
  expect_gt(length(observed_batch_rows), 1L)
})

test_that("compiled scoring reads DelayedMatrix inputs in bounded row blocks", {
  means <- rbind(c(-2, 0), c(0, 2), c(2, -1))
  model <- structure(list(
    means = means,
    covariances = list(diag(2L), diag(c(1.5, 0.75)), diag(c(0.5, 2))),
    weights = c(0.2, 0.5, 0.3),
    n_features = 2L,
    n_clusters = 3L,
    covariance_regularization = 1e-6
  ), class = c("cpp_gmm", "torch_gmm"))
  set.seed(921)
  x <- matrix(stats::rnorm(514L), ncol = 2L)
  delayed <- DelayedArray::DelayedArray(x)

  expected_labels <- cpp_gmm_predict(model, x, batch_size = 23L)
  expected_probability <- cpp_gmm_predict_proba(model, x, batch_size = 23L)
  expected_likelihood <- cpp_gmm_log_likelihood(model, x, batch_size = 23L)
  observed_batch_rows <- integer()
  original_get_block <- CellCharterR:::.cpp_gmm_get_block
  testthat::local_mocked_bindings(
    .cpp_gmm_get_block = function(reader, index) {
      observed_batch_rows <<- c(observed_batch_rows, length(index))
      original_get_block(reader, index)
    },
    .package = "CellCharterR"
  )

  expect_identical(
    cpp_gmm_predict(model, delayed, batch_size = 23L), expected_labels
  )
  expect_equal(
    cpp_gmm_predict_proba(model, delayed, batch_size = 23L),
    expected_probability,
    tolerance = 0
  )
  expect_equal(
    cpp_gmm_log_likelihood(model, delayed, batch_size = 23L),
    expected_likelihood,
    tolerance = 0
  )
  expect_lte(max(observed_batch_rows), 23L)
  expect_identical(
    length(observed_batch_rows),
    as.integer(3L * ceiling(nrow(x) / 23L))
  )

  delayed_bad <- DelayedArray::DelayedArray(rbind(x[1:3, ], c(NA_real_, 0)))
  expect_error(
    cpp_gmm_predict(model, delayed_bad, batch_size = 2L),
    "finite numeric matrix blocks"
  )
})

test_that("generic prediction preserves double precision for compiled models", {
  model <- structure(list(
    means = matrix(c(1e8, 1e8 + 4), ncol = 1L),
    covariances = list(matrix(1), matrix(1)),
    weights = c(0.5, 0.5),
    n_features = 1L,
    n_clusters = 2L,
    covariance_regularization = 1e-8
  ), class = c("cpp_gmm", "torch_gmm"))
  x <- matrix(c(1e8 + 0.25, 1e8 + 3.75), ncol = 1L)

  testthat::local_mocked_bindings(
    torch_gmm_predict = function(...) stop("Torch prediction was dispatched"),
    .package = "CellCharterR"
  )
  expect_identical(predict(model, x, batch_size = 1L), c(0L, 1L))
  expect_identical(
    predict(model, x, batch_size = 1L),
    cpp_gmm_predict(model, x, batch_size = 2L)
  )
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
    artifact <- fit$embedded_best_artifacts[[as.character(k)]]
    expect_s3_class(model, "torch_gmm")
    expect_identical(model$labels, cpp_gmm_predict(model, x, batch_size = 23L))
    expect_equal(model$selection_nll, utils::tail(model$nll_history, 1L),
                 tolerance = 0)
    expect_equal(artifact$nll, model$selection_nll, tolerance = 0)
    expect_equal(artifact$final_nll, model$nll, tolerance = 0)
    expect_equal(fit$best_index[as.character(k), "nll"],
                 model$selection_nll, tolerance = 0)
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

test_that("incompatible C++ campaigns fail closed without overwriting artifacts", {
  x <- cpp_gmm_fixture()
  checkpoint_dir <- tempfile("cpp-autok-fail-closed-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  arguments <- list(
    x = x, n_clusters = 2:3, checkpoint_dir = checkpoint_dir,
    max_runs = 1L, convergence_tol = 0, seed = 59L,
    backend = "cpp", workers = 2L, checkpoint_mode = "round",
    fit_args = list(max_iter = 2L, tolerance = 1e-6), progress = FALSE
  )
  invisible(do.call(cellcharter_autok_fit, arguments))
  files <- list.files(checkpoint_dir, recursive = TRUE, full.names = TRUE)
  before <- vapply(files, .autok_file_sha256, character(1L))

  arguments$x[[1L]] <- arguments$x[[1L]] + 0.125
  expect_error(
    do.call(cellcharter_autok_fit, arguments),
    "Existing artifacts were not modified"
  )
  expect_identical(
    unname(vapply(files, .autok_file_sha256, character(1L))),
    unname(before)
  )
})

test_that("C++ round checkpoints survive an observability-only wrapper upgrade", {
  x <- cpp_gmm_fixture()
  checkpoint_dir <- tempfile("cpp-autok-progress-upgrade-")
  arguments <- list(
    x = x, n_clusters = 2:3, checkpoint_dir = checkpoint_dir,
    max_runs = 2L, convergence_tol = 0, seed = 63L,
    backend = "cpp", workers = 2L, checkpoint_mode = "round",
    fit_args = list(max_iter = 20L, tolerance = 1e-6), progress = FALSE
  )
  original <- do.call(cellcharter_autok_fit, arguments)

  state_path <- file.path(checkpoint_dir, "autok_state.rds")
  state <- readRDS(state_path)
  previous_contract <- state$resume_contract
  previous_contract$algorithm <- paste0(
    "cpp-full-covariance-em-v1.3.0:", paste(rep("a", 64L), collapse = "")
  )
  previous_id <- digest::digest(
    .cc_contract_canonical(previous_contract), algo = "sha256", serialize = TRUE
  )
  checkpoint_paths <- unique(state$fit_index$path)
  for (path in checkpoint_paths) {
    bundle <- readRDS(path)
    bundle$data_id <- previous_id
    bundle$artifacts <- lapply(bundle$artifacts, function(artifact) {
      artifact$data_id <- previous_id
      artifact
    })
    saveRDS(bundle, path)
  }
  state$resume_contract <- previous_contract
  state$data_id <- previous_id
  state$fit_index$bytes <- as.numeric(file.info(state$fit_index$path)$size)
  # This fixture emulates a checkpoint written before cryptographic artifact
  # guards were introduced; the resume path must validate and add them once.
  state$fit_index$sha256 <- NULL
  saveRDS(state, state_path)

  expect_message(
    resumed <- do.call(cellcharter_autok_fit, arguments),
    "observability-only package upgrade"
  )
  expect_true(all(resumed$fit_index$resumed))
  expect_identical(resumed$best_k, original$best_k)
  expect_equal(resumed$stability, original$stability, tolerance = 0)

  resumed_again <- do.call(cellcharter_autok_fit, arguments)
  expect_true(all(resumed_again$fit_index$resumed))
  expect_identical(resumed_again$best_k, original$best_k)
})

test_that("SPE clustering adapter attaches labels from C++ memory and round artifacts", {
  for (mode in c("none", "round")) {
    spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
    checkpoint_dir <- tempfile(paste0("cpp-sce-", mode, "-"))
    spe <- cellcharter_cluster(
      spe,
      representation = "X_cellcharter_tiny",
      candidate_k = c(2L, 2L),
      checkpoint_dir = checkpoint_dir,
      domain_col = paste0("cpp_domain_", mode),
      max_runs = 2L,
      convergence_tol = 0,
      seed = 71L,
      backend = "cpp",
      workers = 2L,
      threads_per_fit = 1L,
      checkpoint_mode = mode,
      resume = FALSE,
      fit_args = list(max_iter = 10L, tolerance = 1e-5),
      progress = FALSE
    )

    expect_s4_class(spe, "SpatialExperiment")
    expect_false(anyNA(spe[[paste0("cpp_domain_", mode)]]))
    expect_length(spe[[paste0("cpp_domain_", mode)]], ncol(spe))
    stage <- S4Vectors::metadata(spe)$CellCharterR$stages[[
      paste0("clustering:cpp_domain_", mode)
    ]]
    if (identical(mode, "none")) {
      expect_null(stage$parameters$checkpoint_dir)
      expect_null(stage$parameters$selected_artifact)
      expect_false(dir.exists(checkpoint_dir))
    } else {
      expect_true(dir.exists(checkpoint_dir))
      expect_true(file.exists(stage$parameters$selected_artifact))
    }
  }
})
