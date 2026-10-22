test_that("Fowlkes-Mallows is label-invariant and handles singleton partitions", {
  a <- c(0, 0, 1, 1, 2, 2)
  expect_equal(cellcharter_fowlkes_mallows(a, c("x", "x", "z", "z", "q", "q")), 1)
  expect_equal(cellcharter_fowlkes_mallows(a, seq_along(a)), 0)
  expect_error(cellcharter_fowlkes_mallows(a, a[-1]), "same length")
  expect_error(cellcharter_fowlkes_mallows(a, replace(a, 1, NA)), "missing")
})

test_that("compact integer FMI matches the general contingency calculation", {
  set.seed(814)
  labels_a <- sample.int(11L, 50000L, replace = TRUE) - 1L
  labels_b <- sample.int(13L, 50000L, replace = TRUE) - 1L
  expected <- cellcharter_fowlkes_mallows(
    as.character(labels_a), as.character(labels_b)
  )
  expect_equal(
    cellcharter_fowlkes_mallows(labels_a, labels_b),
    expected,
    tolerance = 0
  )
})

test_that("sparse integer FMI labels avoid a range-sized contingency allocation", {
  labels_a <- c(0L, 100000000L, 0L, 100000000L)
  labels_b <- c(0L, 0L, 1L, 1L)
  expect_equal(
    cellcharter_fowlkes_mallows(labels_a, labels_b),
    cellcharter_fowlkes_mallows(as.character(labels_a), as.character(labels_b)),
    tolerance = 0
  )
})

test_that("regularization retry diagnostics retain one AutoK index row per fit", {
  checkpoint_dir <- tempfile("autok-retry-diagnostics-")
  x <- matrix(seq_len(30), ncol = 3)
  calls <- 0L
  fit_hook <- function(x, n_clusters, seed, ...) {
    calls <<- calls + 1L
    list(
      k = n_clusters, nll = as.numeric(n_clusters),
      covariance_regularization = 1e-5,
      regularization_attempts = c(1e-6, 1e-5)
    )
  }
  predict_hook <- function(model, newdata, ...) {
    (seq_len(nrow(newdata)) - 1L) %% model$k
  }

  first <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
    max_runs = 1L, seed = 17L, fit_fun = fit_hook,
    predict_fun = predict_hook, nll_fun = function(model) model$nll,
    data_id = "retry-diagnostics-v1"
  )
  expect_equal(nrow(first$fit_index), length(first$fit_k))
  expect_false(anyDuplicated(first$fit_index[c("k", "run")]) > 0L)
  expect_identical(first$fit_index$regularization_attempts, rep.int(2L, length(first$fit_k)))

  resumed <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
    max_runs = 1L, seed = 17L, fit_fun = fit_hook,
    predict_fun = predict_hook, nll_fun = function(model) model$nll,
    data_id = "retry-diagnostics-v1", resume = TRUE
  )
  expect_equal(calls, length(first$fit_k))
  expect_true(all(resumed$fit_index$resumed))
  expect_false(anyDuplicated(resumed$fit_index[c("k", "run")]) > 0L)
})

test_that("legacy ridge-valued retry diagnostics migrate to an integer count", {
  checkpoint_dir <- tempfile("autok-legacy-retry-value-")
  x <- matrix(seq_len(30), ncol = 3)
  fit_hook <- function(x, n_clusters, seed, ...) {
    list(
      k = n_clusters, nll = as.numeric(n_clusters),
      covariance_regularization = 1e-6,
      regularization_attempts = 1e-6
    )
  }
  predict_hook <- function(model, newdata, ...) {
    (seq_len(nrow(newdata)) - 1L) %% model$k
  }
  first <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
    max_runs = 1L, seed = 23L, fit_fun = fit_hook,
    predict_fun = predict_hook, nll_fun = function(model) model$nll,
    data_id = "legacy-retry-value-v1"
  )
  state_path <- file.path(checkpoint_dir, "autok_state.rds")
  state <- readRDS(state_path)
  state$fit_index$regularization_attempts[] <- 1e-6
  saveRDS(state, state_path)

  resumed <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
    max_runs = 1L, seed = 23L, fit_fun = fit_hook,
    predict_fun = predict_hook, nll_fun = function(model) model$nll,
    data_id = "legacy-retry-value-v1", resume = TRUE
  )
  expect_identical(
    resumed$fit_index$regularization_attempts,
    rep.int(1L, length(first$fit_k))
  )
  expect_true(all(resumed$fit_index$resumed))
})

test_that("neighboring-K edge scores are mirrored like CellCharter 0.2.0", {
  edges <- data.frame(
    run = rep(2L, 3), previous_run = rep(1L, 3),
    lower_k = 1:3, upper_k = 2:4, score = c(0.1, 0.2, 0.3)
  )
  mirrored <- .autok_mirror_stability(edges, 1:4)
  expect_equal(unname(mirrored), rbind(c(0.2, 0.1), c(0.3, 0.2)))
  expect_equal(rownames(mirrored), c("2", "3"))
})

test_that("AutoK selectable bounds start at K equals two", {
  expect_error(
    .autok_k_grid(c(1L, 4L)),
    "lower bound of at least 2"
  )
  expect_identical(.autok_k_grid(c(2L, 4L))$fit_k, 1:5)
  expect_identical(.autok_k_grid(as.list(1L:5L))$requested_k, 2:4)
})

test_that("AutoK converges, selects lowest NLL, and resumes atomically", {
  checkpoint_dir <- tempfile("autok-")
  x <- matrix(seq_len(60), ncol = 3)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L

  fit_hook <- function(x, n_clusters, seed, ...) {
    calls$n <- calls$n + 1L
    list(k = n_clusters, seed = seed, run = seed - 40L, nll = 100 - (seed - 40L))
  }
  predict_hook <- function(model, newdata, ...) {
    (seq_len(nrow(newdata)) - 1L) %% model$k
  }

  fit <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 4L), checkpoint_dir = checkpoint_dir,
    max_runs = 10L, convergence_tol = 0.01, seed = 41L,
    fit_fun = fit_hook, predict_fun = predict_hook,
    nll_fun = function(model) model$nll, data_id = "fixture-v1"
  )

  # Identical labels across repetitions make the stability curve stop at the
  # first permitted check, after three complete runs.
  expect_true(fit$converged)
  expect_equal(fit$runs_completed, 3L)
  expect_equal(fit$fit_k, 1:5)
  expect_equal(fit$candidate_k, 2:4)
  expect_equal(calls$n, 15L)
  expect_equal(nrow(fit$edge_stability), choose(3, 2) * 4)
  expect_equal(dim(fit$stability), c(3L, 2L * choose(3, 2)))
  expect_equal(fit$best_k, cellcharter_autok_best_k(fit))
  expect_true(all(fit$best_index$run == 3L))
  expect_equal(cellcharter_autok_model(fit, 3L)$run, 3L)
  expect_length(list.files(checkpoint_dir, pattern = "[.]tmp$", recursive = TRUE), 0L)

  resumed <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 4L), checkpoint_dir = checkpoint_dir,
    max_runs = 10L, convergence_tol = 0.01, seed = 41L,
    fit_fun = fit_hook, predict_fun = predict_hook,
    nll_fun = function(model) model$nll, data_id = "fixture-v1", resume = TRUE
  )
  expect_equal(resumed$runs_completed, 3L)
  expect_true(all(resumed$fit_index$resumed))
  expect_equal(resumed$stability, fit$stability)
})

test_that("AutoK fit identity includes implementation and fit/predict parameters", {
  checkpoint_dir <- tempfile("autok-contract-")
  x <- matrix(seq_len(24), ncol = 3)
  fit_hook <- function(x, n_clusters, seed, scale = 1) {
    list(k = n_clusters, nll = n_clusters * scale)
  }
  predict_hook <- function(model, newdata, offset = 0L) {
    ((seq_len(nrow(newdata)) - 1L + offset) %% model$k)
  }
  first <- cellcharter_autok_fit(
    x, c(2L, 3L), checkpoint_dir, max_runs = 1L,
    fit_fun = fit_hook, predict_fun = predict_hook,
    nll_fun = function(model) model$nll,
    fit_args = list(scale = 1), predict_args = list(offset = 0L),
    data_id = "same-input"
  )
  expect_identical(first$input_data_id, "same-input")
  expect_match(first$data_id, "^[0-9a-f]{64}$")

  expect_error(
    cellcharter_autok_fit(
      x, c(2L, 3L), checkpoint_dir, max_runs = 1L,
      fit_fun = fit_hook, predict_fun = predict_hook,
      nll_fun = function(model) model$nll,
      fit_args = list(scale = 2), predict_args = list(offset = 0L),
      data_id = "same-input"
    ),
    "does not match current `data_id`"
  )
})

test_that("default AutoK checkpoints survive compatible implementation upgrades", {
  skip_if_not_installed("torch")
  checkpoint_dir <- tempfile("autok-transitive-")
  x <- matrix(c(-2, -2, -1.8, -2.1, 2, 2, 2.1, 1.8), ncol = 2, byrow = TRUE)
  first <- cellcharter_autok_fit(
    x, c(2L, 3L), checkpoint_dir, max_runs = 1L, seed = 9L,
    fit_args = list(batch_size = 4L, max_iter = 2L, device = "cpu",
                    num_threads = 1L),
    predict_args = list(batch_size = 4L, device = "cpu"),
    data_id = "same-input"
  )

  original <- CellCharterR:::.rc_logprob
  testthat::local_mocked_bindings(
    .rc_logprob = function(x, means, chol, weights) {
      original(x, means, chol, weights)
    },
    .package = "CellCharterR"
  )
  resumed <- cellcharter_autok_fit(
      x, c(2L, 3L), checkpoint_dir, max_runs = 1L, seed = 9L,
      fit_args = list(batch_size = 4L, max_iter = 2L, device = "cpu",
                      num_threads = 1L),
      predict_args = list(batch_size = 4L, device = "cpu"),
      data_id = "same-input"
  )
  expect_true(all(resumed$fit_index$resumed))
  expect_identical(resumed$data_id, first$data_id)
  expect_match(first$implementation_digest, "^[0-9a-f]{64}$")
  expect_error(
    cellcharter_autok_fit(
      x, c(2L, 3L), checkpoint_dir, max_runs = 1L, seed = 9L,
      fit_args = list(batch_size = 2L, max_iter = 2L, device = "cpu",
                      num_threads = 1L),
      predict_args = list(batch_size = 4L, device = "cpu"),
      data_id = "same-input"
    ),
    "numerical execution contract"
  )
})

test_that("default AutoK resume identity ignores runtime-only controls", {
  first <- .autok_resume_contract(
    "input", list(max_iter = 10L, batch_size = 4096L, device = "cuda", num_threads = 2L),
    list(batch_size = 4096L, device = "cuda"), TRUE, "old"
  )
  second <- .autok_resume_contract(
    "input", list(max_iter = 10L, batch_size = 32768L, device = "cpu", num_threads = 8L),
    list(batch_size = 65536L, device = "cpu"), TRUE, "new"
  )
  expect_identical(first, second)
  changed <- second
  changed$fit_args$max_iter <- 11L
  expect_false(identical(first, changed))

  torch_a <- .autok_execution_contract(
    "torch", list(batch_size = 4096L, device = "cpu", num_threads = 1L),
    list(batch_size = 4096L, device = "cpu"), TRUE, 1L
  )
  torch_b <- .autok_execution_contract(
    "torch", list(batch_size = 8192L, device = "cpu", num_threads = 2L),
    list(batch_size = 4096L, device = "cpu"), TRUE, 1L
  )
  expect_false(identical(torch_a, torch_b))
  expect_false(identical(
    .autok_execution_contract("cpp", list(), list(), FALSE, 1L),
    .autok_execution_contract("cpp", list(), list(), FALSE, 2L)
  ))
})

test_that("Torch legacy migration rejects C++ checkpoint identities", {
  checkpoint_dir <- tempfile("autok-cross-backend-")
  dir.create(checkpoint_dir)
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  representation_sha256 <- paste(rep.int("a", 64L), collapse = "")
  cpp_contract <- .autok_resume_contract(
    "input", list(max_iter = 10L), list(), default_torch = FALSE,
    implementation_digest = "unused",
    algorithm = "cpp-full-covariance-em-v1.3.0",
    representation_sha256 = representation_sha256
  )
  cpp_id <- digest::digest(
    .cc_contract_canonical(cpp_contract), algo = "sha256", serialize = TRUE
  )
  state <- list(
    data_id = cpp_id, compatible_data_ids = character(),
    input_data_id = "input", fit_args = list(max_iter = 10L),
    predict_args = list(), fit_k = 1:4, candidate_k = 2:3,
    seed = 42L, n_obs = 8L, n_features = 2L,
    representation_sha256 = representation_sha256,
    backend = "cpp", resume_contract = cpp_contract
  )
  saveRDS(state, file.path(checkpoint_dir, "autok_state.rds"))
  current <- state[c(
    "input_data_id", "fit_args", "predict_args", "fit_k", "candidate_k",
    "seed", "n_obs", "n_features", "representation_sha256"
  )]
  expect_identical(
    .autok_legacy_data_id(checkpoint_dir, current), character()
  )

  torch_contract <- .autok_resume_contract(
    "input", list(max_iter = 10L), list(), default_torch = TRUE,
    implementation_digest = "unused",
    representation_sha256 = representation_sha256
  )
  state$backend <- "torch"
  state$resume_contract <- torch_contract
  state$data_id <- digest::digest(
    .cc_contract_canonical(torch_contract), algo = "sha256", serialize = TRUE
  )
  saveRDS(state, file.path(checkpoint_dir, "autok_state.rds"))
  expect_identical(
    .autok_legacy_data_id(checkpoint_dir, current), state$data_id
  )
})

test_that("AutoK duration formatter is readable", {
  expect_identical(.autok_format_duration(9), "9s")
  expect_identical(.autok_format_duration(125), "2m 05s")
  expect_identical(.autok_format_duration(9000), "2h 30m")
})

test_that("AutoK retries CUDA OOM with a smaller batch", {
  observed <- integer()
  result <- .autok_oom_retry(
    function(args) {
      observed <<- c(observed, args$batch_size)
      if (args$batch_size > 1024L) stop("CUDA out of memory")
      "ok"
    },
    list(batch_size = 4096L), stage = "test"
  )
  expect_identical(observed, c(4096L, 2048L, 1024L))
  expect_identical(result$value, "ok")
  expect_identical(result$args$batch_size, 1024L)
})

test_that("AutoK CUDA OOM retry can reach batch size one", {
  observed <- integer()
  result <- .autok_oom_retry(
    function(args) {
      observed <<- c(observed, args$batch_size)
      if (args$batch_size > 1L) stop("CUDA out of memory")
      "ok"
    },
    list(batch_size = 8L), stage = "minimum-batch-test"
  )
  expect_identical(observed, c(8L, 4L, 2L, 1L))
  expect_identical(result$args$batch_size, 1L)
})

test_that("first-fit OOM batches become the durable Torch execution contract", {
  skip_if_not_installed("torch")
  checkpoint_dir <- tempfile("autok-oom-contract-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  fake_fit <- function(x, n_clusters, seed, batch_size, ...) {
    if (batch_size > 2L) stop("CUDA out of memory")
    list(k = n_clusters, seed = seed, nll = n_clusters)
  }
  fake_predict <- function(model, newdata, batch_size, ...) {
    if (batch_size > 4L) stop("CUDA out of memory")
    rep.int(0L, nrow(newdata))
  }
  testthat::local_mocked_bindings(
    .autok_default_fit = fake_fit,
    .autok_default_predict = fake_predict,
    .package = "CellCharterR"
  )
  x <- matrix(seq_len(24), ncol = 3L)
  first <- cellcharter_autok_fit(
    x, c(2L, 2L), checkpoint_dir, max_runs = 1L,
    fit_args = list(batch_size = 8L, device = "cpu", num_threads = 1L),
    predict_args = list(batch_size = 8L, device = "cpu"),
    progress = FALSE
  )
  expect_identical(first$execution_contract$fit_batch_size, 2L)
  expect_identical(first$execution_contract$predict_batch_size, 4L)
  first_artifact <- readRDS(.autok_checkpoint_path(checkpoint_dir, 1L, 1L))
  expect_identical(first_artifact$execution_contract, first$execution_contract)

  resumed <- cellcharter_autok_fit(
    x, c(2L, 2L), checkpoint_dir, max_runs = 1L,
    fit_args = list(batch_size = 2L, device = "cpu", num_threads = 1L),
    predict_args = list(batch_size = 4L, device = "cpu"),
    progress = FALSE
  )
  expect_true(all(resumed$fit_index$resumed))
  expect_error(
    cellcharter_autok_fit(
      x, c(2L, 2L), checkpoint_dir, max_runs = 1L,
      fit_args = list(batch_size = 8L, device = "cpu", num_threads = 1L),
      predict_args = list(batch_size = 8L, device = "cpu"),
      progress = FALSE
    ),
    "numerical execution contract"
  )
})

test_that("a later OOM cannot publish a mixed-batch AutoK campaign", {
  skip_if_not_installed("torch")
  checkpoint_dir <- tempfile("autok-later-oom-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  fake_fit <- function(x, n_clusters, seed, batch_size, ...) {
    if (n_clusters == 2L && batch_size > 4L) stop("CUDA out of memory")
    list(k = n_clusters, seed = seed, nll = n_clusters)
  }
  fake_predict <- function(model, newdata, batch_size, ...) {
    rep.int(0L, nrow(newdata))
  }
  testthat::local_mocked_bindings(
    .autok_default_fit = fake_fit,
    .autok_default_predict = fake_predict,
    .package = "CellCharterR"
  )
  expect_error(
    cellcharter_autok_fit(
      matrix(seq_len(24), ncol = 3L), c(2L, 2L), checkpoint_dir,
      max_runs = 1L,
      fit_args = list(batch_size = 8L, device = "cpu", num_threads = 1L),
      predict_args = list(batch_size = 8L, device = "cpu"),
      progress = FALSE
    ),
    "changed the numerical execution contract"
  )
  expect_true(file.exists(.autok_checkpoint_path(checkpoint_dir, 1L, 1L)))
  expect_false(file.exists(.autok_checkpoint_path(checkpoint_dir, 2L, 1L)))
})

test_that("custom AutoK hooks never receive silent OOM batch retries", {
  checkpoint_dir <- tempfile("autok-custom-oom-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  fit_hook <- function(x, n_clusters, seed, batch_size, ...) {
    if (n_clusters == 2L) stop("CUDA out of memory")
    list(k = n_clusters, nll = n_clusters)
  }
  predict_hook <- function(model, newdata, batch_size, ...) {
    rep.int(0L, nrow(newdata))
  }
  expect_error(
    cellcharter_autok_fit(
      matrix(seq_len(24), ncol = 3L), c(2L, 2L), checkpoint_dir,
      max_runs = 1L, fit_fun = fit_hook, predict_fun = predict_hook,
      nll_fun = function(model) model$nll,
      fit_args = list(batch_size = 8L),
      predict_args = list(batch_size = 8L), progress = FALSE
    ),
    "CUDA out of memory"
  )
  first <- readRDS(.autok_checkpoint_path(checkpoint_dir, 1L, 1L))
  expect_identical(first$effective_fit_batch_size, 8L)
  expect_false(file.exists(.autok_checkpoint_path(checkpoint_dir, 2L, 1L)))
})

test_that("hashed legacy artifacts still undergo one-time execution validation", {
  skip_if_not_installed("torch")
  checkpoint_dir <- tempfile("autok-hashed-legacy-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  fake_fit <- function(x, n_clusters, seed, batch_size, ...) {
    list(k = n_clusters, seed = seed, nll = n_clusters)
  }
  fake_predict <- function(model, newdata, batch_size, ...) {
    rep.int(0L, nrow(newdata))
  }
  testthat::local_mocked_bindings(
    .autok_default_fit = fake_fit,
    .autok_default_predict = fake_predict,
    .package = "CellCharterR"
  )
  x <- matrix(seq_len(24), ncol = 3L)
  invisible(cellcharter_autok_fit(
    x, c(2L, 2L), checkpoint_dir, max_runs = 1L,
    fit_args = list(batch_size = 8L, device = "cpu", num_threads = 1L),
    predict_args = list(batch_size = 8L, device = "cpu"), progress = FALSE
  ))

  state_path <- file.path(checkpoint_dir, "autok_state.rds")
  state <- readRDS(state_path)
  artifact_path <- state$fit_index$path[[1L]]
  artifact <- readRDS(artifact_path)
  artifact$execution_contract <- NULL
  artifact$effective_fit_batch_size <- 4L
  saveRDS(artifact, artifact_path)
  state$fit_index$bytes[[1L]] <- as.numeric(file.info(artifact_path)$size)
  state$fit_index$sha256[[1L]] <- .autok_file_sha256(artifact_path)
  state$fit_index$execution_contract_sha256 <- NULL
  saveRDS(state, state_path)

  expect_error(
    cellcharter_autok_fit(
      x, c(2L, 2L), checkpoint_dir, max_runs = 1L,
      fit_args = list(batch_size = 8L, device = "cpu", num_threads = 1L),
      predict_args = list(batch_size = 8L, device = "cpu"), progress = FALSE
    ),
    "cannot prove the requested numerical execution contract"
  )
})

test_that("Torch legacy state without an execution contract fails closed", {
  state <- list(
    backend = "torch",
    resume_contract = list(algorithm = "torch-full-covariance-em-v1"),
    fit_args = list(batch_size = 8L, device = "auto"),
    predict_args = list(batch_size = 8L, device = "auto")
  )
  expect_null(.autok_saved_execution_contract(state))
})

test_that("direct AutoK model loading verifies integrity and repairs legacy C++ class", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  legacy_model <- structure(
    list(runtime = list(backend = "compiled-cpp-double"), n_clusters = 2L),
    class = "torch_gmm"
  )
  saveRDS(list(model = legacy_model, backend = "cpp"), path)
  row <- data.frame(
    k = 2L, path = path, bytes = as.numeric(file.info(path)$size),
    sha256 = .autok_file_sha256(path), stringsAsFactors = FALSE
  )
  object <- structure(list(
    candidate_k = 2L, stability_mean = 1, best_k = 2L, best_index = row
  ), class = "cellcharter_autok")

  restored <- cellcharter_autok_model(object)
  expect_s3_class(restored, "cpp_gmm")
  expect_true(inherits(restored, "torch_gmm"))

  object$best_index$sha256 <- paste(rep.int("0", 64L), collapse = "")
  expect_error(cellcharter_autok_model(object), "SHA-256 integrity")
})

test_that("C++ fit checkpoint mode publishes each K before starting the next", {
  checkpoint_dir <- tempfile("autok-cpp-per-fit-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  calls <- 0L
  testthat::local_mocked_bindings(
    cpp_gmm_fit = function(x, n_clusters, seed, num_threads, ...) {
      calls <<- calls + 1L
      if (calls == 2L) stop("synthetic interruption after first K")
      structure(list(
        labels = rep.int(0L, nrow(x)), nll = as.numeric(n_clusters),
        nll_ = as.numeric(n_clusters), runtime = list(backend = "compiled-cpp-double")
      ), class = c("cpp_gmm", "torch_gmm"))
    },
    .package = "CellCharterR"
  )

  expect_error(
    cellcharter_autok_fit(
      matrix(seq_len(24), ncol = 3L), c(2L, 2L), checkpoint_dir,
      max_runs = 1L, backend = "cpp", checkpoint_mode = "fit",
      progress = FALSE
    ),
    "synthetic interruption"
  )
  first <- .autok_checkpoint_path(checkpoint_dir, 1L, 1L)
  second <- .autok_checkpoint_path(checkpoint_dir, 2L, 1L)
  expect_true(file.exists(first))
  expect_false(file.exists(second))
  expect_identical(readRDS(first)$k, 1L)
})

test_that("ordinary AutoK calls default to the supported compiled backend", {
  testthat::local_mocked_bindings(
    cpp_gmm_fit = function(x, n_clusters, seed, num_threads, ...) {
      structure(list(
        labels = (seq_len(nrow(x)) - 1L) %% n_clusters,
        selection_nll = as.numeric(n_clusters), nll = as.numeric(n_clusters),
        nll_history = as.numeric(n_clusters), converged = TRUE,
        iterations = 1L, covariance_regularization = 1e-6,
        regularization_attempts = 1e-6,
        runtime = list(device = "cpu", backend = "compiled-cpp-double")
      ), class = c("cpp_gmm", "torch_gmm"))
    },
    .package = "CellCharterR"
  )
  fit <- cellcharter_autok_fit(
    matrix(seq_len(30), ncol = 3L), c(2L, 2L),
    tempfile("autok-default-cpp-"), max_runs = 1L, progress = FALSE
  )
  expect_identical(fit$backend, "cpp")
  expect_error(
    cellcharter_autok_fit(
      matrix(seq_len(30), ncol = 3L), c(2L, 2L),
      tempfile("autok-cpp-hooks-"), max_runs = 1L,
      backend = "cpp", fit_fun = function(...) NULL, progress = FALSE
    ),
    "does not call R hooks"
  )
})

test_that("C++ AutoK ranks repetitions by the reference monitored NLL", {
  x <- matrix(seq_len(48), ncol = 4L)
  testthat::local_mocked_bindings(
    cpp_gmm_fit_many = function(x, n_clusters, seeds, ...) {
      lapply(seq_along(n_clusters), function(i) {
        k <- n_clusters[[i]]
        selection_nll <- as.numeric(seeds[[i]] - 9L)
        structure(list(
          labels = (seq_len(nrow(x)) - 1L) %% k,
          selection_nll = selection_nll,
          nll = 3 - selection_nll,
          nll_history = selection_nll,
          runtime = list(device = "cpu", backend = "compiled-cpp-double")
        ), class = c("cpp_gmm", "torch_gmm"))
      })
    },
    .package = "CellCharterR"
  )

  fit <- cellcharter_autok_fit(
    x, c(2L, 2L), tempfile("autok-selection-nll-"),
    max_runs = 2L, convergence_tol = 0, seed = 10L,
    backend = "cpp", workers = 1L, checkpoint_mode = "none",
    progress = FALSE
  )

  expect_true(all(fit$best_index$run == 1L))
  expect_true(all(fit$best_index$nll == 1))
  for (artifact in fit$embedded_best_artifacts) {
    expect_equal(artifact$nll, 1, tolerance = 0)
    expect_equal(artifact$final_nll, 2, tolerance = 0)
    expect_equal(artifact$model$selection_nll, 1, tolerance = 0)
    expect_equal(artifact$model$nll, 2, tolerance = 0)
  }
})

test_that("C++ numerical ABI rejects older checkpoint identities", {
  for (version in c("v1.0.1", "v1.1.0", "v1.2.0")) {
    previous <- list(
      algorithm = paste0(
        "cpp-full-covariance-em-", version, ":",
        paste(rep.int("a", 64L), collapse = "")
      ),
      representation_sha256 = "representation",
      input_data_id = "input",
      fit_args = list(),
      predict_args = list()
    )
    state <- list(resume_contract = previous, data_id = "legacy-data-id")
    expect_identical(
      .autok_cpp_observability_compatible_id(state, previous),
      character()
    )
  }
})

test_that("allocation-only FMI hardening preserves pre-hardening resume identity", {
  old_digest <- "8b479382a95c7a4ef4354fd590cbaece42c3923d6397b12fe7caeae65c55591c"
  hardened_digest <- "19b99ee2dfbf401d643225618a28a9683bae77c9a9c2a208a585afad0c172688"
  current <- list(
    fit_k = 1:4, candidate_k = 2:3, seed = 1L,
    n_obs = 10L, n_features = 3L
  )
  state <- c(current, list(
    data_id = "compatible-id", runs_completed = 1L,
    edge_stability = data.frame(), convergence_history = data.frame(),
    similarity_digest = old_digest
  ))
  expect_true(.autok_compatible_state(
    state, current, accepted_data_ids = "compatible-id",
    similarity_digest = hardened_digest, default_similarity = TRUE
  ))
  expect_false(.autok_compatible_state(
    state, current, accepted_data_ids = "compatible-id",
    similarity_digest = hardened_digest, default_similarity = FALSE
  ))
})

test_that("AutoK disables nested direct-fit OOM retries", {
  observed <- NULL
  testthat::local_mocked_bindings(
    torch_gmm_fit = function(...) {
      observed <<- list(...)
      "model"
    },
    .package = "CellCharterR"
  )
  expect_identical(
    .autok_default_fit(matrix(1, 1L, 1L), 1L, 2L, batch_size = 64L),
    "model"
  )
  expect_identical(observed$.oom_retries, 0L)
  expect_identical(observed$batch_size, 64L)
})

test_that("AutoK rejects a checkpoint from another representation", {
  checkpoint_dir <- tempfile("autok-mismatch-")
  x <- matrix(seq_len(24), ncol = 3)
  fit_hook <- function(x, n_clusters, seed, ...) list(k = n_clusters, nll = n_clusters)
  predict_hook <- function(model, newdata, ...) rep.int(0L, nrow(newdata))

  invisible(cellcharter_autok_fit(
    x, c(2L, 3L), checkpoint_dir, max_runs = 1L,
    fit_fun = fit_hook, predict_fun = predict_hook,
    nll_fun = function(model) model$nll, data_id = "first"
  ))
  expect_error(
    cellcharter_autok_fit(
      x, c(2L, 3L), checkpoint_dir, max_runs = 1L,
      fit_fun = fit_hook, predict_fun = predict_hook,
      nll_fun = function(model) model$nll, data_id = "second"
    ),
    "does not match current `data_id`"
  )
})

test_that("AutoK binds resume identity to same-shaped representation contents", {
  checkpoint_dir <- tempfile("autok-content-identity-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  x <- matrix(seq_len(24), ncol = 3)
  fit_hook <- function(x, n_clusters, seed, ...) list(k = n_clusters, nll = n_clusters)
  predict_hook <- function(model, newdata, ...) rep.int(0L, nrow(newdata))

  first <- cellcharter_autok_fit(
    x, c(2L, 3L), checkpoint_dir, max_runs = 1L,
    fit_fun = fit_hook, predict_fun = predict_hook,
    nll_fun = function(model) model$nll, data_id = "same-provenance",
    progress = FALSE
  )
  changed <- x
  changed[[1L]] <- changed[[1L]] + 0.5
  expect_false(identical(
    first$representation_sha256,
    .autok_representation_sha256(changed)
  ))
  expect_error(
    cellcharter_autok_fit(
      changed, c(2L, 3L), checkpoint_dir, max_runs = 1L,
      fit_fun = fit_hook, predict_fun = predict_hook,
      nll_fun = function(model) model$nll, data_id = "same-provenance",
      progress = FALSE
    ),
    "representation contents"
  )
})

test_that("representation fingerprints are storage-backend independent", {
  x <- matrix(seq_len(30), ncol = 3,
              dimnames = list(paste0("spot", 1:10), paste0("f", 1:3)))
  delayed <- DelayedArray::DelayedArray(x)
  expect_identical(
    .autok_representation_sha256(x, target_block_values = 7L),
    .autok_representation_sha256(delayed, target_block_values = 7L)
  )
})

test_that("committed fit artifacts have cryptographic integrity guards", {
  checkpoint_dir <- tempfile("autok-artifact-integrity-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  x <- matrix(seq_len(24), ncol = 3)
  fit_hook <- function(x, n_clusters, seed, ...) list(k = n_clusters, nll = n_clusters)
  predict_hook <- function(model, newdata, ...) rep.int(0L, nrow(newdata))
  first <- cellcharter_autok_fit(
    x, c(2L, 3L), checkpoint_dir, max_runs = 1L,
    fit_fun = fit_hook, predict_fun = predict_hook,
    nll_fun = function(model) model$nll, progress = FALSE
  )
  expect_true(all(grepl("^[0-9a-f]{64}$", first$fit_index$sha256)))

  path <- first$fit_index$path[[1L]]
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  position <- max(1L, length(bytes) - 10L)
  bytes[[position]] <- as.raw(bitwXor(as.integer(bytes[[position]]), 1L))
  writeBin(bytes, path)
  expect_error(
    cellcharter_autok_fit(
      x, c(2L, 3L), checkpoint_dir, max_runs = 1L,
      fit_fun = fit_hook, predict_fun = predict_hook,
      nll_fun = function(model) model$nll, progress = FALSE
    ),
    "SHA-256 integrity"
  )
})

test_that("AutoK stability diagnostics preserve scores and expose peaks", {
  stability <- rbind(
    `2` = c(0.50, 0.70, 0.60),
    `3` = c(0.75, 0.85, 0.80),
    `4` = c(0.55, 0.65, 0.60),
    `5` = c(0.68, 0.72, 0.70)
  )
  fit <- structure(list(
    candidate_k = 2:5,
    stability = stability,
    stability_mean = rowMeans(stability),
    best_k = 3L
  ), class = "cellcharter_autok")

  diagnostic <- cellcharter_autok_stability_data(
    fit, k_range = 1:5, highlight_k = c(4L, 5L)
  )
  expect_equal(diagnostic$summary$k, 1:5)
  expect_false(diagnostic$summary$estimable[[1L]])
  expect_true(is.na(diagnostic$summary$mean_fowlkes_mallows[[1L]]))
  expect_equal(nrow(diagnostic$scores), length(stability))
  expect_equal(
    diagnostic$scores$fowlkes_mallows,
    as.numeric(t(stability))
  )
  expect_equal(diagnostic$summary$k[diagnostic$summary$local_peak], 3L)
  expect_equal(diagnostic$summary$k[diagnostic$summary$selected], 3L)
  expect_equal(diagnostic$summary$k[diagnostic$summary$highlighted], c(4L, 5L))
})

test_that("AutoK stability plot carries its exact diagnostic data", {
  skip_if_not_installed("ggplot2")
  stability <- rbind(`2` = c(0.5, 0.6), `3` = c(0.7, 0.8), `4` = c(0.6, 0.7))
  fit <- structure(list(
    candidate_k = 2:4, stability = stability,
    stability_mean = rowMeans(stability), best_k = 3L
  ), class = "cellcharter_autok")
  plot <- cellcharter_plot_autok_stability(fit, 1:4, highlight_k = c(3L, 4L))
  expect_s3_class(plot, "ggplot")
  expect_equal(
    attr(plot, "cellcharter_stability_data")$summary$k,
    1:4
  )
})
test_that("AutoK safely accepts producer-stage identity migration", {
  representation_sha256 <- paste(rep("a", 64L), collapse = "")
  current <- .autok_resume_contract(
    input_data_id = representation_sha256,
    fit_args = list(max_iter = 10L), predict_args = list(),
    default_torch = FALSE, implementation_digest = "implementation",
    algorithm = "cpp-full-covariance-em-v1.5.0",
    representation_sha256 = representation_sha256
  )
  previous <- current
  previous$input_data_id <- "historical-producer-stage-id"
  state <- list(
    resume_contract = previous,
    data_id = digest::digest(
      .cc_contract_canonical(previous), algo = "sha256", serialize = TRUE
    )
  )

  expect_identical(
    .autok_input_identity_compatible_id(state, current),
    state$data_id
  )
  changed <- current
  changed$representation_sha256 <- paste(rep("b", 64L), collapse = "")
  expect_length(.autok_input_identity_compatible_id(state, changed), 0L)
})
