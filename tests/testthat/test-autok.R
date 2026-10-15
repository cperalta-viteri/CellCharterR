test_that("Fowlkes-Mallows is label-invariant and handles singleton partitions", {
  a <- c(0, 0, 1, 1, 2, 2)
  expect_equal(cellcharter_fowlkes_mallows(a, c("x", "x", "z", "z", "q", "q")), 1)
  expect_equal(cellcharter_fowlkes_mallows(a, seq_along(a)), 0)
  expect_error(cellcharter_fowlkes_mallows(a, a[-1]), "same length")
  expect_error(cellcharter_fowlkes_mallows(a, replace(a, 1, NA)), "missing")
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
