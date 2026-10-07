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

  no_refit <- function(...) stop("a completed checkpoint was unexpectedly refitted")
  resumed <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 4L), checkpoint_dir = checkpoint_dir,
    max_runs = 10L, convergence_tol = 0.01, seed = 41L,
    fit_fun = no_refit, predict_fun = no_refit,
    nll_fun = function(model) model$nll, data_id = "fixture-v1", resume = TRUE
  )
  expect_equal(resumed$runs_completed, 3L)
  expect_true(all(resumed$fit_index$resumed))
  expect_equal(resumed$stability, fit$stability)
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
