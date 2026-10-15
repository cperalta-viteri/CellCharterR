autok_counting_hooks <- function(counter, delay = 0) {
  list(
    fit = function(x, n_clusters, seed, ...) {
      if (!is.null(counter$fail_after) && counter$fit >= counter$fail_after) {
        stop("synthetic interruption")
      }
      counter$fit <- counter$fit + 1L
      list(k = n_clusters, seed = seed, nll = n_clusters + seed / 1000)
    },
    predict = function(model, newdata, ...) {
      counter$predict <- counter$predict + 1L
      ((seq_len(nrow(newdata)) - 1L + model$seed) %% model$k) + 1L
    },
    similarity = function(labels_a, labels_b) {
      counter$similarity <- counter$similarity + 1L
      if (delay > 0) Sys.sleep(delay)
      cellcharter_fowlkes_mallows(labels_a, labels_b)
    }
  )
}

autok_test_call <- function(checkpoint_dir, counter, max_runs = 3L,
                            convergence_tol = 0, hooks = NULL) {
  if (is.null(hooks)) hooks <- autok_counting_hooks(counter)
  cellcharter_autok_fit(
    matrix(seq_len(72), ncol = 3),
    n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
    max_runs = max_runs, convergence_tol = convergence_tol, seed = 19L,
    fit_fun = hooks$fit, predict_fun = hooks$predict,
    nll_fun = function(model) model$nll,
    similarity_fun = hooks$similarity, data_id = "safe-resume-fixture"
  )
}

test_that("a completed AutoK campaign resumes without computational hooks", {
  checkpoint_dir <- tempfile("autok-complete-resume-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$fit <- counter$predict <- counter$similarity <- 0L

  original <- autok_test_call(checkpoint_dir, counter)
  counter$fit <- counter$predict <- counter$similarity <- 0L
  resumed <- autok_test_call(checkpoint_dir, counter)

  expect_identical(
    c(fit = counter$fit, predict = counter$predict, similarity = counter$similarity),
    c(fit = 0L, predict = 0L, similarity = 0L)
  )
  expect_equal(resumed$edge_stability, original$edge_stability)
  expect_equal(resumed$convergence_history, original$convergence_history)
  expect_identical(resumed$best_k, original$best_k)
  expect_true(all(resumed$fit_index$resumed))
})

test_that("legacy state metadata is migrated once even when no run remains", {
  checkpoint_dir <- tempfile("autok-legacy-state-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$fit <- counter$predict <- counter$similarity <- 0L

  hooks <- autok_counting_hooks(counter)
  legacy_call <- function() {
    cellcharter_autok_fit(
      matrix(seq_len(72), ncol = 3),
      n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
      max_runs = 3L, convergence_tol = 0, seed = 19L,
      fit_fun = hooks$fit, predict_fun = hooks$predict,
      nll_fun = function(model) model$nll,
      data_id = "safe-resume-legacy-fixture"
    )
  }
  original <- legacy_call()
  state_path <- file.path(checkpoint_dir, "autok_state.rds")
  legacy <- readRDS(state_path)
  legacy$fit_index <- NULL
  legacy$similarity_digest <- NULL
  saveRDS(legacy, state_path)

  counter$fit <- counter$predict <- counter$similarity <- 0L
  migrated <- legacy_call()
  upgraded <- readRDS(state_path)
  expect_identical(
    c(fit = counter$fit, predict = counter$predict, similarity = counter$similarity),
    c(fit = 0L, predict = 0L, similarity = 0L)
  )
  expect_equal(migrated$edge_stability, original$edge_stability)
  expect_equal(nrow(upgraded$fit_index), nrow(original$fit_index))
  expect_match(upgraded$similarity_digest, "^[0-9a-f]{64}$")

  state_time <- file.info(state_path)$mtime
  counter$fit <- counter$predict <- counter$similarity <- 0L
  second <- legacy_call()
  expect_equal(second$edge_stability, original$edge_stability)
  expect_identical(counter$similarity, 0L)
  expect_gte(file.info(state_path)$mtime, state_time)
})

test_that("a validated legacy input identity transitions across repeated resumes", {
  skip_if_not_installed("torch")
  checkpoint_dir <- tempfile("autok-input-transition-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  x <- matrix(c(-2, -2, -1.8, -2.1, 2, 2, 2.1, 1.8),
              ncol = 2L, byrow = TRUE)
  call <- function(input_id) {
    cellcharter_autok_fit(
      x, c(2L, 2L), checkpoint_dir, max_runs = 1L, seed = 37L,
      fit_args = list(batch_size = 4L, max_iter = 2L, device = "cpu",
                      num_threads = 1L, initialization_method = "torch"),
      predict_args = list(batch_size = 4L, device = "cpu"),
      data_id = input_id
    )
  }
  original <- call("historical-workflow-contract")
  state_path <- file.path(checkpoint_dir, "autok_state.rds")
  state <- readRDS(state_path)
  state$compatible_data_ids <- unique(c(state$compatible_data_ids, state$data_id))
  state$legacy_input_data_ids <- state$input_data_id
  state$input_data_id <- "stable-representation-id"
  saveRDS(state, state_path)

  migrated <- call("stable-representation-id")
  second <- call("stable-representation-id")
  expect_true(all(migrated$fit_index$resumed))
  expect_true(all(second$fit_index$resumed))
  expect_identical(second$input_data_id, "stable-representation-id")
  expect_true(original$data_id %in% second$compatible_data_ids)
})

test_that("an incomplete final round is quarantined and restarted atomically", {
  checkpoint_dir <- tempfile("autok-partial-resume-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$fit <- counter$predict <- counter$similarity <- 0L
  hooks <- autok_counting_hooks(counter)
  counter$fail_after <- 6L # four K in round 1, then two artifacts in round 2

  expect_error(
    autok_test_call(checkpoint_dir, counter, hooks = hooks),
    "synthetic interruption"
  )
  round_one <- vapply(1:4, function(k) {
    digest::digest(file = .autok_checkpoint_path(checkpoint_dir, k, 1L),
                   algo = "sha256")
  }, character(1))
  expect_true(file.exists(.autok_checkpoint_path(checkpoint_dir, 1L, 2L)))
  expect_true(file.exists(.autok_checkpoint_path(checkpoint_dir, 2L, 2L)))

  counter$fit <- counter$predict <- counter$similarity <- 0L
  counter$fail_after <- NULL
  completed <- autok_test_call(checkpoint_dir, counter)

  # Round 2 must restart all four K; round 3 is then fitted normally. Round 1
  # remains reusable and must not be rewritten.
  expect_identical(counter$fit, 8L)
  expect_identical(counter$predict, 8L)
  expect_identical(counter$similarity, 9L)
  expect_true(all(completed$fit_index$resumed[completed$fit_index$run == 1L]))
  expect_false(any(completed$fit_index$resumed[completed$fit_index$run == 2L]))
  expect_identical(vapply(1:4, function(k) {
    digest::digest(file = .autok_checkpoint_path(checkpoint_dir, k, 1L),
                   algo = "sha256")
  }, character(1)), round_one)

  quarantined <- list.files(checkpoint_dir, recursive = TRUE, full.names = TRUE)
  quarantined <- quarantined[
    grepl("quarant", quarantined, ignore.case = TRUE) &
      grepl("run_002[.]rds$", quarantined)
  ]
  expect_gte(length(quarantined), 2L)
})

test_that("a changed convergence policy is evaluated from completed rounds", {
  checkpoint_dir <- tempfile("autok-policy-resume-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$fit <- counter$predict <- counter$similarity <- 0L

  strict <- autok_test_call(
    checkpoint_dir, counter, max_runs = 5L, convergence_tol = 0
  )
  expect_identical(strict$runs_completed, 5L)
  counter$fit <- counter$predict <- counter$similarity <- 0L
  relaxed <- autok_test_call(
    checkpoint_dir, counter, max_runs = 5L, convergence_tol = 1
  )

  expect_identical(counter$fit, 0L)
  expect_identical(counter$predict, 0L)
  expect_true(all(relaxed$fit_index$resumed))
  expect_identical(relaxed$runs_completed, 3L)
  expect_true(isTRUE(relaxed$converged))
  expect_identical(relaxed$convergence_tol, 1)
})
