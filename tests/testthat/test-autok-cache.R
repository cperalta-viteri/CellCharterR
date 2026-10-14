test_that("AutoK caches current-round labels without changing stability", {
  checkpoint_dir <- tempfile("autok-label-cache-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  x <- matrix(seq_len(48), ncol = 3)
  fit_hook <- function(x, n_clusters, seed, ...) {
    list(k = n_clusters, seed = seed, nll = n_clusters + seed / 1000)
  }
  predict_hook <- function(model, newdata, ...) {
    ((seq_len(nrow(newdata)) - 1L + model$seed) %% model$k)
  }

  reads <- 0L
  original_read_labels <- CellCharterR:::.autok_read_labels
  testthat::local_mocked_bindings(
    .autok_read_labels = function(path, n_obs) {
      reads <<- reads + 1L
      original_read_labels(path, n_obs)
    },
    .package = "CellCharterR"
  )
  fit <- cellcharter_autok_fit(
    x, n_clusters = c(2L, 3L), checkpoint_dir = checkpoint_dir,
    max_runs = 3L, convergence_tol = 0, seed = 11L,
    fit_fun = fit_hook, predict_fun = predict_hook,
    nll_fun = function(model) model$nll, data_id = "cache-fixture"
  )

  # Four fitted K values give three neighboring edges. Current-round labels
  # are retained from the fit loop; only 3 * (r - 1) prior labels are read.
  expect_identical(reads, 9L)

  expected <- vector("list", 9L)
  index <- 0L
  for (run in 2:3) {
    for (previous_run in seq_len(run - 1L)) {
      for (edge in 1:3) {
        lower <- readRDS(.autok_checkpoint_path(checkpoint_dir, edge, run))$labels
        upper <- readRDS(.autok_checkpoint_path(checkpoint_dir, edge + 1L, previous_run))$labels
        index <- index + 1L
        expected[[index]] <- data.frame(
          run = run, previous_run = previous_run, lower_k = edge,
          upper_k = edge + 1L,
          score = cellcharter_fowlkes_mallows(lower, upper)
        )
      }
    }
  }
  expect_equal(fit$edge_stability, do.call(rbind, expected))
})

test_that("cached stability label reader rejects malformed checkpoints", {
  checkpoint_dir <- tempfile("autok-label-invalid-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  path <- .autok_checkpoint_path(checkpoint_dir, 2L, 1L)
  dir.create(dirname(path), recursive = TRUE)
  saveRDS(list(labels = c(1L, NA_integer_)), path)

  expect_error(
    .autok_read_labels(path, n_obs = 2L),
    "invalid labels"
  )
})
