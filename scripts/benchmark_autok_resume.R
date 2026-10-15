#!/usr/bin/env Rscript

# Reproducible, data-free benchmark for AutoK resume overhead. It compares a
# legacy-style resume that replays all similarity calculations with the public
# optimized resume path. All artifacts live under a fresh temporary directory.

if (file.exists("DESCRIPTION")) {
  # Source only the implementation needed by this isolated benchmark. This
  # avoids depending on a preinstalled copy of the package under development.
  source(file.path("R", "provenance.R"), local = globalenv())
  source(file.path("R", "torch_gmm.R"), local = globalenv())
  source(file.path("R", "autok.R"), local = globalenv())
} else if (requireNamespace("CellCharterR", quietly = TRUE)) {
  library(CellCharterR)
} else {
  stop("Run from the source repository with `devtools` installed, or install CellCharterR first.")
}

counter <- new.env(parent = emptyenv())
counter$fit <- counter$predict <- counter$similarity <- 0L
fit_hook <- function(x, n_clusters, seed, ...) {
  counter$fit <- counter$fit + 1L
  list(k = n_clusters, seed = seed, nll = n_clusters + seed / 1000)
}
predict_hook <- function(model, newdata, ...) {
  counter$predict <- counter$predict + 1L
  ((seq_len(nrow(newdata)) - 1L + model$seed) %% model$k) + 1L
}
similarity_hook <- function(labels_a, labels_b) {
  counter$similarity <- counter$similarity + 1L
  Sys.sleep(0.01) # make replay cost measurable without a large dataset
  cellcharter_fowlkes_mallows(labels_a, labels_b)
}

root <- tempfile("CellCharterR-autok-resume-benchmark-")
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
checkpoint_dir <- file.path(root, "isolated-checkpoints")
checkpoint_path <- function(k, run) {
  file.path(checkpoint_dir, "fits", sprintf("k_%03d", k),
            sprintf("run_%03d.rds", run))
}
x <- matrix(seq_len(1200), ncol = 6)
args <- list(
  x = x, n_clusters = c(2L, 5L), checkpoint_dir = checkpoint_dir,
  max_runs = 5L, convergence_tol = 0, seed = 41L,
  fit_fun = fit_hook, predict_fun = predict_hook,
  nll_fun = function(model) model$nll, similarity_fun = similarity_hook,
  data_id = "benchmark-synthetic-v1"
)
invisible(do.call(cellcharter_autok_fit, args))

measure <- function(expression) {
  counter$fit <- counter$predict <- counter$similarity <- 0L
  elapsed <- unname(system.time(force(expression))[["elapsed"]])
  data.frame(elapsed_seconds = elapsed, fit_calls = counter$fit,
             predict_calls = counter$predict,
             similarity_calls = counter$similarity)
}

# This explicitly models the old behavior: reopen completed label artifacts and
# replay every neighboring-K comparison on every resume.
legacy_replay <- function() {
  fit_k <- 1:6
  for (run in 2:5) {
    for (previous_run in seq_len(run - 1L)) {
      for (edge in seq_len(length(fit_k) - 1L)) {
        lower <- readRDS(checkpoint_path(fit_k[[edge]], run))$labels
        upper <- readRDS(checkpoint_path(
          fit_k[[edge + 1L]], previous_run
        ))$labels
        similarity_hook(lower, upper)
      }
    }
  }
}

old <- measure(legacy_replay())
new <- measure(do.call(cellcharter_autok_fit, args))
result <- rbind(
  cbind(implementation = "old_similarity_replay", old),
  cbind(implementation = "optimized_completed_resume", new)
)
result$speedup <- result$elapsed_seconds[[1L]] /
  pmax(result$elapsed_seconds, .Machine$double.eps)
print(result, row.names = FALSE)
stopifnot(new$fit_calls == 0L, new$predict_calls == 0L,
          new$similarity_calls == 0L)
