#!/usr/bin/env Rscript

# Fresh compiled-CPU AutoK run from the already validated, disk-backed mouse
# spleen aggregation. This intentionally writes no intermediate fit or round
# checkpoints; only completed final deliverables are saved.
suppressPackageStartupMessages({
  library(CellCharterR)
  library(SingleCellExperiment)
})

input_path <- Sys.getenv(
  "CODEX_AGGREGATED_SCE",
  "results/CellCharterR_mouse_spleen_v2/aggregation_pinned_global/codex_sce_aggregated.rds"
)
output_dir <- Sys.getenv(
  "CODEX_CPP_RESULTS_DIR", "results/CellCharterR_mouse_spleen_cpp_v1"
)
workers <- as.integer(Sys.getenv("CODEX_AUTOK_WORKERS", "8"))
threads_per_fit <- as.integer(Sys.getenv("CODEX_AUTOK_THREADS_PER_FIT", "4"))
seed <- as.integer(Sys.getenv("CODEX_SEED", "0"))
max_runs <- as.integer(Sys.getenv("CODEX_MAX_RUNS", "10"))

if (!file.exists(input_path)) stop("Missing aggregated mouse-spleen SCE: ", input_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sce <- readRDS(input_path)
x <- SingleCellExperiment::reducedDim(sce, "X_cellcharter")
input_sha256 <- CellCharterR::cellcharter_file_sha256(input_path)

started <- Sys.time()
autok <- CellCharterR::cellcharter_autok_fit(
  x,
  n_clusters = c(2L, 20L),
  checkpoint_dir = file.path(output_dir, "unused-no-checkpoints"),
  max_runs = max_runs,
  convergence_tol = 0.01,
  seed = seed,
  backend = "cpp",
  workers = workers,
  threads_per_fit = threads_per_fit,
  checkpoint_mode = "none",
  resume = FALSE,
  fit_args = list(
    max_iter = 100L,
    tolerance = 0.001,
    covariance_regularization = 1e-6
  ),
  progress = TRUE,
  data_id = paste0("sha256:", input_sha256)
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

selected_model <- CellCharterR::cellcharter_autok_model(autok)
selected_labels <- selected_model$labels
if (length(selected_labels) != ncol(sce)) {
  stop("Selected model did not return one label per mouse-spleen spot.")
}

# Fit the tutorial-matched K=11 sensitivity separately with its published seed.
k11_started <- Sys.time()
k11_model <- CellCharterR::cpp_gmm_fit(
  as.matrix(x), n_clusters = 11L, max_iter = 100L, tolerance = 0.001,
  covariance_regularization = 1e-6, seed = 12345L,
  num_threads = threads_per_fit
)
k11_elapsed <- as.numeric(difftime(Sys.time(), k11_started, units = "secs"))

provenance <- list(
  backend = "compiled-cpp-double",
  backend_algorithm = autok$resume_contract$algorithm,
  reference_version = autok$reference_version,
  reference_commit = autok$reference_commit,
  initializer = selected_model$initialization$method,
  initializer_limitation = paste(
    "Deterministic C++ K-means++/Lloyd; not sklearn RNG-identical.",
    "The EM and scoring equations are the pinned two-pass full-covariance target."
  ),
  input_path = normalizePath(input_path),
  input_sha256 = input_sha256,
  dimensions = dim(x),
  workers = workers,
  threads_per_fit = threads_per_fit,
  seed = seed,
  max_runs = max_runs,
  runs_completed = autok$runs_completed,
  elapsed_seconds = elapsed,
  selected_k = autok$best_k,
  selected_nll = selected_model$nll,
  k11_seed = 12345L,
  k11_nll = k11_model$nll,
  k11_elapsed_seconds = k11_elapsed,
  finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  session_info = utils::sessionInfo()
)

# These are final analysis deliverables, not resumable intermediate checkpoints.
saveRDS(autok, file.path(output_dir, "autok_result.rds"), compress = FALSE)
saveRDS(selected_model, file.path(output_dir, "selected_model.rds"), compress = FALSE)
saveRDS(list(k = autok$best_k, labels = selected_labels, seed = seed),
        file.path(output_dir, "selected_labels.rds"), compress = FALSE)
saveRDS(k11_model, file.path(output_dir, "tutorial_k11_model.rds"), compress = FALSE)
saveRDS(list(k = 11L, labels = k11_model$labels, seed = 12345L),
        file.path(output_dir, "tutorial_k11_labels.rds"), compress = FALSE)
saveRDS(provenance, file.path(output_dir, "provenance.rds"), compress = FALSE)

cat(sprintf(
  "Compiled AutoK complete: selected K=%d after %d runs in %.1f minutes; K=11 sensitivity %.1f minutes.\n",
  autok$best_k, autok$runs_completed, elapsed / 60, k11_elapsed / 60
))
