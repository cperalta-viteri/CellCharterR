#!/usr/bin/env Rscript

# Reproducible compiled-CPU mouse-spleen AutoK benchmark.
#
# Run this script from the repository root. Set BLAS/OpenMP controls before R
# starts; changing them from inside a running R process can be too late after a
# BLAS library has initialized:
#
#   OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
#   VECLIB_MAXIMUM_THREADS=1 GOTO_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#   CODEX_AUTOK_WORKERS=8 CODEX_AUTOK_THREADS_PER_FIT=4 \
#   CODEX_CPP_RESULTS_DIR=results/CellCharterR_mouse_spleen_cpp_reproduction \
#   Rscript scripts/run_cpp_mouse_spleen_autok.R
#
# This runner always uses the compiled C++ CPU backend. It does not initialize
# torch or CUDA. AutoK round checkpoints make an interrupted run resumable;
# completed deliverables are published atomically under results/.

.env_positive_integer <- function(name, default) {
  raw <- Sys.getenv(name, as.character(default))
  if (!grepl("^[0-9]+$", raw)) {
    stop(name, " must be a positive integer; got `", raw, "`.", call. = FALSE)
  }
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer; got `", raw, "`.", call. = FALSE)
  }
  value
}

.atomic_save_rds <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), dirname(path), ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  base::saveRDS(object, temporary, compress = compress)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically publish `", path, "`.", call. = FALSE)
  }
  invisible(path)
}

.peak_r_memory_mb <- function() {
  usage <- gc()
  as.numeric(
    usage["Ncells", "max used"] * 56 + usage["Vcells", "max used"] * 8
  ) / 1024^2
}

.peak_process_memory_mb <- function() {
  # Linux VmHWM covers the complete process, including native allocations and
  # all C++ worker threads. It is unavailable on some operating systems.
  if (!file.exists("/proc/self/status")) return(NA_real_)
  status <- readLines("/proc/self/status", warn = FALSE)
  line <- status[grepl("^VmHWM:", status)]
  if (length(line) != 1L) return(NA_real_)
  value_kb <- suppressWarnings(as.numeric(sub(
    "^VmHWM:[[:space:]]*([0-9]+)[[:space:]]+kB.*$", "\\1", line
  )))
  if (length(value_kb) != 1L || !is.finite(value_kb)) NA_real_ else value_kb / 1024
}

.git_output <- function(args) {
  tryCatch(
    system2("git", args, stdout = TRUE, stderr = FALSE),
    warning = function(e) character(), error = function(e) character()
  )
}

.sha256_if_file <- function(path) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    return(NA_character_)
  }
  CellCharterR::cellcharter_file_sha256(path)
}

# These variables control BLAS/OpenMP implementations rather than the native
# C++ AutoK topology below. Keeping them at one avoids hidden oversubscription.
.thread_environment_names <- c(
  "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "GOTO_NUM_THREADS", "OMP_NUM_THREADS"
)
.thread_environment <- Sys.getenv(.thread_environment_names, unset = NA_character_)
names(.thread_environment) <- .thread_environment_names
.unsafe_thread_environment <- is.na(.thread_environment) | .thread_environment != "1"
if (any(.unsafe_thread_environment)) {
  warning(
    paste0(
      "For a reproducible CPU benchmark, set these variables to 1 before R starts: ",
      paste(.thread_environment_names[.unsafe_thread_environment], collapse = ", "),
      ". See the invocation at the top of this script."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(CellCharterR)
  library(SingleCellExperiment)
})

runner_started <- Sys.time()
process_peak_at_start_mb <- .peak_process_memory_mb()

runner_path <- "scripts/run_cpp_mouse_spleen_autok.R"
if (!file.exists(runner_path) || !file.exists("DESCRIPTION")) {
  stop("Run this script from the CellCharterR repository root.", call. = FALSE)
}

repository_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
dir.create(file.path(repository_root, "results"), recursive = TRUE, showWarnings = FALSE)
results_root <- normalizePath("results", winslash = "/", mustWork = TRUE)
input_path <- Sys.getenv(
  "CODEX_AGGREGATED_SCE",
  "results/CellCharterR_mouse_spleen_v2/aggregation_pinned_global/codex_sce_aggregated.rds"
)
output_dir <- Sys.getenv(
  "CODEX_CPP_RESULTS_DIR", "results/CellCharterR_mouse_spleen_cpp_v1"
)
manifest_path <- Sys.getenv(
  "CODEX_AGGREGATION_MANIFEST",
  file.path(dirname(input_path), "aggregation_manifest.rds")
)

if (!file.exists(input_path)) {
  stop("Missing aggregated mouse-spleen SCE: ", input_path, call. = FALSE)
}
input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)
if (!startsWith(paste0(input_path, "/"), paste0(repository_root, "/"))) {
  stop("CODEX_AGGREGATED_SCE must resolve inside this repository.", call. = FALSE)
}
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
if (identical(output_dir, results_root) ||
    !startsWith(paste0(output_dir, "/"), paste0(results_root, "/"))) {
  stop("CODEX_CPP_RESULTS_DIR must be a subdirectory of results/.", call. = FALSE)
}
if (file.exists(file.path(output_dir, "provenance.rds"))) {
  stop(
    "A completed benchmark already exists in `", output_dir,
    "`. Choose a new CODEX_CPP_RESULTS_DIR to avoid overwriting it.",
    call. = FALSE
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Scientific settings are deliberately fixed for the primary benchmark.
candidate_k <- 2:20
max_runs <- 10L
convergence_tol <- 0.01
seed <- 0L
backend <- "cpp"
workers <- .env_positive_integer("CODEX_AUTOK_WORKERS", 8L)
threads_per_fit <- .env_positive_integer("CODEX_AUTOK_THREADS_PER_FIT", 4L)
if (workers > 32L || threads_per_fit > 32L ||
    as.double(workers) * as.double(threads_per_fit) > 32) {
  stop(
    "The compiled backend supports at most 32 requested native fit threads; ",
    "require CODEX_AUTOK_WORKERS * CODEX_AUTOK_THREADS_PER_FIT <= 32.",
    call. = FALSE
  )
}
logical_cores <- parallel::detectCores(logical = TRUE)
physical_cores <- parallel::detectCores(logical = FALSE)
if (is.finite(logical_cores) && workers * threads_per_fit > logical_cores) {
  warning(
    "Requested native fit topology exceeds the detected logical CPU count.",
    call. = FALSE
  )
}

input_sha256 <- CellCharterR::cellcharter_file_sha256(input_path)
aggregation_manifest <- NULL
aggregation_manifest_sha256 <- NA_character_
if (file.exists(manifest_path)) {
  manifest_path <- normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  if (!startsWith(paste0(manifest_path, "/"), paste0(repository_root, "/"))) {
    stop("CODEX_AGGREGATION_MANIFEST must resolve inside this repository.", call. = FALSE)
  }
  aggregation_manifest_sha256 <- CellCharterR::cellcharter_file_sha256(manifest_path)
  aggregation_manifest <- readRDS(manifest_path)
} else {
  warning(
    "No aggregation manifest was found; the input RDS SHA-256 will remain the provenance anchor.",
    call. = FALSE
  )
}

sce <- readRDS(input_path)
if (!inherits(sce, "SingleCellExperiment")) {
  stop("The input must be a SingleCellExperiment.", call. = FALSE)
}
if (!"X_cellcharter" %in% SingleCellExperiment::reducedDimNames(sce)) {
  stop("The input lacks reducedDim(sce, 'X_cellcharter').", call. = FALSE)
}
x <- SingleCellExperiment::reducedDim(sce, "X_cellcharter")
if (length(dim(x)) != 2L || nrow(x) != ncol(sce) || nrow(x) < max(candidate_k) + 1L) {
  stop("The X_cellcharter representation has incompatible dimensions.", call. = FALSE)
}
cell_ids <- colnames(sce)
if (is.null(cell_ids) || anyNA(cell_ids) || any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) {
  stop("The input SCE must have complete, unique cell identifiers.", call. = FALSE)
}
cell_order_sha256 <- digest::digest(as.character(cell_ids), algo = "sha256")
if (!is.null(aggregation_manifest)) {
  if (!is.null(aggregation_manifest$dimensions) &&
      !identical(as.integer(aggregation_manifest$dimensions), as.integer(dim(x)))) {
    stop("Aggregation manifest dimensions do not match X_cellcharter.", call. = FALSE)
  }
  if (!is.null(aggregation_manifest$cell_order_sha256) &&
      !identical(aggregation_manifest$cell_order_sha256, cell_order_sha256)) {
    stop("Aggregation manifest cell order does not match the input SCE.", call. = FALSE)
  }
}

native_dll_path <- tryCatch(
  getLoadedDLLs()[["CellCharterR"]][["path"]],
  error = function(e) NA_character_
)
source_sha256 <- stats::setNames(
  vapply(
    c(runner_path, "R/autok.R", "R/cpp_gmm.R", "src/cpp_gmm.cpp"),
    .sha256_if_file,
    character(1)
  ),
  c("runner", "autok_r", "cpp_gmm_r", "cpp_gmm_cpp")
)

message(sprintf(
  paste0(
    "Starting CPU-only AutoK: selectable K=2:20, max_runs=10, ",
    "convergence_tol=0.01, seed=0, workers=%d, threads_per_fit=%d."
  ),
  workers, threads_per_fit
))

autok_started <- Sys.time()
autok <- CellCharterR::cellcharter_autok_fit(
  x,
  n_clusters = c(2L, 20L),
  checkpoint_dir = file.path(output_dir, "autok_checkpoints"),
  max_runs = max_runs,
  convergence_tol = convergence_tol,
  seed = seed,
  backend = backend,
  workers = workers,
  threads_per_fit = threads_per_fit,
  checkpoint_mode = "round",
  resume = TRUE,
  fit_args = list(
    max_iter = 100L,
    tolerance = 0.001,
    covariance_regularization = 1e-6
  ),
  progress = TRUE,
  data_id = paste0("sha256:", input_sha256)
)
autok_finished <- Sys.time()
autok_elapsed <- as.numeric(difftime(autok_finished, autok_started, units = "secs"))

if (!identical(autok$backend, backend) ||
    !identical(as.integer(autok$requested_k), c(2L, 20L)) ||
    !identical(as.integer(autok$max_runs), max_runs) ||
    !identical(as.numeric(autok$convergence_tol), convergence_tol) ||
    !identical(as.integer(autok$seed), seed)) {
  stop("Returned AutoK metadata does not match the fixed benchmark contract.", call. = FALSE)
}
selected_model <- CellCharterR::cellcharter_autok_model(autok)
selected_labels <- selected_model$labels
if (length(selected_labels) != ncol(sce)) {
  stop("Selected model did not return one label per mouse-spleen cell.", call. = FALSE)
}

# Fit the tutorial-matched K=11 sensitivity separately with its published seed.
k11_seed <- 12345L
k11_started <- Sys.time()
k11_model <- CellCharterR::cpp_gmm_fit(
  x,
  n_clusters = 11L,
  max_iter = 100L,
  tolerance = 0.001,
  covariance_regularization = 1e-6,
  seed = k11_seed,
  num_threads = threads_per_fit
)
k11_finished <- Sys.time()
k11_elapsed <- as.numeric(difftime(k11_finished, k11_started, units = "secs"))
if (length(k11_model$labels) != ncol(sce)) {
  stop("The K=11 sensitivity fit did not return one label per cell.", call. = FALSE)
}

best_row <- autok$best_index[autok$best_index$k == autok$best_k, , drop = FALSE]
if (nrow(best_row) != 1L) stop("AutoK selected-model index is invalid.", call. = FALSE)
fits_completed <- nrow(autok$fit_index)
expected_fits_completed <- length(autok$fit_k) * autok$runs_completed
if (!identical(as.integer(fits_completed), as.integer(expected_fits_completed))) {
  stop("AutoK fit index is incomplete for its reported rounds.", call. = FALSE)
}
resumed_fits <- sum(autok$fit_index$resumed %in% TRUE)
if (resumed_fits > 0L) {
  warning(
    "AutoK resumed durable fits; elapsed and memory measurements cover only this invocation.",
    call. = FALSE
  )
}
final_stability_change <- if (nrow(autok$convergence_history)) {
  utils::tail(autok$convergence_history$stability_change, 1L)
} else {
  NA_real_
}

runner_finished <- Sys.time()
process_peak_at_end_mb <- .peak_process_memory_mb()
provenance <- list(
  schema_version = 2L,
  benchmark_role = "compiled C++ CPU AutoK; no GPU or torch execution",
  backend = "compiled-cpp-double",
  backend_requested = backend,
  backend_algorithm = autok$resume_contract$algorithm,
  implementation_digest = autok$implementation_digest,
  reference_version = autok$reference_version,
  reference_commit = autok$reference_commit,
  # Retain these top-level fields for the existing downstream artifact builder.
  input_path = input_path,
  input_sha256 = input_sha256,
  representation_sha256 = autok$representation_sha256,
  dimensions = dim(x),
  workers = workers,
  threads_per_fit = threads_per_fit,
  seed = seed,
  max_runs = max_runs,
  runs_completed = autok$runs_completed,
  elapsed_seconds = autok_elapsed,
  selected_k = autok$best_k,
  selected_nll = selected_model$nll,
  k11_seed = k11_seed,
  k11_nll = k11_model$nll,
  k11_elapsed_seconds = k11_elapsed,
  implementation = list(
    autok_sha256 = autok$implementation_digest,
    similarity_sha256 = autok$similarity_digest,
    source_sha256 = source_sha256,
    native_library_path = native_dll_path,
    native_library_sha256 = .sha256_if_file(native_dll_path),
    package_version = as.character(utils::packageVersion("CellCharterR")),
    git_commit = paste(.git_output(c("rev-parse", "HEAD")), collapse = "\n"),
    git_status_porcelain = .git_output(c("status", "--short"))
  ),
  initializer = selected_model$initialization$method,
  initializer_limitation = paste(
    "Deterministic C++ K-means++/Lloyd; not sklearn RNG-identical.",
    "The EM and scoring equations are the pinned two-pass full-covariance target."
  ),
  input = list(
    path = input_path,
    sha256 = input_sha256,
    numeric_representation_sha256 = autok$representation_sha256,
    dimensions = dim(x),
    cell_order_sha256 = cell_order_sha256,
    aggregation_manifest_path = if (file.exists(manifest_path)) manifest_path else NA_character_,
    aggregation_manifest_sha256 = aggregation_manifest_sha256,
    aggregation_data_id = aggregation_manifest$data_id,
    aggregation_graph_source = aggregation_manifest$analysis_graph_source,
    aggregation_sample_sha256 = aggregation_manifest$sample_sha256,
    aggregation_declared_sce_path = aggregation_manifest$sce_checkpoint,
    aggregation_declared_sce_path_matches = if (
      is.null(aggregation_manifest$sce_checkpoint)
    ) {
      NA
    } else {
      identical(
        normalizePath(
          aggregation_manifest$sce_checkpoint,
          winslash = "/",
          mustWork = FALSE
        ),
        input_path
      )
    }
  ),
  scientific_parameters = list(
    candidate_k = candidate_k,
    fitted_k = autok$fit_k,
    max_runs = max_runs,
    convergence_tol = convergence_tol,
    seed = seed,
    covariance_type = "full",
    max_iter = 100L,
    tolerance = 0.001,
    covariance_regularization = 1e-6
  ),
  worker_topology = list(
    outer_workers = workers,
    threads_per_fit = threads_per_fit,
    maximum_concurrent_fit_threads = workers * threads_per_fit,
    detected_logical_cores = logical_cores,
    detected_physical_cores = physical_cores,
    blas_openmp_environment_before_package_load = .thread_environment,
    blas_openmp_environment_all_one = !any(.unsafe_thread_environment),
    device = "cpu",
    gpu_used = FALSE,
    torch_namespace_loaded = "torch" %in% loadedNamespaces()
  ),
  autok = list(
    converged = autok$converged,
    runs_completed = autok$runs_completed,
    fits_per_run = length(autok$fit_k),
    fits_completed = fits_completed,
    fits_resumed = resumed_fits,
    fits_new_this_invocation = fits_completed - resumed_fits,
    convergence_history = autok$convergence_history,
    final_stability_change = final_stability_change,
    selected_k = autok$best_k,
    selected_run = best_row$run[[1L]],
    selected_seed = best_row$seed[[1L]],
    selected_nll = selected_model$nll
  ),
  sensitivity_k11 = list(
    k = 11L,
    seed = k11_seed,
    nll = k11_model$nll,
    elapsed_seconds = k11_elapsed
  ),
  timing = list(
    measurement_scope = if (resumed_fits) {
      "current invocation; not a complete campaign benchmark because durable fits were resumed"
    } else {
      "complete fresh campaign invocation"
    },
    runner_started_utc = format(runner_started, tz = "UTC", usetz = TRUE),
    autok_started_utc = format(autok_started, tz = "UTC", usetz = TRUE),
    autok_finished_utc = format(autok_finished, tz = "UTC", usetz = TRUE),
    autok_elapsed_seconds = autok_elapsed,
    k11_elapsed_seconds = k11_elapsed,
    total_elapsed_seconds = as.numeric(
      difftime(runner_finished, runner_started, units = "secs")
    ),
    finished_utc = format(runner_finished, tz = "UTC", usetz = TRUE)
  ),
  memory = list(
    measurement = "Linux VmHWM process high-water mark and R gc() high-water estimate",
    process_peak_at_runner_start_mb = process_peak_at_start_mb,
    process_peak_mb = process_peak_at_end_mb,
    process_peak_increment_lower_bound_mb = if (
      is.finite(process_peak_at_start_mb) && is.finite(process_peak_at_end_mb)
    ) max(0, process_peak_at_end_mb - process_peak_at_start_mb) else NA_real_,
    r_heap_peak_mb = .peak_r_memory_mb()
  ),
  session_info = utils::sessionInfo()
)

artifacts <- c(
  autok = file.path(output_dir, "autok_result.rds"),
  selected_model = file.path(output_dir, "selected_model.rds"),
  selected_labels = file.path(output_dir, "selected_labels.rds"),
  tutorial_k11_model = file.path(output_dir, "tutorial_k11_model.rds"),
  tutorial_k11_labels = file.path(output_dir, "tutorial_k11_labels.rds")
)
.atomic_save_rds(autok, artifacts[["autok"]], compress = FALSE)
.atomic_save_rds(selected_model, artifacts[["selected_model"]], compress = FALSE)
.atomic_save_rds(
  list(
    k = autok$best_k,
    labels = selected_labels,
    seed = seed,
    selected_fit_seed = best_row$seed[[1L]],
    cell_ids = cell_ids,
    cell_order_sha256 = cell_order_sha256,
    input_sha256 = input_sha256,
    role_id = "r_autok_primary"
  ),
  artifacts[["selected_labels"]],
  compress = FALSE
)
.atomic_save_rds(k11_model, artifacts[["tutorial_k11_model"]], compress = FALSE)
.atomic_save_rds(
  list(
    k = 11L,
    labels = k11_model$labels,
    seed = k11_seed,
    cell_ids = cell_ids,
    cell_order_sha256 = cell_order_sha256,
    input_sha256 = input_sha256,
    role_id = "r_tutorial_k11_sensitivity"
  ),
  artifacts[["tutorial_k11_labels"]],
  compress = FALSE
)
provenance$artifacts <- lapply(artifacts, function(path) list(
  path = path,
  bytes = unname(file.info(path)$size),
  sha256 = CellCharterR::cellcharter_file_sha256(path)
))
# Include publication/checksum overhead and any serialization high-water mark
# in the runner totals, while retaining the isolated AutoK fit elapsed time.
runner_finished <- Sys.time()
process_peak_at_end_mb <- .peak_process_memory_mb()
provenance$timing$total_elapsed_seconds <- as.numeric(
  difftime(runner_finished, runner_started, units = "secs")
)
provenance$timing$finished_utc <- format(runner_finished, tz = "UTC", usetz = TRUE)
provenance$memory$process_peak_mb <- process_peak_at_end_mb
provenance$memory$process_peak_increment_lower_bound_mb <- if (
  is.finite(process_peak_at_start_mb) && is.finite(process_peak_at_end_mb)
) max(0, process_peak_at_end_mb - process_peak_at_start_mb) else NA_real_
provenance$memory$r_heap_peak_mb <- .peak_r_memory_mb()
# Publishing provenance last is the completion marker for this result bundle.
.atomic_save_rds(
  provenance,
  file.path(output_dir, "provenance.rds"),
  compress = FALSE
)

cat(sprintf(
  paste0(
    "Compiled CPU AutoK complete: selected K=%d; %d fits across %d runs; ",
    "converged=%s; AutoK %.1f minutes; K=11 sensitivity %.1f minutes.\n",
    "Results: %s\n"
  ),
  autok$best_k,
  fits_completed,
  autok$runs_completed,
  autok$converged,
  autok_elapsed / 60,
  k11_elapsed / 60,
  output_dir
))
