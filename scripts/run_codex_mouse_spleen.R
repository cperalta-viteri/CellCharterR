#!/usr/bin/env Rscript

# Checkpointed R reproduction of the CellCharter CODEX mouse-spleen workflow.
# Large matrices stay HDF5-backed and the analysis state is orchestrated by a
# SingleCellExperiment. Set CODEX_STAGES to a comma-separated subset to resume
# particular stages
# (representation,graph,aggregation,autok,tutorial_k11,downstream,report).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(HDF5Array)
  library(Matrix)
})
source("R/sce_workflow.R")
source("R/aggregation.R")
source("R/torch_gmm.R")
source("R/provenance.R")
source("R/autok.R")
source("R/enrichment.R")

input_dir <- normalizePath(
  Sys.getenv("CODEX_INPUT_DIR", file.path("results", "codex_mouse_spleen", "input")),
  winslash = "/", mustWork = FALSE
)
result_dir <- normalizePath(
  Sys.getenv(
    "CODEX_RESULTS_DIR",
    file.path("results", "CellCharterR_mouse_spleen_v2")
  ),
  winslash = "/", mustWork = FALSE
)
frozen_result_roots <- normalizePath(c(
  file.path("results", "codex_mouse_spleen"),
  file.path("results", "CellCharterR_mouse_spleen")
), winslash = "/", mustWork = FALSE)
if (any(result_dir == frozen_result_roots |
        startsWith(result_dir, paste0(frozen_result_roots, "/")))) {
  stop(
    "The requested CODEX_RESULTS_DIR is a frozen first-round/Phase-0 result ",
    "root. Use a new path such as results/CellCharterR_mouse_spleen_v2."
  )
}
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

requested_stages <- strsplit(
  Sys.getenv(
    "CODEX_STAGES",
    "representation,graph,aggregation,autok,tutorial_k11,downstream,report"
  ),
  ",", fixed = TRUE
)[[1L]]
requested_stages <- trimws(requested_stages)
known_stages <- c(
  "representation", "graph", "aggregation", "autok", "tutorial_k11",
  "downstream", "report"
)
if (any(!requested_stages %in% known_stages)) {
  stop("Unknown CODEX_STAGES value(s): ",
       paste(setdiff(requested_stages, known_stages), collapse = ", "))
}

# Load Torch before opening the large collection of HDF5 and graph checkpoints.
# On some Linux builds processx/ps initialization can fail late in a long-lived
# R process even though the same namespace loads reliably at process startup.
# Loading only for stages that require fitted models preserves lightweight
# graph/aggregation-only resumes and gives an actionable error before any fit.
torch_needed <- any(requested_stages %in% c(
  "autok", "tutorial_k11", "downstream", "report"
))
torch_load_error <- NULL
torch_loaded <- !torch_needed || tryCatch({
  loadNamespace("torch")
  TRUE
}, error = function(error) {
  torch_load_error <<- conditionMessage(error)
  FALSE
})
if (!torch_loaded) {
  stop(
    "R `torch` cannot be loaded before the model stages. Earlier checkpoints ",
    "remain safe. Load error: ", torch_load_error
  )
}

raw_h5ad <- file.path(input_dir, "codex_mouse_spleen.h5ad")
clustered_h5ad <- file.path(input_dir, "codex_mouse_spleen_clustered.h5ad")
official_autok_zip <- file.path(input_dir, "BALBc-MRL_trVAE_CellCharter_autok.zip")
if (!all(file.exists(c(raw_h5ad, clustered_h5ad)))) {
  stop(
    "Official inputs are absent. Run `Rscript scripts/download_codex_mouse_spleen.R` ",
    "or set CODEX_INPUT_DIR before restarting this workflow."
  )
}

seed <- as.integer(Sys.getenv("CODEX_SEED", "0"))
tutorial_k11_seed <- 12345L
batch_size <- as.integer(Sys.getenv("CODEX_BATCH_SIZE", "8192"))
aggregation_block_rows <- .cc_aggregation_block_rows(suppressWarnings(as.numeric(
  Sys.getenv("CODEX_AGGREGATION_BLOCK_ROWS", "8192")
)))
num_threads <- as.integer(Sys.getenv(
  "CODEX_NUM_THREADS",
  as.character(max(1L, parallel::detectCores(logical = FALSE)))
))
autok_bounds <- c(
  as.integer(Sys.getenv("CODEX_MIN_K", "2")),
  as.integer(Sys.getenv("CODEX_MAX_K", "20"))
)
max_runs <- as.integer(Sys.getenv("CODEX_MAX_RUNS", "10"))
convergence_tol <- as.numeric(Sys.getenv("CODEX_AUTOK_CONVERGENCE_TOL", "0.01"))
n_perms <- as.integer(Sys.getenv("CODEX_N_PERMS", "100"))
gmm_max_iter <- as.integer(Sys.getenv("CODEX_GMM_MAX_ITER", "100"))
gmm_tolerance <- as.numeric(Sys.getenv("CODEX_GMM_TOLERANCE", "0.001"))
covariance_regularization <- as.numeric(Sys.getenv("CODEX_COVARIANCE_REGULARIZATION", "0.000001"))
graph_source <- Sys.getenv("CODEX_GRAPH_SOURCE", "r_pinned_global")
if (!graph_source %in% c("r_published", "r_pinned_global", "official")) {
  stop("CODEX_GRAPH_SOURCE must be `r_published`, `r_pinned_global`, or `official`.")
}

.atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), dirname(path), ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  # The SCE stores only HDF5 seeds and metadata here. Base serialization keeps
  # those lazy references without asking S4Vectors to realize out-of-memory
  # assays; the referenced HDF5 files are separately checksummed in metadata.
  base::saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("Could not checkpoint: ", path)
  invisible(path)
}

.peak_r_memory_mb <- function() {
  current <- gc()
  as.numeric(current["Ncells", "max used"] * 56 + current["Vcells", "max used"] * 8) / 1024^2
}

.peak_process_memory_mb <- function() {
  # Linux VmHWM is the measured peak resident set for the complete process,
  # including native/HDF5 allocations that R's garbage collector cannot see.
  if (!file.exists("/proc/self/status")) return(NA_real_)
  status <- readLines("/proc/self/status", warn = FALSE)
  line <- status[grepl("^VmHWM:", status)]
  if (length(line) != 1L) return(NA_real_)
  value_kb <- suppressWarnings(as.numeric(sub(
    "^VmHWM:[[:space:]]*([0-9]+)[[:space:]]+kB.*$", "\\1", line
  )))
  if (length(value_kb) != 1L || !is.finite(value_kb)) NA_real_ else value_kb / 1024
}

.record_stage <- function(stage, started, extra = list()) {
  path <- file.path(result_dir, "provenance", paste0(stage, ".rds"))
  payload <- c(list(
    stage = stage,
    started_utc = format(started, tz = "UTC", usetz = TRUE),
    finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    peak_r_memory_mb = .peak_r_memory_mb(),
    session_info = utils::sessionInfo()
  ), extra)
  .atomic_save_rds(payload, path)
}

.workflow_source <- "scripts/run_codex_mouse_spleen.R"
.contract_path <- function(stage) {
  file.path(result_dir, "provenance", "contracts", paste0(stage, ".rds"))
}
.installed_version <- function(package) {
  if (identical(package, "CellCharterR")) {
    return(as.character(read.dcf("DESCRIPTION", fields = "Version")[[1L]]))
  }
  if (!requireNamespace(package, quietly = TRUE)) return("not-installed")
  as.character(utils::packageVersion(package))
}
.stage_versions <- function(packages = character()) {
  values <- c(
    list(
      CellCharterR = .installed_version("CellCharterR"),
      R = paste(R.version$major, R.version$minor, sep = "."),
      reference_cellcharter = "0.2.0",
      reference_commit = "2a8fdde8213e32e5a00f1f9571bffd3dbb825a79",
      reference_squidpy = "1.3.0",
      reference_squidpy_commit = "e2ffc4d03d02e1b02c9bf15ae87fb7091a64df4b",
      reference_pycave = "3.2.1",
      reference_pycave_commit = "6df56fa15797940a2063f3e70b3a7641f592a9f6",
      reference_pytorch = "1.12.1",
      reference_python = "3.10"
    ),
    stats::setNames(lapply(packages, .installed_version), packages)
  )
  values[!duplicated(names(values))]
}
.stage_implementation <- function(stage, functions) {
  cellcharter_implementation_digest(
    functions = functions,
    source_regions = stats::setNames(list(list(
      path = .workflow_source,
      begin = paste0("# CHECKPOINT CONTRACT REGION: ", stage, " BEGIN"),
      end = paste0("# CHECKPOINT CONTRACT REGION: ", stage, " END")
    )), stage)
  )
}
.require_contract <- function(stage, expected, artifacts, allow_verify = FALSE,
                              evidence = list()) {
  path <- .contract_path(stage)
  if (cellcharter_checkpoint_reusable(path, expected, artifacts)) return(TRUE)
  if (!file.exists(path) && isTRUE(allow_verify) && all(file.exists(artifacts))) {
    cellcharter_record_checkpoint(
      expected, artifacts, status = "verified_existing_checkpoint",
      evidence = evidence, path = path
    )
    return(TRUE)
  }
  FALSE
}

.hdf5_seed_identity <- function(value) {
  seed <- try(DelayedArray::seed(value), silent = TRUE)
  if (inherits(seed, "try-error") || !inherits(seed, "HDF5ArraySeed")) {
    return(NULL)
  }
  list(
    filepath = normalizePath(methods::slot(seed, "filepath"), winslash = "/",
                             mustWork = TRUE),
    dataset = sub("^/", "", methods::slot(seed, "name")),
    dimensions = as.integer(dim(value))
  )
}

.representation_checkpoint_matches_inputs <- function(sce, raw_h5ad,
                                                        clustered_h5ad) {
  inputs <- metadata(sce)$CellCharterR$inputs
  required_inputs <- c(
    "raw_h5ad", "raw_sha256", "clustered_h5ad", "clustered_sha256",
    "representation_dataset"
  )
  if (!is.list(inputs) || !all(required_inputs %in% names(inputs)) ||
      any(vapply(inputs[required_inputs], is.null, logical(1)))) {
    return(FALSE)
  }
  expected <- list(
    raw = list(
      filepath = normalizePath(raw_h5ad, winslash = "/", mustWork = TRUE),
      sha256 = cellcharter_file_sha256(raw_h5ad)
    ),
    clustered = list(
      filepath = normalizePath(clustered_h5ad, winslash = "/", mustWork = TRUE),
      sha256 = cellcharter_file_sha256(clustered_h5ad)
    )
  )
  protein <- .hdf5_seed_identity(assay(sce, "protein"))
  representation <- .hdf5_seed_identity(reducedDim(sce, "X_trVAE"))
  spatial <- .hdf5_seed_identity(reducedDim(sce, "spatial"))
  identical(normalizePath(inputs$raw_h5ad, winslash = "/", mustWork = TRUE),
              expected$raw$filepath) &&
    identical(inputs$raw_sha256, expected$raw$sha256) &&
    identical(normalizePath(inputs$clustered_h5ad, winslash = "/", mustWork = TRUE),
              expected$clustered$filepath) &&
    identical(inputs$clustered_sha256, expected$clustered$sha256) &&
    identical(inputs$representation_dataset, "obsm/X_trVAE") &&
    identical(protein, list(
      filepath = expected$raw$filepath, dataset = "X",
      dimensions = as.integer(dim(sce))
    )) &&
    identical(representation, list(
      filepath = expected$clustered$filepath, dataset = "obsm/X_trVAE",
      dimensions = c(ncol(sce), 10L)
    )) &&
    identical(spatial, list(
      filepath = expected$raw$filepath, dataset = "obsm/spatial",
      dimensions = c(ncol(sce), 2L)
    ))
}

.sample_path <- function(directory, sample_id, suffix = ".rds") {
  file.path(directory, paste0(gsub("[^A-Za-z0-9_.-]", "_", sample_id), suffix))
}

.read_h5ad_csr_block <- function(filepath, group, rows) {
  rows <- as.integer(rows)
  if (!length(rows) || !identical(rows, seq.int(rows[[1L]], rows[[length(rows)]]))) {
    stop("H5AD sparse block rows must be one non-empty contiguous run.")
  }
  pointer <- as.integer(rhdf5::h5read(
    filepath, paste0(group, "/indptr"),
    index = list(seq.int(rows[[1L]], rows[[length(rows)]] + 1L))
  ))
  first_entry <- pointer[[1L]]
  last_entry <- pointer[[length(pointer)]]
  if (last_entry == first_entry) {
    column <- integer()
  } else {
    column <- as.integer(rhdf5::h5read(
      filepath, paste0(group, "/indices"),
      index = list(seq.int(first_entry + 1L, last_entry))
    )) - rows[[1L]] + 2L
  }
  row <- rep.int(seq_along(rows), diff(pointer))
  if (length(column) && (any(column < 1L) || any(column > length(rows)))) {
    stop("Published graph contains an edge across physical samples.")
  }
  graph <- Matrix::sparseMatrix(
    i = row, j = column, x = 1, dims = rep.int(length(rows), 2L),
    dimnames = list(colnames(sce)[rows], colnames(sce)[rows])
  )
  Matrix::drop0(graph > 0)
}

.validate_graph_identity <- function(graph, expected_ids, context) {
  .cc_align_observation_ids(
    rownames(graph), expected_ids, paste0(context, " row identifiers"),
    "SCE sample cell identifiers"
  )
  .cc_align_observation_ids(
    colnames(graph), expected_ids, paste0(context, " column identifiers"),
    "SCE sample cell identifiers"
  )
  invisible(TRUE)
}

.aggregation_checkpoint_valid <- function(path, expected_ids, sample_id, block_rows,
                                          contract = NULL) {
  complete_path <- paste0(path, ".complete.rds")
  if (!file.exists(path) || !file.exists(complete_path)) return(FALSE)
  marker <- try(readRDS(complete_path), silent = TRUE)
  existing <- try(HDF5Array::HDF5Array(path, "X_cellcharter"), silent = TRUE)
  if (inherits(marker, "try-error") || inherits(existing, "try-error")) return(FALSE)
  marker_ids <- try(
    .cc_validate_observation_ids(
      marker$cell_ids, paste0("Aggregation checkpoint `", sample_id, "` cell identifiers")
    ), silent = TRUE
  )
  if (inherits(marker_ids, "try-error")) return(FALSE)
  aligned <- try(
    .cc_align_observation_ids(
      marker_ids, expected_ids,
      paste0("Aggregation checkpoint `", sample_id, "` cell identifiers"),
      "SCE sample cell identifiers"
    ), silent = TRUE
  )
  valid_memory_provenance <-
    is.numeric(marker$estimated_dense_workspace_bytes) &&
    length(marker$estimated_dense_workspace_bytes) == 1L &&
    is.finite(marker$estimated_dense_workspace_bytes) &&
    marker$estimated_dense_workspace_bytes >= 0 &&
    is.numeric(marker$measured_r_heap_peak_after_mb) &&
    length(marker$measured_r_heap_peak_after_mb) == 1L &&
    is.finite(marker$measured_r_heap_peak_after_mb) &&
    marker$measured_r_heap_peak_after_mb >= 0 &&
    is.numeric(marker$measured_process_peak_after_mb) &&
    length(marker$measured_process_peak_after_mb) == 1L &&
    (is.na(marker$measured_process_peak_after_mb) ||
       marker$measured_process_peak_after_mb >= 0)
  valid <- !inherits(aligned, "try-error") && identical(marker$sample_id, sample_id) &&
    identical(marker$cell_order_sha256, .cc_observation_order_sha256(expected_ids)) &&
    identical(marker$block_rows, block_rows) &&
    valid_memory_provenance &&
    identical(as.integer(dim(existing)), c(length(expected_ids), 40L)) &&
    identical(marker$sha256, cellcharter_file_sha256(path))
  if (!valid || is.null(contract)) return(valid)
  cellcharter_checkpoint_reusable(
    paste0(path, ".contract.rds"), contract,
    c(hdf5 = path, completion_marker = complete_path)
  )
}

# CHECKPOINT CONTRACT REGION: representation BEGIN
representation_dir <- file.path(result_dir, "representation")
sce_path <- file.path(representation_dir, "codex_sce.rds")
representation_ids <- .cc_h5ad_index(raw_h5ad)
representation_order_sha256 <- .cc_observation_order_sha256(representation_ids)
representation_contract <- cellcharter_stage_contract(
  stage = "representation",
  implementation_digest = .stage_implementation(
    "representation", list(
      read_codex_sce = read_codex_sce,
      read_h5ad_dataframe = .cc_h5ad_dataframe,
      read_h5ad_index = .cc_h5ad_index,
      validate_observation_ids = .cc_validate_observation_ids,
      observation_order_sha256 = .cc_observation_order_sha256,
      align_observation_ids = .cc_align_observation_ids,
      validate_sce_identity = .cc_validate_sce_identity,
      codex_h5ad_identity = .cc_codex_h5ad_identity,
      hdf5_seed_identity = .hdf5_seed_identity,
      verify_legacy_representation = .representation_checkpoint_matches_inputs
    )
  ),
  upstream = list(
    raw_h5ad_sha256 = cellcharter_file_sha256(raw_h5ad),
    clustered_h5ad_sha256 = cellcharter_file_sha256(clustered_h5ad)
  ),
  cell_order_sha256 = representation_order_sha256,
  parameters = list(
    representation_dataset = "obsm/X_trVAE",
    spatial_dataset = "obsm/spatial",
    raw_cluster_identity = "sample+local_index",
    transpose_to_feature_by_cell = TRUE
  ),
  seeds = list(seed = NA_integer_),
  versions = .stage_versions(c(
    "SingleCellExperiment", "HDF5Array", "rhdf5", "DelayedArray"
  )),
  runtime = list(storage = "HDF5-backed reduced dimensions")
)
representation_created <- FALSE
if ("representation" %in% requested_stages && !file.exists(sce_path)) {
  started <- Sys.time()
  sce <- read_codex_sce(raw_h5ad, clustered_h5ad)
  metadata(sce)$CellCharterR$parameters <- list(
    seed = seed, autok_bounds = autok_bounds, max_runs = max_runs,
    convergence_tol = convergence_tol, batch_size = batch_size,
    aggregation_block_rows = aggregation_block_rows,
    num_threads = num_threads, n_perms = n_perms,
    gmm_max_iter = gmm_max_iter, gmm_tolerance = gmm_tolerance,
    covariance_regularization = covariance_regularization,
    graph_source = graph_source
  )
  .atomic_save_rds(sce, sce_path)
  representation_created <- TRUE
  .record_stage("representation", started, list(
    checkpoint = sce_path,
    cells = ncol(sce), markers = nrow(sce), trvae_dimensions = ncol(reducedDim(sce, "X_trVAE"))
  ))
}
if (!file.exists(sce_path)) stop("Representation checkpoint is absent: ", sce_path)
sce <- readRDS(sce_path)
# Migrate first-round checkpoints in memory after the package rename. The old
# key remains readable so durable expensive artifacts do not need regeneration.
if (is.null(metadata(sce)$CellCharterR) && !is.null(metadata(sce)$rCellCharter)) {
  metadata(sce)$CellCharterR <- metadata(sce)$rCellCharter
}
if (is.null(metadata(sce)$CellCharterR$inputs$cell_order_sha256)) {
  stop(
    "Representation SCE checkpoint lacks `inputs$cell_order_sha256`. Regenerate ",
    "with `CODEX_STAGES=representation` under a new CODEX_RESULTS_DIR."
  )
}
.cc_validate_sce_identity(
  sce, "Representation SCE checkpoint", c("X_trVAE", "spatial"),
  metadata(sce)$CellCharterR$inputs$cell_order_sha256
)
invisible(.cc_align_observation_ids(
  representation_ids, colnames(sce),
  "Current raw CODEX H5AD observation identifiers",
  "representation SCE checkpoint cell identifiers"
))
# Older representation checkpoints are upgraded only when fields are absent.
# Never replace creation-time parameters with values from a later resume or
# report-only invocation; stage-specific provenance records the latter.
requested_parameters <- list(
  seed = seed, autok_bounds = autok_bounds, max_runs = max_runs,
  convergence_tol = convergence_tol, batch_size = batch_size,
  aggregation_block_rows = aggregation_block_rows,
  num_threads = num_threads, n_perms = n_perms,
  gmm_max_iter = gmm_max_iter, gmm_tolerance = gmm_tolerance,
  covariance_regularization = covariance_regularization,
  graph_source = graph_source
)
reference_updated <- is.null(metadata(sce)$CellCharterR$reference$scipy_version)
if (reference_updated) {
  metadata(sce)$CellCharterR$reference$scipy_version <-
    "not declared by cellcharter_analyses environment"
}
autok_provenance_updated <- file.exists(official_autok_zip) &&
  is.null(metadata(sce)$CellCharterR$inputs$official_autok_sha256)
if (autok_provenance_updated) {
  metadata(sce)$CellCharterR$inputs$official_autok_zip <-
    normalizePath(official_autok_zip)
  metadata(sce)$CellCharterR$inputs$official_autok_sha256 <-
    cellcharter_file_sha256(official_autok_zip)
}
labels_updated <- !"reference_cluster_cellcharter" %in% colnames(colData(sce))
if (labels_updated) {
  clustered_obs <- .cc_h5ad_dataframe(clustered_h5ad, "obs")
  current_identity <- .cc_codex_h5ad_identity(
    raw_h5ad, clustered_h5ad,
    .cc_h5ad_dataframe(raw_h5ad, "obs"), clustered_obs
  )
  .cc_align_observation_ids(
    current_identity$raw_ids, colnames(sce),
    "Current raw CODEX H5AD observation identifiers",
    "representation SCE checkpoint cell identifiers"
  )
  if (nrow(clustered_obs) != ncol(sce) ||
      !"cluster_cellcharter" %in% colnames(clustered_obs)) {
    stop("Cannot attach the published H5AD `cluster_cellcharter` reference labels.")
  }
  colData(sce)$reference_cluster_cellcharter <-
    clustered_obs$cluster_cellcharter
}
parameters_updated <- is.null(metadata(sce)$CellCharterR$parameters)
if (parameters_updated) {
  metadata(sce)$CellCharterR$parameters <- requested_parameters
  metadata(sce)$CellCharterR$parameters_provenance <-
    "Backfilled from the first invocation that upgraded this legacy checkpoint."
}
# Legacy additions above are deliberately in-memory only. The immutable SCE
# checkpoint is never rewritten during validation/migration.
representation_artifacts <- c(sce = sce_path)
if (representation_created) {
  cellcharter_record_checkpoint(
    representation_contract, representation_artifacts,
    status = "created_under_contract",
    evidence = list(
      cell_identity_validation = "raw and clustered H5AD sample+local_index exact alignment",
      reduced_dimensions = c("X_trVAE", "spatial")
    ),
    path = .contract_path("representation")
  )
} else if (!.require_contract(
  "representation", representation_contract, representation_artifacts,
  allow_verify = !file.exists(.contract_path("representation")) &&
    .representation_checkpoint_matches_inputs(sce, raw_h5ad, clustered_h5ad),
  evidence = list(
    verification_kind = "read-only validation of pre-contract checkpoint",
    original_creation_contract = "unknown",
    verified_cell_order_sha256 = representation_order_sha256,
    verified_input_alignment = "raw H5AD IDs exactly equal checkpoint column IDs",
    verified_reduced_dimensions = c("X_trVAE", "spatial"),
    checkpoint_mtime_utc = format(file.info(sce_path)$mtime, tz = "UTC", usetz = TRUE)
  )
)) {
  stop(
    "Representation checkpoint contract is stale or invalid. Its data were not ",
    "modified; use a new CODEX_RESULTS_DIR to regenerate it."
  )
}
representation_provenance_path <- file.path(result_dir, "provenance", "representation.rds")
if (!file.exists(representation_provenance_path)) {
  # A resumed legacy checkpoint may predate stage-level provenance. Record the
  # verification event truthfully instead of inventing the original duration.
  verified_at <- Sys.time()
  .record_stage("representation", verified_at, list(
    checkpoint = sce_path,
    checkpoint_mtime_utc = format(file.info(sce_path)$mtime, tz = "UTC", usetz = TRUE),
    provenance_recovered_from_existing_checkpoint = TRUE,
    original_elapsed_seconds = NA_real_,
    cells = ncol(sce), markers = nrow(sce),
    trvae_dimensions = ncol(reducedDim(sce, "X_trVAE"))
  ))
}
if (!any(requested_stages %in% c(
  "graph", "aggregation", "autok", "tutorial_k11", "downstream", "report"
))) {
  message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
  quit(save = "no", status = 0L)
}
# CHECKPOINT CONTRACT REGION: representation END

# CHECKPOINT CONTRACT REGION: graph BEGIN
sample_vector <- .cc_validate_observation_ids(
  as.character(colData(sce)$sample), "Representation SCE sample identifiers",
  allow_duplicates = TRUE
)
sample_ids <- if (is.factor(colData(sce)$sample)) {
  levels(droplevels(colData(sce)$sample))
} else {
  unique(sample_vector)
}
rows_by_sample <- stats::setNames(
  lapply(sample_ids, function(id) which(sample_vector == id)), sample_ids
)
rows_path <- file.path(result_dir, "graph", "rows_by_sample.rds")

graph_dir <- file.path(result_dir, "graph")
untrimmed_dir <- file.path(graph_dir, "untrimmed")
trimmed_dir <- file.path(graph_dir, "trimmed")
graph_manifest_path <- file.path(graph_dir, "graph_manifest.rds")
graph_manifest_created <- FALSE
graph_implementation <- .stage_implementation("graph", list(
  delaunay_graph = delaunay_graph,
  remove_long_links = remove_long_links,
  global_distance_cutoff = cellcharter_global_distance_cutoff,
  coordinates = .cc_coordinates,
  adjacency_topology = .cc_adjacency_topology,
  read_h5ad_csr_block = .read_h5ad_csr_block
))
.graph_equal <- function(left, right) {
  identical(dim(left), dim(right)) && identical(dimnames(left), dimnames(right)) &&
    Matrix::nnzero(Matrix::drop0(left - right)) == 0L &&
    identical(attr(left, "distance_cutoff"), attr(right, "distance_cutoff"))
}
.graph_contract <- function(sample_id, kind, expected_ids, coordinates,
                            upstream, cutoff = NA_real_) {
  cellcharter_stage_contract(
    stage = paste("graph", kind, sample_id, sep = "/"),
    implementation_digest = graph_implementation,
    upstream = upstream,
    cell_order_sha256 = .cc_observation_order_sha256(expected_ids),
    parameters = list(
      sample_id = sample_id, kind = kind,
      coordinates_sha256 = digest::digest(coordinates, algo = "sha256"),
      qhull_options = "Qbb Qc Qz Q12 Qt", percentile = 99,
      cutoff = cutoff, cross_sample_edges = FALSE
    ),
    seeds = list(seed = NA_integer_),
    versions = .stage_versions(c("geometry", "Matrix")), runtime = list()
  )
}
.record_or_verify_graph <- function(path, contract, candidate, context, created) {
  record_path <- paste0(path, ".contract.rds")
  artifacts <- c(graph = path)
  if (created) {
    cellcharter_record_checkpoint(
      contract, artifacts, status = "created_under_contract",
      evidence = list(validation = context), path = record_path
    )
  } else if (!cellcharter_checkpoint_reusable(record_path, contract, artifacts)) {
    if (file.exists(record_path) || !.graph_equal(readRDS(path), candidate)) {
      stop(context, " has a stale contract or fails exact read-only reconstruction.")
    }
    cellcharter_record_checkpoint(
      contract, artifacts, status = "verified_existing_checkpoint",
      evidence = list(
        verification_kind = "read-only exact algorithmic reconstruction",
        original_creation_contract = "unknown",
        exact_sparse_matrix_match = TRUE,
        verified_artifact_sha256 = cellcharter_file_sha256(path)
      ), path = record_path
    )
  }
  invisible(TRUE)
}

graph_checkpoint_complete <- file.exists(graph_manifest_path)
if (graph_checkpoint_complete) {
  existing_graph_manifest <- try(readRDS(graph_manifest_path), silent = TRUE)
  graph_checkpoint_complete <- !inherits(existing_graph_manifest, "try-error") &&
    is.character(existing_graph_manifest$trimmed_paths) &&
    identical(names(existing_graph_manifest$trimmed_paths), sample_ids) &&
    all(file.exists(existing_graph_manifest$trimmed_paths)) &&
    is.character(existing_graph_manifest$r_published_paths) &&
    all(file.exists(existing_graph_manifest$r_published_paths)) &&
    is.character(existing_graph_manifest$r_pinned_global_paths) &&
    all(file.exists(existing_graph_manifest$r_pinned_global_paths)) &&
    is.character(existing_graph_manifest$official_paths) &&
    all(file.exists(existing_graph_manifest$official_paths))
}

if ("graph" %in% requested_stages && !graph_checkpoint_complete) {
  started <- Sys.time()
  untrimmed_dir <- file.path(graph_dir, "r_qhull_untrimmed_v1")
  global_trimmed_dir <- file.path(graph_dir, "r_qhull_trimmed_v1")
  published_trimmed_dir <- file.path(graph_dir, "r_qhull_per_sample_trimmed_v1")
  official_dir <- file.path(graph_dir, "official_published")
  dir.create(untrimmed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(global_trimmed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(published_trimmed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(official_dir, recursive = TRUE, showWarnings = FALSE)
  coordinates_by_sample <- stats::setNames(vector("list", length(sample_ids)), sample_ids)
  graphs <- stats::setNames(vector("list", length(sample_ids)), sample_ids)
  per_sample <- vector("list", length(sample_ids))

  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[[i]]
    rows <- rows_by_sample[[sample_id]]
    coordinates <- as.matrix(reducedDim(sce, "spatial")[rows, , drop = FALSE])
    rownames(coordinates) <- colnames(sce)[rows]
    coordinates_by_sample[[sample_id]] <- coordinates
    path <- .sample_path(untrimmed_dir, sample_id)
    sample_started <- Sys.time()
    graph_created <- !file.exists(path)
    if (file.exists(path)) {
      graph <- readRDS(path)
    } else {
      # Qhull construction is strictly sample-local, matching Squidpy's
      # scipy.spatial.Delaunay algorithm without cross-section edges.
      graph <- delaunay_graph(coordinates)
      .atomic_save_rds(graph, path)
    }
    .validate_graph_identity(
      graph, colnames(sce)[rows], paste0("Untrimmed graph checkpoint `", sample_id, "`")
    )
    untrimmed_contract <- .graph_contract(
      sample_id, "untrimmed", colnames(sce)[rows], coordinates,
      upstream = list(
        representation_sha256 = cellcharter_file_sha256(sce_path),
        representation_contract_id = representation_contract$contract_id
      )
    )
    verification_candidate <- if (graph_created) graph else delaunay_graph(coordinates)
    .record_or_verify_graph(
      path, untrimmed_contract, verification_candidate,
      paste0("Untrimmed graph checkpoint `", sample_id, "`"), graph_created
    )
    graphs[[sample_id]] <- graph
    per_sample[[i]] <- data.frame(
      sample = sample_id, cells = length(rows), directed_edges_before = Matrix::nnzero(graph),
      delaunay_seconds = as.numeric(difftime(Sys.time(), sample_started, units = "secs")),
      untrimmed_path = path, stringsAsFactors = FALSE
    )
    message("Graph checkpoint ready: ", sample_id)
  }

  # v0.2.0 uses one percentile over the combined block-diagonal distance
  # matrix. This computes the same statistic from sample lists without ever
  # materialising a cohort-wide graph.
  cutoff <- cellcharter_global_distance_cutoff(graphs, coordinates_by_sample, 99)
  r_global_paths <- stats::setNames(character(length(sample_ids)), sample_ids)
  r_published_paths <- stats::setNames(character(length(sample_ids)), sample_ids)
  official_paths <- stats::setNames(character(length(sample_ids)), sample_ids)
  parity_rows <- vector("list", length(sample_ids))
  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[[i]]
    global_path <- .sample_path(global_trimmed_dir, sample_id)
    expected_global <- remove_long_links(
      graphs[[sample_id]], coordinates_by_sample[[sample_id]],
      percentile = 99, cutoff = cutoff
    )
    global_created <- !file.exists(global_path)
    if (global_created) {
      trimmed <- expected_global
      .atomic_save_rds(trimmed, global_path)
    }
    r_global_paths[[sample_id]] <- global_path
    r_global <- readRDS(global_path)
    .validate_graph_identity(
      r_global, colnames(sce)[rows_by_sample[[sample_id]]],
      paste0("Pinned-global graph checkpoint `", sample_id, "`")
    )
    global_contract <- .graph_contract(
      sample_id, "pinned_global_trimmed",
      colnames(sce)[rows_by_sample[[sample_id]]], coordinates_by_sample[[sample_id]],
      upstream = list(
        untrimmed_sha256 = cellcharter_file_sha256(
          .sample_path(untrimmed_dir, sample_id)
        ),
        representation_contract_id = representation_contract$contract_id
      ), cutoff = cutoff
    )
    .record_or_verify_graph(
      global_path, global_contract, expected_global,
      paste0("Pinned-global graph checkpoint `", sample_id, "`"), global_created
    )
    published_path <- .sample_path(published_trimmed_dir, sample_id)
    expected_published <- remove_long_links(
      graphs[[sample_id]], coordinates_by_sample[[sample_id]], percentile = 99
    )
    published_created <- !file.exists(published_path)
    if (published_created) {
      # The released H5AD empirically used a separate 99th-percentile cutoff
      # within each library. Reconstruct it in R and verify it edge-for-edge
      # against the published adjacency below.
      published_graph <- expected_published
      .atomic_save_rds(published_graph, published_path)
    } else {
      published_graph <- readRDS(published_path)
    }
    .validate_graph_identity(
      published_graph, colnames(sce)[rows_by_sample[[sample_id]]],
      paste0("Published-cutoff graph checkpoint `", sample_id, "`")
    )
    published_contract <- .graph_contract(
      sample_id, "published_per_sample_trimmed",
      colnames(sce)[rows_by_sample[[sample_id]]], coordinates_by_sample[[sample_id]],
      upstream = list(
        untrimmed_sha256 = cellcharter_file_sha256(
          .sample_path(untrimmed_dir, sample_id)
        ),
        representation_contract_id = representation_contract$contract_id
      ), cutoff = attr(expected_published, "distance_cutoff")
    )
    .record_or_verify_graph(
      published_path, published_contract, expected_published,
      paste0("Published-cutoff graph checkpoint `", sample_id, "`"), published_created
    )
    r_published_paths[[sample_id]] <- published_path
    official_path <- .sample_path(official_dir, sample_id)
    expected_official <- .read_h5ad_csr_block(
      clustered_h5ad, "obsp/spatial_connectivities", rows_by_sample[[sample_id]]
    )
    official_created <- !file.exists(official_path)
    if (official_created) {
      official <- expected_official
      .atomic_save_rds(official, official_path)
    } else {
      official <- readRDS(official_path)
    }
    .validate_graph_identity(
      official, colnames(sce)[rows_by_sample[[sample_id]]],
      paste0("Official graph checkpoint `", sample_id, "`")
    )
    official_contract <- .graph_contract(
      sample_id, "official_published",
      colnames(sce)[rows_by_sample[[sample_id]]], coordinates_by_sample[[sample_id]],
      upstream = list(
        clustered_h5ad_sha256 = cellcharter_file_sha256(clustered_h5ad),
        representation_contract_id = representation_contract$contract_id
      )
    )
    .record_or_verify_graph(
      official_path, official_contract, expected_official,
      paste0("Official graph checkpoint `", sample_id, "`"), official_created
    )
    official_paths[[sample_id]] <- official_path
    published_delta <- Matrix::drop0(published_graph - official)
    global_delta <- Matrix::drop0(r_global - official)
    parity_rows[[i]] <- data.frame(
      sample = sample_id,
      published_cutoff = attr(published_graph, "distance_cutoff"),
      r_published_directed_edges = Matrix::nnzero(published_graph),
      official_directed_edges = Matrix::nnzero(official),
      published_r_only_edges = sum(published_delta@x > 0),
      published_official_only_edges = sum(published_delta@x < 0),
      pinned_global_directed_edges = Matrix::nnzero(r_global),
      global_r_only_edges = sum(global_delta@x > 0),
      global_official_only_edges = sum(global_delta@x < 0),
      stringsAsFactors = FALSE
    )
  }
  if (identical(graph_source, "r_published")) {
    selected_paths <- r_published_paths
  } else if (identical(graph_source, "r_pinned_global")) {
    selected_paths <- r_global_paths
  } else {
    selected_paths <- official_paths
  }
  graph_parity <- do.call(rbind, parity_rows)
  if (identical(graph_source, "r_published") &&
      any(graph_parity$published_r_only_edges != 0L |
          graph_parity$published_official_only_edges != 0L)) {
    stop("R per-sample graph reconstruction does not match the published adjacency.")
  }
  manifest <- list(
    reference = "CellCharter 0.2.0 / Squidpy 1.3.0",
    r_qhull = "geometry::delaunayn with scipy-compatible Qbb Qc Qz Q12 Qt options",
    pinned_cutoff_scope = "global over directed positive distances in sample-local block graphs",
    published_cutoff_scope = "independently within each physical sample (empirically verified from the released H5AD)",
    percentile = 99, cutoff = cutoff,
    analysis_graph_source = if (identical(graph_source, "r_published")) {
      "R-native Qhull with per-sample published-artifact cutoffs"
    } else if (identical(graph_source, "r_pinned_global")) {
      "R-native Qhull with pinned CellCharter 0.2.0 global cutoff"
    } else {
      "official published H5AD adjacency (validation-only fixed input)"
    },
    per_sample = do.call(rbind, per_sample),
    graph_parity = graph_parity,
    cell_order_sha256 = .cc_observation_order_sha256(colnames(sce)),
    sample_cell_order_sha256 = vapply(
      rows_by_sample,
      function(rows) .cc_observation_order_sha256(colnames(sce)[rows]), character(1)
    ),
    r_published_paths = r_published_paths,
    r_pinned_global_paths = r_global_paths,
    official_paths = official_paths,
    trimmed_paths = selected_paths
  )
  .atomic_save_rds(rows_by_sample, rows_path)
  .atomic_save_rds(manifest, graph_manifest_path)
  graph_manifest_created <- TRUE
  .record_stage("graph", started, list(
    checkpoint = graph_manifest_path, distance_cutoff = cutoff,
    analysis_graph_source = manifest$analysis_graph_source,
    graph_parity = manifest$graph_parity,
    in_memory_graph_bytes = as.numeric(utils::object.size(graphs)),
    in_memory_coordinate_bytes = as.numeric(utils::object.size(coordinates_by_sample))
  ))
  graph_checkpoint_complete <- TRUE
}
if (!graph_checkpoint_complete) {
  stop("Graph checkpoint bundle is incomplete; run with CODEX_STAGES=graph (or the default workflow) to regenerate missing artifacts.")
}
if (!file.exists(graph_manifest_path)) stop("Graph checkpoint is absent: ", graph_manifest_path)
graph_manifest <- readRDS(graph_manifest_path)
expected_graph_source <- if (identical(graph_source, "official")) {
  "official published H5AD adjacency (validation-only fixed input)"
} else if (identical(graph_source, "r_published")) {
  "R-native Qhull with per-sample published-artifact cutoffs"
} else {
  "R-native Qhull with pinned CellCharter 0.2.0 global cutoff"
}
if (!identical(graph_manifest$analysis_graph_source, expected_graph_source)) {
  stop("Existing graph manifest uses a different CODEX_GRAPH_SOURCE; use a new CODEX_RESULTS_DIR.")
}
graphs <- lapply(graph_manifest$trimmed_paths, readRDS)
names(graphs) <- names(graph_manifest$trimmed_paths)
for (sample_id in sample_ids) {
  .validate_graph_identity(
    graphs[[sample_id]], colnames(sce)[rows_by_sample[[sample_id]]],
    paste0("Selected graph checkpoint `", sample_id, "`")
  )
}
if (!is.null(graph_manifest$cell_order_sha256) &&
    !identical(
      graph_manifest$cell_order_sha256,
      .cc_observation_order_sha256(colnames(sce))
    )) {
  stop("Graph manifest cell order does not match the representation SCE checkpoint.")
}
if (!file.exists(rows_path)) .atomic_save_rds(rows_by_sample, rows_path)
graph_bundle_artifacts <- c(
  manifest = graph_manifest_path,
  rows_by_sample = rows_path,
  stats::setNames(
    unname(graph_manifest$r_published_paths),
    paste0("published_", seq_along(graph_manifest$r_published_paths))
  ),
  stats::setNames(
    unname(graph_manifest$r_pinned_global_paths),
    paste0("pinned_global_", seq_along(graph_manifest$r_pinned_global_paths))
  ),
  stats::setNames(
    unname(graph_manifest$official_paths),
    paste0("official_", seq_along(graph_manifest$official_paths))
  )
)
graph_contract <- cellcharter_stage_contract(
  stage = "graph",
  implementation_digest = graph_implementation,
  upstream = list(
    representation_sha256 = cellcharter_file_sha256(sce_path),
    representation_contract_id = representation_contract$contract_id,
    clustered_h5ad_sha256 = cellcharter_file_sha256(clustered_h5ad)
  ),
  cell_order_sha256 = .cc_observation_order_sha256(colnames(sce)),
  parameters = list(
    samples = sample_ids, graph_source = graph_source,
    percentile = 99, cutoff_scope = if (identical(graph_source, "r_pinned_global")) {
      "global directed positive distances"
    } else if (identical(graph_source, "r_published")) {
      "per physical sample"
    } else {
      "official fixed input"
    },
    qhull_options = "Qbb Qc Qz Q12 Qt", cross_sample_edges = FALSE
  ),
  seeds = list(seed = NA_integer_),
  versions = .stage_versions(c("geometry", "Matrix", "rhdf5")), runtime = list()
)
if (graph_manifest_created) {
  cellcharter_record_checkpoint(
    graph_contract, graph_bundle_artifacts, status = "created_under_contract",
    evidence = list(
      sample_local_identity_validated = TRUE,
      published_graph_edge_parity = graph_manifest$graph_parity
    ), path = .contract_path("graph")
  )
} else if (!.require_contract("graph", graph_contract, graph_bundle_artifacts)) {
  stop(
    "Graph bundle has no matching creation contract. Pre-contract partial graph ",
    "files can be verified while completing CODEX_STAGES=graph, but a complete ",
    "uncontracted manifest is not accepted; use a new CODEX_RESULTS_DIR."
  )
}
if (!any(requested_stages %in% c(
  "aggregation", "autok", "tutorial_k11", "downstream", "report"
))) {
  message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
  quit(save = "no", status = 0L)
}
# CHECKPOINT CONTRACT REGION: graph END

# CHECKPOINT CONTRACT REGION: aggregation BEGIN
aggregation_dir <- file.path(
  result_dir,
  if (identical(graph_source, "r_pinned_global")) {
    "aggregation_pinned_global"
  } else {
    "aggregation"
  }
)
aggregation_manifest_path <- file.path(aggregation_dir, "aggregation_manifest.rds")
aggregation_implementation <- .stage_implementation("aggregation", list(
  aggregate_hdf5 = cellcharter_aggregate_hdf5,
  aggregate_rows = .cc_aggregate_rows_block,
  coordinates = .cc_coordinates,
  adjacency_topology = .cc_adjacency_topology,
  binary_graph = .cc_binary_graph,
  layers = .cc_layers,
  row_normalise = .cc_row_normalise,
  hop = .cc_hop,
  contiguous_runs = .cc_contiguous_runs,
  block_rows = .cc_aggregation_block_rows,
  block_indices = .cc_block_indices
))
.aggregation_sample_contract <- function(sample_id) {
  rows <- rows_by_sample[[sample_id]]
  graph_path <- graph_manifest$trimmed_paths[[sample_id]]
  cellcharter_stage_contract(
    stage = paste0("aggregation/", sample_id),
    implementation_digest = aggregation_implementation,
    upstream = list(
      representation_sha256 = cellcharter_file_sha256(sce_path),
      representation_contract_id = representation_contract$contract_id,
      graph_sha256 = cellcharter_file_sha256(graph_path),
      graph_contract_id = graph_contract$contract_id
    ),
    cell_order_sha256 = .cc_observation_order_sha256(colnames(sce)[rows]),
    parameters = list(
      sample_id = sample_id, dataset = "X_cellcharter", n_layers = 3L,
      hop_semantics = "exclusive", zero_neighbor = "zeros",
      block_rows = aggregation_block_rows, output_features = 40L
    ),
    seeds = list(seed = NA_integer_),
    versions = .stage_versions(c("HDF5Array", "DelayedArray", "Matrix")),
    runtime = list(storage = "HDF5", bounded_source_and_output_rows = TRUE)
  )
}
aggregation_contract <- cellcharter_stage_contract(
  stage = "aggregation",
  implementation_digest = aggregation_implementation,
  upstream = list(
    representation_sha256 = cellcharter_file_sha256(sce_path),
    representation_contract_id = representation_contract$contract_id,
    graph_manifest_sha256 = cellcharter_file_sha256(graph_manifest_path),
    graph_contract_id = graph_contract$contract_id
  ),
  cell_order_sha256 = .cc_observation_order_sha256(colnames(sce)),
  parameters = list(
    samples = sample_ids, dataset = "X_cellcharter", n_layers = 3L,
    hop_semantics = "exclusive", zero_neighbor = "zeros",
    block_rows = aggregation_block_rows, output_features = 40L,
    graph_source = graph_source, parity_rows_seed = 20230816L,
    parity_rows = 2048L, parity_tolerance = 1e-5
  ),
  seeds = list(parity_seed = 20230816L),
  versions = .stage_versions(c("HDF5Array", "DelayedArray", "Matrix")),
  runtime = list(storage = "HDF5", bounded_source_and_output_rows = TRUE)
)
aggregation_manifest_created <- FALSE
aggregation_checkpoint_complete <- file.exists(aggregation_manifest_path)
if (aggregation_checkpoint_complete) {
  existing_aggregation_manifest <- try(readRDS(aggregation_manifest_path), silent = TRUE)
  aggregation_checkpoint_complete <-
    !inherits(existing_aggregation_manifest, "try-error") &&
    is.character(existing_aggregation_manifest$sample_files) &&
    identical(names(existing_aggregation_manifest$sample_files), sample_ids) &&
    all(file.exists(existing_aggregation_manifest$sample_files)) &&
    all(file.exists(paste0(existing_aggregation_manifest$sample_files, ".complete.rds"))) &&
    file.exists(existing_aggregation_manifest$sce_checkpoint) &&
    file.exists(existing_aggregation_manifest$parity_checkpoint) &&
    is.character(existing_aggregation_manifest$published_graph_parity_checkpoint) &&
    length(existing_aggregation_manifest$published_graph_parity_checkpoint) == 1L &&
    file.exists(existing_aggregation_manifest$published_graph_parity_checkpoint) &&
    identical(
      existing_aggregation_manifest$cell_order_sha256,
      .cc_observation_order_sha256(colnames(sce))
    ) &&
    all(vapply(sample_ids, function(sample_id) {
      .aggregation_checkpoint_valid(
        existing_aggregation_manifest$sample_files[[sample_id]],
        colnames(sce)[rows_by_sample[[sample_id]]], sample_id,
        aggregation_block_rows, .aggregation_sample_contract(sample_id)
      )
    }, logical(1)))
}
if ("aggregation" %in% requested_stages && !aggregation_checkpoint_complete) {
  started <- Sys.time()
  dir.create(aggregation_dir, recursive = TRUE, showWarnings = FALSE)
  message(
    "CellCharterR aggregation: block_rows=", aggregation_block_rows,
    "; dense representation/output workspace is fixed by this row budget."
  )
  sample_files <- stats::setNames(character(length(sample_ids)), sample_ids)
  for (sample_id in sample_ids) {
    rows <- rows_by_sample[[sample_id]]
    path <- .sample_path(aggregation_dir, sample_id, ".h5")
    sample_files[[sample_id]] <- path
    complete_path <- paste0(path, ".complete.rds")
    expected_ids <- colnames(sce)[rows]
    valid_checkpoint <- .aggregation_checkpoint_valid(
      path, expected_ids, sample_id, aggregation_block_rows,
      .aggregation_sample_contract(sample_id)
    )
    if (!valid_checkpoint) {
      if (file.exists(path) || file.exists(complete_path) ||
          file.exists(paste0(path, ".contract.rds"))) {
        stop(
          "Aggregation sample checkpoint `", sample_id,
          "` exists without its matching contract. It was not overwritten; ",
          "use a new CODEX_RESULTS_DIR."
        )
      }
      x_sample <- reducedDim(sce, "X_trVAE")[rows, , drop = FALSE]
      temporary <- tempfile(
        paste0(".", basename(path), "-"), dirname(path), fileext = ".h5"
      )
      sample_r_heap_peak_before_mb <- .peak_r_memory_mb()
      sample_process_peak_before_mb <- .peak_process_memory_mb()
      aggregation_output <- cellcharter_aggregate_hdf5(
        x_sample, stats::setNames(list(graphs[[sample_id]]), sample_id),
        stats::setNames(list(seq_along(rows)), sample_id), temporary,
        name = "X_cellcharter", n_layers = 3L,
        block_rows = aggregation_block_rows
      )
      aggregation_provenance <- attr(aggregation_output, "aggregation_provenance")
      sample_r_heap_peak_after_mb <- .peak_r_memory_mb()
      sample_process_peak_after_mb <- .peak_process_memory_mb()
      produced <- HDF5Array::HDF5Array(temporary, "X_cellcharter")
      if (!identical(as.integer(dim(produced)), c(length(rows), 40L))) {
        stop("Incomplete aggregation output for sample ", sample_id, ".")
      }
      if (!file.rename(temporary, path)) {
        stop("Could not atomically publish aggregation checkpoint: ", path)
      }
      .atomic_save_rds(list(
        dimensions = dim(produced), sha256 = cellcharter_file_sha256(path),
        sample_id = sample_id, cell_ids = expected_ids,
        cell_order_sha256 = .cc_observation_order_sha256(expected_ids),
        block_rows = aggregation_block_rows,
        estimated_dense_workspace_bytes =
          aggregation_provenance$estimated_dense_workspace_bytes,
        estimated_dense_workspace_scope = aggregation_provenance$estimate_scope,
        measured_r_heap_peak_before_mb = sample_r_heap_peak_before_mb,
        measured_r_heap_peak_after_mb = sample_r_heap_peak_after_mb,
        measured_r_heap_peak_increment_mb = max(
          0, sample_r_heap_peak_after_mb - sample_r_heap_peak_before_mb
        ),
        measured_process_peak_before_mb = sample_process_peak_before_mb,
        measured_process_peak_after_mb = sample_process_peak_after_mb,
        measured_process_peak_increment_mb = if (
          is.na(sample_process_peak_before_mb) || is.na(sample_process_peak_after_mb)
        ) NA_real_ else max(
          0, sample_process_peak_after_mb - sample_process_peak_before_mb
        ),
        measured_process_peak_source = if (is.na(sample_process_peak_after_mb)) {
          "unavailable on this platform"
        } else {
          "Linux /proc/self/status VmHWM"
        },
        finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
      ), complete_path)
      cellcharter_record_checkpoint(
        .aggregation_sample_contract(sample_id),
        c(hdf5 = path, completion_marker = complete_path),
        status = "created_under_contract",
        evidence = list(
          cell_identity_validated = TRUE,
          dimensions = c(length(expected_ids), 40L),
          bounded_block_rows = aggregation_block_rows
        ), path = paste0(path, ".contract.rds")
      )
    }
    message("Aggregation checkpoint ready: ", sample_id)
  }
  if (!identical(unname(unlist(rows_by_sample)), seq_len(ncol(sce)))) {
    stop("CODEX rows are not sample-contiguous; aggregation checkpoint assembly needs an explicit reorder.")
  }
  matrices <- lapply(sample_files, HDF5Array::HDF5Array, name = "X_cellcharter")
  x_cellcharter <- do.call(base::rbind, matrices)
  rownames(x_cellcharter) <- colnames(sce)
  colnames(x_cellcharter) <- paste0("cellcharter_", seq_len(ncol(x_cellcharter)))
  reducedDim(sce, "X_cellcharter") <- x_cellcharter
  aggregated_sce_path <- file.path(aggregation_dir, "codex_sce_aggregated.rds")
  .atomic_save_rds(sce, aggregated_sce_path)

  # Compare a deterministic sample against the official float32 aggregate.
  official <- DelayedArray::t(HDF5Array::HDF5Array(clustered_h5ad, "obsm/X_cellcharter"))
  base::set.seed(20230816L)
  parity_rows <- sort(sample.int(nrow(x_cellcharter), 2048L))
  difference <- as.matrix(x_cellcharter[parity_rows, , drop = FALSE]) -
    as.matrix(official[parity_rows, , drop = FALSE])
  parity <- list(
    rows = parity_rows, max_absolute_error = max(abs(difference)),
    mean_absolute_error = mean(abs(difference)), tolerance = 1e-5,
    passed = max(abs(difference)) <= 1e-5,
    expected_to_match = !identical(graph_source, "r_pinned_global"),
    comparison_target = "published H5AD X_cellcharter",
    interpretation = if (identical(graph_source, "r_pinned_global")) {
      paste(
        "A mismatch is expected because the published artifact used per-sample",
        "cutoffs while pinned CellCharter 0.2.0 applies one global cutoff."
      )
    } else {
      "The selected graph reproduces the published artifact and must match."
    }
  )
  .atomic_save_rds(parity, file.path(aggregation_dir, "python_r_parity.rds"))
  if (isTRUE(parity$expected_to_match) && !isTRUE(parity$passed)) {
    stop(
      "R aggregation does not match the official X_cellcharter fixture (max error ",
      signif(parity$max_absolute_error, 6), "). AutoK was not started."
    )
  }
  published_graph_parity_path <- file.path(
    aggregation_dir, "published_graph_python_r_parity.rds"
  )
  if (identical(graph_source, "r_pinned_global")) {
    # Validate the released per-sample-cutoff aggregate independently while
    # realizing only one physical sample at a time. This keeps the default
    # source-faithful workflow self-contained without a second full HDF5 copy.
    published_differences <- vector("list", length(sample_ids))
    published_row_ids <- vector("list", length(sample_ids))
    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[[i]]
      sample_rows <- rows_by_sample[[sample_id]]
      selected_global_rows <- parity_rows[parity_rows %in% sample_rows]
      if (!length(selected_global_rows)) next
      selected_local_rows <- match(selected_global_rows, sample_rows)
      published_graph <- readRDS(graph_manifest$r_published_paths[[sample_id]])
      published_aggregate <- .cc_aggregate_rows_block(
        reducedDim(sce, "X_trVAE")[sample_rows, , drop = FALSE],
        .cc_adjacency_topology(published_graph, length(sample_rows)),
        selected_local_rows, .cc_layers(3L), aggregation_block_rows
      )
      published_differences[[i]] <-
        published_aggregate -
        as.matrix(official[selected_global_rows, , drop = FALSE])
      published_row_ids[[i]] <- selected_global_rows
    }
    published_difference <- do.call(base::rbind, published_differences)
    published_rows <- unlist(published_row_ids, use.names = FALSE)
    published_difference <- published_difference[order(published_rows), , drop = FALSE]
    published_rows <- sort(published_rows)
    published_parity <- list(
      rows = published_rows,
      max_absolute_error = max(abs(published_difference)),
      mean_absolute_error = mean(abs(published_difference)),
      tolerance = 1e-5,
      passed = max(abs(published_difference)) <= 1e-5,
      expected_to_match = TRUE,
      comparison_target = "published H5AD X_cellcharter",
      analysis_graph_source =
        "R-native Qhull with per-sample published-artifact cutoffs",
      interpretation = paste(
        "Bounded-memory per-sample validation of the released aggregation;",
        "no cohort-wide sensitivity matrix is materialised."
      )
    )
  } else {
    published_parity <- parity
  }
  .atomic_save_rds(published_parity, published_graph_parity_path)
  if (!isTRUE(published_parity$passed)) {
    stop(
      "R per-sample-cutoff aggregation does not match the published fixture (max error ",
      signif(published_parity$max_absolute_error, 6), ")."
    )
  }
  sample_sha256 <- vapply(sample_files, cellcharter_file_sha256, character(1))
  sample_aggregation_provenance <- lapply(
    sample_files, function(path) readRDS(paste0(path, ".complete.rds"))[
      c(
        "block_rows", "estimated_dense_workspace_bytes",
        "estimated_dense_workspace_scope", "measured_r_heap_peak_before_mb",
        "measured_r_heap_peak_after_mb", "measured_r_heap_peak_increment_mb",
        "measured_process_peak_before_mb", "measured_process_peak_after_mb",
        "measured_process_peak_increment_mb", "measured_process_peak_source"
      )
    ]
  )
  data_id <- aggregation_contract$contract_id
  manifest <- list(
    sample_files = sample_files, dataset = "X_cellcharter", n_layers = 3L,
    dimensions = dim(x_cellcharter), sce_checkpoint = aggregated_sce_path,
    parity_checkpoint = file.path(aggregation_dir, "python_r_parity.rds"),
    published_graph_parity_checkpoint = published_graph_parity_path,
    sample_sha256 = sample_sha256, data_id = data_id,
    block_rows = aggregation_block_rows,
    memory_provenance_by_sample = sample_aggregation_provenance,
    cell_order_sha256 = .cc_observation_order_sha256(colnames(sce)),
    sample_cell_order_sha256 = vapply(
      rows_by_sample,
      function(rows) .cc_observation_order_sha256(colnames(sce)[rows]), character(1)
    ),
    analysis_graph_source = graph_manifest$analysis_graph_source
  )
  .atomic_save_rds(manifest, aggregation_manifest_path)
  aggregation_manifest_created <- TRUE
  .record_stage("aggregation", started, list(
    checkpoint = aggregation_manifest_path, parity = parity,
    block_rows = aggregation_block_rows,
    estimated_dense_workspace_bytes = max(vapply(
      sample_aggregation_provenance,
      function(value) value$estimated_dense_workspace_bytes,
      numeric(1)
    )),
    measured_r_heap_peak_increment_mb = max(vapply(
      sample_aggregation_provenance,
      function(value) value$measured_r_heap_peak_increment_mb,
      numeric(1)
    )),
    measured_process_peak_increment_mb = {
      values <- vapply(
        sample_aggregation_provenance,
        function(value) value$measured_process_peak_increment_mb,
        numeric(1)
      )
      if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE)
    }
  ))
  aggregation_checkpoint_complete <- TRUE
}
if (!aggregation_checkpoint_complete) {
  stop(
    "Aggregation checkpoint bundle lacks immutable `cell_ids`/`cell_order_sha256` ",
    "metadata or is incomplete. Regenerate it with `CODEX_STAGES=aggregation` ",
    "under a new CODEX_RESULTS_DIR; legacy length-only checkpoints are not accepted."
  )
}
if (!file.exists(aggregation_manifest_path)) stop("Aggregation checkpoint is absent: ", aggregation_manifest_path)
aggregation_manifest <- readRDS(aggregation_manifest_path)
aggregation_bundle_artifacts <- c(
  manifest = aggregation_manifest_path,
  aggregated_sce = aggregation_manifest$sce_checkpoint,
  parity = aggregation_manifest$parity_checkpoint,
  published_graph_parity = aggregation_manifest$published_graph_parity_checkpoint,
  stats::setNames(
    unname(aggregation_manifest$sample_files),
    paste0("sample_", seq_along(aggregation_manifest$sample_files))
  )
)
if (aggregation_manifest_created) {
  cellcharter_record_checkpoint(
    aggregation_contract, aggregation_bundle_artifacts,
    status = "created_under_contract",
    evidence = list(
      cell_identity_validated = TRUE,
      parity_passed = readRDS(aggregation_manifest$published_graph_parity_checkpoint)$passed
    ), path = .contract_path("aggregation")
  )
} else if (!.require_contract(
  "aggregation", aggregation_contract, aggregation_bundle_artifacts
)) {
  stop(
    "Aggregation checkpoint contract is missing or stale. It was not accepted ",
    "by existence/dimensions; use a new CODEX_RESULTS_DIR."
  )
}
if (!identical(aggregation_manifest$block_rows, aggregation_block_rows)) {
  stop(
    "Aggregation manifest was produced with a different `block_rows` budget. ",
    "Regenerate with `CODEX_STAGES=aggregation` under a new CODEX_RESULTS_DIR."
  )
}
if (!identical(
      aggregation_manifest$cell_order_sha256,
      .cc_observation_order_sha256(colnames(sce))
    )) {
  stop(
    "Aggregation manifest lacks or mismatches `cell_order_sha256`. Regenerate ",
    "with `CODEX_STAGES=aggregation` under a new CODEX_RESULTS_DIR."
  )
}
if (!identical(
      aggregation_manifest$analysis_graph_source,
      graph_manifest$analysis_graph_source
    )) {
  stop(
    "Aggregation checkpoint was produced from a different graph source. ",
    "A content-addressed artifact is immutable; regenerate under a new ",
    "CODEX_RESULTS_DIR."
  )
}
representation_sce <- sce
sce <- readRDS(aggregation_manifest$sce_checkpoint)
.cc_validate_sce_identity(
  sce, "Aggregated SCE checkpoint", "X_cellcharter",
  aggregation_manifest$cell_order_sha256
)
invisible(.cc_align_observation_ids(
  colnames(sce), colnames(representation_sce),
  "Aggregated SCE checkpoint cell identifiers",
  "representation SCE checkpoint cell identifiers"
))
if ("reference_cluster_cellcharter" %in% colnames(colData(representation_sce))) {
  colData(sce)$reference_cluster_cellcharter <-
    colData(representation_sce)$reference_cluster_cellcharter
}
metadata(sce)$CellCharterR <- metadata(representation_sce)$CellCharterR
x_cellcharter <- reducedDim(sce, "X_cellcharter")
if (!any(requested_stages %in% c("autok", "tutorial_k11", "downstream", "report"))) {
  message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
  quit(save = "no", status = 0L)
}
# CHECKPOINT CONTRACT REGION: aggregation END

# CHECKPOINT CONTRACT REGION: autok BEGIN
autok_dir <- file.path(
  result_dir,
  if (identical(graph_source, "r_pinned_global")) "autok_pinned_global" else "autok"
)
autok_path <- file.path(autok_dir, "autok_result.rds")
selected_dir <- file.path(
  result_dir,
  if (identical(graph_source, "r_pinned_global")) "selected" else "selected_published"
)
labels_path <- file.path(selected_dir, "selected_labels.rds")
selected_model_path <- file.path(selected_dir, "selected_model.rds")
tutorial_k11_dir <- file.path(result_dir, "tutorial_k11_standalone")
tutorial_k11_model_path <- file.path(tutorial_k11_dir, "model.rds")
tutorial_k11_model_provenance_path <- file.path(
  tutorial_k11_dir, "model_provenance.rds"
)
tutorial_k11_labels_path <- file.path(tutorial_k11_dir, "labels.rds")
reference_dir <- file.path(result_dir, "reference")
reference_labels_path <- file.path(
  reference_dir, "published_h5ad_cluster_cellcharter_labels.rds"
)
if (!"torch" %in% loadedNamespaces()) {
  stop("Internal error: Torch was not preloaded for the requested model stage.")
}
cell_order_sha256 <- .cc_observation_order_sha256(colnames(sce))
gmm_implementation_functions <- list(
  fit = torch_gmm_fit,
  predict = torch_gmm_predict,
  predict_proba = torch_gmm_predict_proba,
  reader = .rc_reader,
  positive_integer = .rc_positive_integer,
  nonnegative_number = .rc_nonnegative_number,
  get = .rc_get,
  tensor = .rc_tensor,
  kmeans_distances = .rc_kmeans_distances,
  kmeans_tolerance = .rc_kmeans_tolerance,
  kmeans_plus_plus = .rc_kmeans_plus_plus,
  stream_kmeans = .rc_stream_kmeans,
  log_probability = .rc_logprob,
  cholesky_failure = .rc_cholesky_failure,
  cholesky = .rc_cholesky,
  estep = .rc_estep_stats,
  mstep_means_weights = .rc_mstep_means_weights,
  hard_initialization = .rc_hard_initialize,
  centered_covariances = .rc_centered_covariances,
  model_device = .rc_model_device,
  runtime_info = torch_gmm_runtime_info
)
autok_runtime <- torch_gmm_runtime_info("auto", batch_size, seed, num_threads)
autok_contract <- cellcharter_stage_contract(
  stage = "autok",
  implementation_digest = .stage_implementation("autok", c(list(
    autok = cellcharter_autok_fit,
    autok_k_grid = .autok_k_grid,
    autok_read_checkpoint = .autok_read_checkpoint,
    autok_mirror_stability = .autok_mirror_stability,
    autok_mape = .autok_mape,
    similarity = cellcharter_fowlkes_mallows
  ), gmm_implementation_functions)),
  upstream = list(
    aggregation_manifest_sha256 = cellcharter_file_sha256(aggregation_manifest_path),
    aggregation_data_id = aggregation_manifest$data_id,
    aggregation_contract_id = aggregation_contract$contract_id
  ),
  cell_order_sha256 = cell_order_sha256,
  parameters = list(
    candidate_bounds = autok_bounds, fitted_bounds = c(autok_bounds[[1L]] - 1L,
                                                       autok_bounds[[2L]] + 1L),
    max_runs = max_runs, convergence_tol = convergence_tol,
    batch_size = batch_size, max_iter = gmm_max_iter,
    tolerance = gmm_tolerance,
    covariance_regularization = covariance_regularization,
    covariance_type = "full", initialization = "hard nearest-centroid",
    prediction_batch_size = batch_size, device_request = "auto",
    num_threads = num_threads,
    model_selection = "lowest final NLL within K; highest FM stability across K"
  ),
  seeds = list(base_seed = seed),
  versions = .stage_versions(c("torch", "Matrix")), runtime = autok_runtime
)

role_metadata <- function(role_id, role, role_seed, model_path) {
  list(
    role_id = role_id,
    role = role,
    seed = role_seed,
    model_path = model_path,
    cell_ids = colnames(sce),
    cell_order_sha256 = cell_order_sha256,
    aggregation_data_id = aggregation_manifest$data_id
  )
}

# The released H5AD labels are a fixed comparison artifact. They are never an
# R fit and must not be described as the tutorial-matched R K=11 sensitivity.
if (!file.exists(reference_labels_path)) {
  if (!"reference_cluster_cellcharter" %in% colnames(colData(sce))) {
    stop("The published H5AD `cluster_cellcharter` labels are unavailable.")
  }
  published_reference <- c(role_metadata(
    role_id = "published_h5ad_reference",
    role = "Published H5AD cluster_cellcharter reference artifact (not an R fit)",
    role_seed = NA_integer_,
    model_path = NA_character_
  ), list(
    k = 11L,
    labels = as.character(colData(sce)$reference_cluster_cellcharter),
    source_path = clustered_h5ad,
    source_sha256 = cellcharter_file_sha256(clustered_h5ad),
    source_column = "cluster_cellcharter"
  ))
  .atomic_save_rds(published_reference, reference_labels_path)
}
published_reference <- readRDS(reference_labels_path)
if (!identical(published_reference$role_id, "published_h5ad_reference") ||
    !identical(published_reference$cell_ids, colnames(sce)) ||
    !identical(published_reference$cell_order_sha256, cell_order_sha256) ||
    !identical(published_reference$aggregation_data_id, aggregation_manifest$data_id) ||
    length(published_reference$labels) != ncol(sce)) {
  stop(
    "Published H5AD reference checkpoint has the wrong role or immutable identity; ",
    "use a new CODEX_RESULTS_DIR."
  )
}

if ("autok" %in% requested_stages) {
  started <- Sys.time()
  runtime <- autok_runtime
  fit_data_id <- autok_contract$contract_id
  autok <- cellcharter_autok_fit(
    x_cellcharter, n_clusters = autok_bounds, checkpoint_dir = autok_dir,
    max_runs = max_runs, convergence_tol = convergence_tol, seed = seed,
    fit_args = list(
      batch_size = batch_size, max_iter = gmm_max_iter,
      tolerance = gmm_tolerance,
      covariance_regularization = covariance_regularization,
      device = "auto", num_threads = num_threads
    ),
    predict_args = list(batch_size = batch_size, device = "auto"),
    data_id = fit_data_id, resume = TRUE
  )
  dir.create(selected_dir, recursive = TRUE, showWarnings = FALSE)
  best_row <- autok$best_index[autok$best_index$k == autok$best_k, , drop = FALSE]
  best_artifact <- readRDS(best_row$path[[1L]])
  if (!identical(best_artifact$data_id, autok$data_id) ||
      length(best_artifact$labels) != ncol(sce)) {
    stop("Selected AutoK fit does not match the immutable aggregation identity.")
  }
  .atomic_save_rds(best_artifact$model, selected_model_path)
  .atomic_save_rds(list(
    k = autok$best_k, labels = best_artifact$labels, model_path = selected_model_path,
    autok_fit_artifact_path = best_row$path[[1L]],
    selected_model_path = selected_model_path,
    cell_ids = colnames(sce),
    cell_order_sha256 = cell_order_sha256,
    aggregation_data_id = aggregation_manifest$data_id,
    aggregation_contract_id = aggregation_contract$contract_id,
    autok_contract_id = autok_contract$contract_id,
    autok_fit_data_id = autok$data_id,
    seed = as.integer(best_artifact$seed),
    autok_base_seed = seed,
    role_id = "r_autok_primary",
    role = "Primary CellCharterR result selected independently by AutoK"
  ), labels_path)
  cellcharter_record_checkpoint(
    autok_contract,
    c(result = autok_path, selected_model = selected_model_path,
      selected_labels = labels_path),
    status = "created_under_contract",
    evidence = list(
      selected_k = autok$best_k, runs_completed = autok$runs_completed,
      effective_per_fit_data_id = autok$data_id
    ), path = .contract_path("autok")
  )
  .record_stage("autok", started, list(
    checkpoint = autok_path, selected_k = autok$best_k, runtime = runtime,
    candidate_bounds = autok_bounds, fitted_grid = autok$fit_k
  ))
}
if (any(requested_stages %in% c("autok", "downstream", "report")) &&
    !.require_contract(
      "autok", autok_contract,
      c(result = autok_path, selected_model = selected_model_path,
        selected_labels = labels_path)
    )) {
  stop(
    "AutoK/selected-role checkpoint is absent or stale for the exact data, ",
    "implementation, parameters, seed, and runtime. Run CODEX_STAGES=autok ",
    "under a new result root if stale."
  )
}
# CHECKPOINT CONTRACT REGION: autok END

# Plotting is deliberately outside the fitted-model contract region: changing
# presentation must not invalidate an otherwise identical expensive AutoK fit.
# Persist both the exact plotted values and the rendered diagnostic so a report
# rebuild never needs to refit AutoK or infer stability from another statistic.
if (any(requested_stages %in% c("autok", "report"))) {
  autok_diagnostic_dir <- file.path(autok_dir, "diagnostics")
  autok_diagnostic_path <- file.path(
    autok_diagnostic_dir, "stability_k1_k20.rds"
  )
  autok_diagnostic_plot_path <- file.path(
    autok_diagnostic_dir, "stability_k1_k20.png"
  )
  if (!exists("autok", inherits = FALSE)) autok <- readRDS(autok_path)
  dir.create(autok_diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
  autok_diagnostic <- cellcharter_autok_stability_data(
    autok, k_range = 1:20, highlight_k = c(4L, 11L)
  )
  autok_diagnostic_plot <- cellcharter_plot_autok_stability(
    autok, k_range = 1:20, highlight_k = c(4L, 11L),
    title = "Mouse spleen CellCharterR AutoK stability"
  )
  .atomic_save_rds(autok_diagnostic, autok_diagnostic_path)
  utils::write.csv(
    autok_diagnostic$summary,
    file.path(autok_diagnostic_dir, "stability_k1_k20_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    autok_diagnostic$scores,
    file.path(autok_diagnostic_dir, "stability_k1_k20_scores.csv"),
    row.names = FALSE
  )
  temporary_plot <- tempfile(
    ".stability_k1_k20-", tmpdir = autok_diagnostic_dir, fileext = ".png"
  )
  ggplot2::ggsave(
    temporary_plot, plot = autok_diagnostic_plot,
    width = 10, height = 6, dpi = 180
  )
  if (!file.rename(temporary_plot, autok_diagnostic_plot_path)) {
    unlink(temporary_plot)
    stop("Could not atomically checkpoint the AutoK stability plot.")
  }
}

# CHECKPOINT CONTRACT REGION: tutorial_k11 BEGIN
tutorial_runtime <- torch_gmm_runtime_info(
  "auto", batch_size, tutorial_k11_seed, num_threads
)
tutorial_contract <- cellcharter_stage_contract(
  stage = "tutorial_k11",
  implementation_digest = .stage_implementation(
    "tutorial_k11", gmm_implementation_functions
  ),
  upstream = list(
    aggregation_manifest_sha256 = cellcharter_file_sha256(aggregation_manifest_path),
    aggregation_data_id = aggregation_manifest$data_id,
    aggregation_contract_id = aggregation_contract$contract_id
  ),
  cell_order_sha256 = cell_order_sha256,
  parameters = list(
    role_id = "r_tutorial_k11_sensitivity", n_clusters = 11L,
    batch_size = batch_size, max_iter = gmm_max_iter,
    tolerance = gmm_tolerance,
    covariance_regularization = covariance_regularization,
    covariance_type = "full", initialization = "hard nearest-centroid",
    prediction_batch_size = batch_size, device_request = "auto",
    num_threads = num_threads
  ),
  seeds = list(seed = tutorial_k11_seed),
  versions = .stage_versions(c("torch", "Matrix")), runtime = tutorial_runtime
)
if ("tutorial_k11" %in% requested_stages) {
  started <- Sys.time()
  dir.create(tutorial_k11_dir, recursive = TRUE, showWarnings = FALSE)
  runtime <- tutorial_runtime
  tutorial_fit_id <- tutorial_contract$contract_id
  tutorial_artifacts <- c(
    model = tutorial_k11_model_path,
    model_provenance = tutorial_k11_model_provenance_path,
    labels = tutorial_k11_labels_path
  )
  tutorial_reusable <- .require_contract(
    "tutorial_k11", tutorial_contract, tutorial_artifacts
  )
  if (!tutorial_reusable && any(file.exists(tutorial_artifacts))) {
    stop(
      "Standalone K=11 artifacts exist without the exact expected contract. ",
      "They were not overwritten; use a new CODEX_RESULTS_DIR."
    )
  }
  tutorial_model_created <- !tutorial_reusable
  if (tutorial_model_created) {
    tutorial_k11_model <- torch_gmm_fit(
      x_cellcharter,
      n_clusters = 11L,
      batch_size = batch_size,
      max_iter = gmm_max_iter,
      tolerance = gmm_tolerance,
      covariance_regularization = covariance_regularization,
      device = "auto",
      seed = tutorial_k11_seed,
      num_threads = num_threads
    )
    .atomic_save_rds(tutorial_k11_model, tutorial_k11_model_path)
    .atomic_save_rds(c(role_metadata(
      role_id = "r_tutorial_k11_sensitivity",
      role = "Standalone CellCharterR tutorial-matched K=11 model checkpoint",
      role_seed = tutorial_k11_seed,
      model_path = tutorial_k11_model_path
    ), list(
      k = 11L,
      fit_id = tutorial_fit_id,
      fit_function = "torch_gmm_fit",
      runtime = runtime,
      tutorial_contract_id = tutorial_contract$contract_id,
      aggregation_contract_id = aggregation_contract$contract_id
    )), tutorial_k11_model_provenance_path)
  } else {
    if (!file.exists(tutorial_k11_model_provenance_path)) {
      stop(
        "Standalone tutorial K=11 model lacks its creation-time provenance; ",
        "use a new CODEX_RESULTS_DIR."
      )
    }
    tutorial_k11_model <- readRDS(tutorial_k11_model_path)
  }
  tutorial_model_provenance <- readRDS(tutorial_k11_model_provenance_path)
  if (!identical(
        tutorial_model_provenance$role_id, "r_tutorial_k11_sensitivity"
      ) || !identical(tutorial_model_provenance$k, 11L) ||
      !identical(tutorial_model_provenance$seed, tutorial_k11_seed) ||
      !identical(tutorial_model_provenance$fit_id, tutorial_fit_id) ||
      !identical(tutorial_model_provenance$cell_ids, colnames(sce)) ||
      !identical(tutorial_model_provenance$cell_order_sha256, cell_order_sha256) ||
      !identical(
        tutorial_model_provenance$aggregation_data_id,
        aggregation_manifest$data_id
      )) {
    stop(
      "Standalone tutorial K=11 model provenance has the wrong role or immutable ",
      "identity; use a new CODEX_RESULTS_DIR."
    )
  }
  if (!inherits(tutorial_k11_model, "torch_gmm") ||
      !identical(tutorial_k11_model$n_clusters, 11L) ||
      !identical(tutorial_k11_model$runtime$seed, tutorial_k11_seed)) {
    stop(
      "Standalone tutorial K=11 model checkpoint has the wrong role, K, or seed; ",
      "use a new CODEX_RESULTS_DIR."
    )
  }
  if (!tutorial_reusable) {
    tutorial_k11_labels <- torch_gmm_predict(
      tutorial_k11_model, x_cellcharter,
      batch_size = batch_size, device = "auto"
    )
    tutorial_checkpoint <- c(role_metadata(
      role_id = "r_tutorial_k11_sensitivity",
      role = paste(
        "Standalone CellCharterR tutorial-matched K=11 sensitivity fit;",
        "not an AutoK run and not the published H5AD labels"
      ),
      role_seed = tutorial_k11_seed,
      model_path = tutorial_k11_model_path
    ), list(
      k = 11L,
      labels = tutorial_k11_labels,
      fit_id = tutorial_fit_id,
      tutorial_contract_id = tutorial_contract$contract_id,
      aggregation_contract_id = aggregation_contract$contract_id,
      fit_function = "torch_gmm_fit",
      reference_call = "Cluster(n_clusters=11, random_state=12345)",
      reference_version = "CellCharter 0.2.0",
      reference_commit = "2a8fdde8213e32e5a00f1f9571bffd3dbb825a79"
    ))
    .atomic_save_rds(tutorial_checkpoint, tutorial_k11_labels_path)
    cellcharter_record_checkpoint(
      tutorial_contract, tutorial_artifacts, status = "created_under_contract",
      evidence = list(
        role_id = "r_tutorial_k11_sensitivity", k = 11L,
        effective_fit_id = tutorial_fit_id
      ), path = .contract_path("tutorial_k11")
    )
  }
  tutorial_checkpoint <- readRDS(tutorial_k11_labels_path)
  if (!identical(tutorial_checkpoint$role_id, "r_tutorial_k11_sensitivity") ||
      !identical(tutorial_checkpoint$k, 11L) ||
      !identical(tutorial_checkpoint$seed, tutorial_k11_seed) ||
      !identical(tutorial_checkpoint$tutorial_contract_id,
                 tutorial_contract$contract_id) ||
      !identical(tutorial_checkpoint$model_path, tutorial_k11_model_path) ||
      !identical(tutorial_checkpoint$cell_ids, colnames(sce)) ||
      !identical(tutorial_checkpoint$cell_order_sha256, cell_order_sha256) ||
      !identical(
        tutorial_checkpoint$aggregation_data_id, aggregation_manifest$data_id
      ) || length(tutorial_checkpoint$labels) != ncol(sce)) {
    stop(
      "Standalone tutorial K=11 label checkpoint has the wrong role or immutable ",
      "identity; use a new CODEX_RESULTS_DIR."
    )
  }
  .record_stage("tutorial_k11", started, list(
    checkpoint = tutorial_k11_labels_path,
    model_path = tutorial_k11_model_path,
    model_provenance_path = tutorial_k11_model_provenance_path,
    fit_id = tutorial_fit_id,
    role_id = "r_tutorial_k11_sensitivity",
    k = 11L,
    seed = tutorial_k11_seed,
    runtime = runtime,
    nll = tutorial_k11_model$nll
  ))
}
if (any(requested_stages %in% c("tutorial_k11", "report")) &&
    !.require_contract(
      "tutorial_k11", tutorial_contract,
      c(model = tutorial_k11_model_path,
        model_provenance = tutorial_k11_model_provenance_path,
        labels = tutorial_k11_labels_path)
    )) {
  stop("Standalone tutorial K=11 checkpoint contract is absent or stale.")
}
# CHECKPOINT CONTRACT REGION: tutorial_k11 END

# CHECKPOINT CONTRACT REGION: downstream BEGIN
downstream_dir <- file.path(
  result_dir,
  if (identical(graph_source, "r_pinned_global")) "downstream" else "downstream_published"
)
downstream_manifest_path <- file.path(downstream_dir, "downstream_manifest.rds")
downstream_needed <- any(requested_stages %in% c("downstream", "report"))
downstream_contract <- NULL
if (downstream_needed) {
  downstream_contract <- cellcharter_stage_contract(
    stage = "downstream",
    implementation_digest = .stage_implementation("downstream", list(
      enrichment = cellcharter_enrichment,
      neighborhood = cellcharter_nhood_enrichment,
      differential = cellcharter_diff_nhood_enrichment,
      complete_factor = .cc_complete_factor,
      neighborhood_inputs = .cc_nhood_inputs,
      neighborhood_from_entries = .cc_nhood_from_entries,
      neighborhood_from_counts = .cc_nhood_from_counts,
      neighborhood_library_stats = .cc_nhood_library_stats,
      neighborhood_from_library_stats = .cc_nhood_from_library_stats,
      differential_display = .cc_diff_nhood_display
    )),
    upstream = list(
      selected_labels_sha256 = cellcharter_file_sha256(labels_path),
      autok_contract_id = autok_contract$contract_id,
      graph_manifest_sha256 = cellcharter_file_sha256(graph_manifest_path),
      graph_contract_id = graph_contract$contract_id,
      aggregation_data_id = aggregation_manifest$data_id,
      aggregation_contract_id = aggregation_contract$contract_id,
      tutorial_k11_labels_sha256 = if (file.exists(tutorial_k11_labels_path)) {
        cellcharter_file_sha256(tutorial_k11_labels_path)
      } else {
        "absent"
      },
      published_reference_labels_sha256 = cellcharter_file_sha256(reference_labels_path)
    ),
    cell_order_sha256 = cell_order_sha256,
    parameters = list(
      cell_type_mapping = "CODEX reduced mapping v1",
      enrichment = "observed/analytical expected log2 fold change",
      neighborhood = "asymmetric analytical expectation by BALBc/MRL",
      differential_contrast = "BALBc-MRL",
      condition_groups = c("BALBc", "MRL"),
      library_unit = "sample", n_perms = n_perms,
      permutation_mode = "cellcharter_v0.2.0",
      tail = "directional inclusive without +1 correction",
      graph_source = graph_source,
      max_graph_bytes = as.numeric(Sys.getenv(
        "CODEX_MAX_GRAPH_BYTES", as.character(2 * 1024^3)
      ))
    ),
    seeds = list(differential_seed = seed),
    versions = .stage_versions(c("Matrix", "SingleCellExperiment")),
    runtime = list(library_order = "pandas first occurrence")
  )
}
downstream_manifest_created <- FALSE
downstream_checkpoint_complete <- file.exists(downstream_manifest_path)
if (downstream_checkpoint_complete) {
  existing_downstream_manifest <- try(readRDS(downstream_manifest_path), silent = TRUE)
  downstream_required <- if (inherits(existing_downstream_manifest, "try-error")) {
    character()
  } else {
    c(
      existing_downstream_manifest$enrichment,
      existing_downstream_manifest$neighborhood_enrichment,
      existing_downstream_manifest$differential,
      existing_downstream_manifest$report_input,
      existing_downstream_manifest$k11_reference_agreement,
      file.path(downstream_dir, "codex_sce_with_domains.rds")
    )
  }
  downstream_checkpoint_complete <-
    !inherits(existing_downstream_manifest, "try-error") &&
    length(downstream_required) == 6L &&
    all(nzchar(downstream_required)) && all(file.exists(downstream_required))
  if (downstream_checkpoint_complete && downstream_needed) {
    downstream_artifacts <- stats::setNames(
      downstream_required,
      c("enrichment", "neighborhood", "differential", "report_input",
        "k11_agreement", "sce")
    )
    downstream_artifacts <- c(manifest = downstream_manifest_path, downstream_artifacts)
    downstream_checkpoint_complete <- .require_contract(
      "downstream", downstream_contract, downstream_artifacts
    )
    if (!downstream_checkpoint_complete) {
      stop(
        "Downstream outputs exist without the exact expected contract. They ",
        "were not overwritten; use a new CODEX_RESULTS_DIR."
      )
    }
  }
}
if ("downstream" %in% requested_stages && !downstream_checkpoint_complete) {
  if (!file.exists(labels_path)) stop("Selected labels are absent; complete AutoK first.")
  started <- Sys.time()
  dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- readRDS(labels_path)
  if (is.null(selected$cell_ids) || is.null(selected$cell_order_sha256) ||
      is.null(selected$aggregation_data_id)) {
    stop(
      "Selected-label checkpoint lacks `cell_ids`, `cell_order_sha256`, or ",
      "`aggregation_data_id`. Regenerate with `CODEX_STAGES=autok` under a ",
      "new CODEX_RESULTS_DIR."
    )
  }
  selected_ids <- .cc_validate_observation_ids(
    selected$cell_ids, "Selected-label checkpoint cell identifiers"
  )
  .cc_align_observation_ids(
    selected_ids, colnames(sce), "Selected-label checkpoint cell identifiers",
    "aggregated SCE checkpoint cell identifiers"
  )
  if (!identical(
        selected$cell_order_sha256,
        .cc_observation_order_sha256(colnames(sce))
      ) || !identical(selected$aggregation_data_id, aggregation_manifest$data_id)) {
    stop(
      "Selected-label checkpoint lacks or mismatches immutable cell-order/data ",
      "identity. Regenerate with `CODEX_STAGES=autok` under a new CODEX_RESULTS_DIR."
    )
  }
  if (length(selected$labels) != ncol(sce)) {
    stop("Selected-label checkpoint has a label-count mismatch after identity validation.")
  }
  labels <- factor(selected$labels, levels = seq.int(0L, selected$k - 1L))
  colData(sce)$r_domain <- labels

  k11_agreement_path <- file.path(
    downstream_dir, "tutorial_sensitivity_k11_reference_agreement.rds"
  )
  if (file.exists(tutorial_k11_labels_path) &&
      "reference_cluster_cellcharter" %in% colnames(colData(sce))) {
    k11_selected <- readRDS(tutorial_k11_labels_path)
    if (is.null(k11_selected$cell_ids) ||
        is.null(k11_selected$cell_order_sha256) ||
        is.null(k11_selected$aggregation_data_id)) {
      stop(
        "K=11 label checkpoint lacks `cell_ids`, `cell_order_sha256`, or ",
        "`aggregation_data_id`. Regenerate with `CODEX_STAGES=tutorial_k11` under a ",
        "new CODEX_RESULTS_DIR."
      )
    }
    k11_ids <- .cc_validate_observation_ids(
      k11_selected$cell_ids, "K=11 label checkpoint cell identifiers"
    )
    .cc_align_observation_ids(
      k11_ids, colnames(sce), "K=11 label checkpoint cell identifiers",
      "aggregated SCE checkpoint cell identifiers"
    )
    if (!identical(
          k11_selected$cell_order_sha256,
          .cc_observation_order_sha256(colnames(sce))
        ) || !identical(k11_selected$aggregation_data_id, aggregation_manifest$data_id) ||
        length(k11_selected$labels) != ncol(sce)) {
      stop(
        "K=11 label checkpoint lacks or mismatches immutable cell-order/data ",
        "identity. Regenerate with `CODEX_STAGES=tutorial_k11` under a new CODEX_RESULTS_DIR."
      )
    }
    agreement <- list(
      metric = "Fowlkes-Mallows (label-permutation invariant)",
      value = cellcharter_fowlkes_mallows(
        k11_selected$labels, colData(sce)$reference_cluster_cellcharter
      ),
      r_k = 11L, reference_k = nlevels(droplevels(
        colData(sce)$reference_cluster_cellcharter
      )),
      reference = "published cluster_cellcharter column in official H5AD"
    )
    .atomic_save_rds(agreement, k11_agreement_path)
  }

  cell_type_reduced <- as.character(colData(sce)$cell_type)
  reduced_mapping <- c(
    "CD4(+)CD8(-)cDC" = "DCs", "CD4(-)CD8(+)cDC" = "DCs",
    "CD4(-)CD8(-) cDC" = "DCs",
    "CD106(+)CD16/32(+)CD31(+) stroma" = "stroma",
    "CD106(+)CD16/32(+)CD31(-)Ly6C(-) stroma" = "stroma",
    "CD106(-)CD16/32(-)Ly6C(+)CD31(+) stroma" = "stroma",
    "ERTR7(+) stroma" = "stroma", "CD3(+) other markers (-)" = "other",
    "CD106(+)CD16/32(-)Ly6C(+)CD31(+)" = "other",
    "CD106(-)CD16/32(+)Ly6C(+)CD31(-)" = "other",
    "noid" = "other", "CD4(+)MHCII(+)" = "other", "capsule" = "other"
  )
  mapped <- unname(reduced_mapping[cell_type_reduced])
  cell_type_reduced[!is.na(mapped)] <- mapped[!is.na(mapped)]
  # pandas.Categorical infers sorted categories for these string annotations;
  # preserve the same matrix row order in the tutorial enrichment output.
  cell_type_reduced <- factor(
    cell_type_reduced, levels = sort(unique(cell_type_reduced))
  )
  colData(sce)$cell_type_reduced <- cell_type_reduced
  enrichment <- cellcharter_enrichment(cell_type_reduced, labels)
  enrichment_path <- file.path(downstream_dir, "cell_type_enrichment.rds")
  .atomic_save_rds(enrichment, enrichment_path)

  if (!identical(unname(unlist(rows_by_sample)), seq_len(ncol(sce)))) {
    stop("Cannot safely assemble block-diagonal graph in SCE order.")
  }
  max_graph_bytes <- as.numeric(Sys.getenv("CODEX_MAX_GRAPH_BYTES", as.character(2 * 1024^3)))
  estimated_graph_bytes <- sum(vapply(
    graphs, function(x) as.numeric(utils::object.size(x)), numeric(1)
  ))
  if (estimated_graph_bytes > max_graph_bytes) {
    stop("Estimated combined graph size exceeds CODEX_MAX_GRAPH_BYTES: ", estimated_graph_bytes)
  }
  combined_graph <- do.call(Matrix::bdiag, unname(graphs))
  graph_bytes <- as.numeric(utils::object.size(combined_graph))
  condition <- droplevels(colData(sce)$condition)
  nhood <- lapply(levels(condition), function(group) {
    keep <- which(condition == group)
    cellcharter_nhood_enrichment(
      combined_graph[keep, keep, drop = FALSE], droplevels(labels[keep])
    )
  })
  names(nhood) <- levels(condition)
  nhood_path <- file.path(downstream_dir, "condition_neighborhood_enrichment.rds")
  .atomic_save_rds(nhood, nhood_path)

  differential <- cellcharter_diff_nhood_enrichment(
    combined_graph, labels, condition, colData(sce)$sample,
    condition_groups = c("BALBc", "MRL"), n_perms = n_perms,
    seed = seed, permutation_mode = "cellcharter_v0.2.0"
  )
  differential_path <- file.path(downstream_dir, "BALBc_vs_MRL_differential_neighborhood.rds")
  .atomic_save_rds(differential, differential_path)

  coordinates <- as.matrix(reducedDim(sce, "spatial"))
  report_input <- data.frame(
    cell_id = colnames(sce), sample = colData(sce)$sample,
    condition = condition, cell_type = colData(sce)$cell_type,
    cell_type_reduced = cell_type_reduced,
    x = coordinates[, 1L], y = coordinates[, 2L], r_domain = labels
  )
  attr(report_input, "artifact_metadata") <- list(
    role_id = selected$role_id,
    role = selected$role,
    seed = selected$seed,
    model_path = selected$model_path,
    cell_id_column = "cell_id",
    cell_order_sha256 = selected$cell_order_sha256,
    aggregation_data_id = selected$aggregation_data_id
  )
  report_input_path <- file.path(downstream_dir, "report_input.rds")
  .atomic_save_rds(report_input, report_input_path)
  .atomic_save_rds(sce, file.path(downstream_dir, "codex_sce_with_domains.rds"))
  manifest <- list(
    selected_k = selected$k, enrichment = enrichment_path,
    neighborhood_enrichment = nhood_path, differential = differential_path,
    report_input = report_input_path, n_perms = n_perms,
    k11_reference_agreement = k11_agreement_path,
    graph_bytes_estimated_before_allocation = estimated_graph_bytes,
    graph_bytes_measured = graph_bytes
  )
  .atomic_save_rds(manifest, downstream_manifest_path)
  downstream_manifest_created <- TRUE
  cellcharter_record_checkpoint(
    downstream_contract,
    c(
      manifest = downstream_manifest_path, enrichment = enrichment_path,
      neighborhood = nhood_path, differential = differential_path,
      report_input = report_input_path, k11_agreement = k11_agreement_path,
      sce = file.path(downstream_dir, "codex_sce_with_domains.rds")
    ),
    status = "created_under_contract",
    evidence = list(
      selected_k = selected$k, library_unit = "sample", n_perms = n_perms,
      differential_seed = seed
    ), path = .contract_path("downstream")
  )
  .record_stage("downstream", started, list(
    checkpoint = downstream_manifest_path, selected_k = selected$k,
    graph_bytes_measured = graph_bytes, n_perms = n_perms
  ))
  downstream_checkpoint_complete <- TRUE
}
# CHECKPOINT CONTRACT REGION: downstream END

# CHECKPOINT CONTRACT REGION: report BEGIN
if ("report" %in% requested_stages) {
  if (!downstream_checkpoint_complete && !file.exists(downstream_manifest_path)) {
    stop("Downstream checkpoints are absent; run with CODEX_STAGES=downstream,report.")
  }
  if (!downstream_checkpoint_complete) {
    stop("Downstream manifest is incomplete; run with CODEX_STAGES=downstream,report to regenerate missing artifacts.")
  }
  if (!identical(graph_source, "r_pinned_global")) {
    stop(
      "The primary/tutorial second-review report requires CODEX_GRAPH_SOURCE=",
      "r_pinned_global. Keep alternate graph-source reports separate."
    )
  }
  if (!file.exists(tutorial_k11_labels_path) || !file.exists(reference_labels_path)) {
    stop(
      "The report requires the standalone tutorial K=11 and published-reference ",
      "roles; run CODEX_STAGES=tutorial_k11,report after AutoK/downstream."
    )
  }
  deliverable_dir <- file.path(result_dir, "second_review_deliverables")
  report_html_path <- file.path(
    result_dir, "report", "CellCharterR_mouse_spleen_AutoK_K11.html"
  )
  report_inputs_path <- file.path(
    deliverable_dir, "CellCharterR_primary_K11_manifest.rds"
  )
  report_implementation <- digest::digest(list(
    workflow_region = .stage_implementation("report", list()),
    builder_source_sha256 = cellcharter_file_sha256(
      "scripts/build_second_review_deliverables.R"
    ),
    report_source_sha256 = cellcharter_file_sha256(
      "reports/CellCharterR_mouse_spleen_AutoK_K11.Rmd"
    )
  ), algo = "sha256")
  report_contract <- cellcharter_stage_contract(
    stage = "report",
    implementation_digest = report_implementation,
    upstream = list(
      downstream_manifest_sha256 = cellcharter_file_sha256(downstream_manifest_path),
      downstream_contract_id = downstream_contract$contract_id,
      selected_labels_sha256 = cellcharter_file_sha256(labels_path),
      autok_contract_id = autok_contract$contract_id,
      tutorial_k11_labels_sha256 = cellcharter_file_sha256(tutorial_k11_labels_path),
      tutorial_k11_contract_id = tutorial_contract$contract_id,
      published_reference_labels_sha256 = cellcharter_file_sha256(reference_labels_path)
    ),
    cell_order_sha256 = cell_order_sha256,
    parameters = list(
      output_file = basename(report_html_path),
      displayed_samples = c("BALBc-1", "MRL-4"),
      dot_size_threshold = 0.58, dot_color_threshold = 0.58,
      significance = 0.05, differential_directions = c("MRL-BALBc", "BALBc-MRL"),
      code_chunks_visible = TRUE
    ),
    seeds = list(downstream_seed = seed, tutorial_k11_seed = tutorial_k11_seed),
    versions = .stage_versions(c("rmarkdown", "knitr", "ggplot2")),
    runtime = list(
      pandoc = as.character(rmarkdown::pandoc_version()),
      renderer = "rmarkdown::render"
    )
  )
  report_artifacts <- c(html = report_html_path, report_inputs = report_inputs_path)
  report_reusable <- .require_contract("report", report_contract, report_artifacts)
  if (!report_reusable &&
      (any(file.exists(report_artifacts)) || file.exists(.contract_path("report")))) {
    stop(
      "Report/render artifacts exist without the exact expected contract. They ",
      "were not overwritten; remove them explicitly or use a new CODEX_RESULTS_DIR."
    )
  }
  if (!report_reusable) {
  started <- Sys.time()
  output_dir <- file.path(result_dir, "report")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  previous_workflow_dir <- Sys.getenv(
    "CELLCHARTERR_WORKFLOW_RESULTS", unset = NA_character_
  )
  previous_report_dir <- Sys.getenv(
    "CELLCHARTERR_REPORT_RESULTS", unset = NA_character_
  )
  tryCatch({
    Sys.setenv(
      CELLCHARTERR_WORKFLOW_RESULTS = result_dir,
      CELLCHARTERR_REPORT_RESULTS = deliverable_dir
    )
    sys.source(
      "scripts/build_second_review_deliverables.R",
      envir = new.env(parent = globalenv())
    )
  }, finally = {
    if (is.na(previous_workflow_dir)) {
      Sys.unsetenv("CELLCHARTERR_WORKFLOW_RESULTS")
    } else {
      Sys.setenv(CELLCHARTERR_WORKFLOW_RESULTS = previous_workflow_dir)
    }
    if (is.na(previous_report_dir)) {
      Sys.unsetenv("CELLCHARTERR_REPORT_RESULTS")
    } else {
      Sys.setenv(CELLCHARTERR_REPORT_RESULTS = previous_report_dir)
    }
  })
  rmarkdown::render(
    "reports/CellCharterR_mouse_spleen_AutoK_K11.Rmd",
    output_file = basename(report_html_path),
    output_dir = output_dir,
    params = list(result_dir = deliverable_dir),
    envir = new.env(parent = globalenv()),
    quiet = FALSE
  )
  .record_stage("report", started, list(
    html = report_html_path,
    report_inputs = file.path(
      deliverable_dir, "CellCharterR_primary_K11_manifest.rds"
    )
  ))
  cellcharter_record_checkpoint(
    report_contract, report_artifacts, status = "created_under_contract",
    evidence = list(
      render_completed = TRUE, code_chunks_visible = TRUE,
      report_inputs_contract_id = downstream_contract$contract_id
    ), path = .contract_path("report")
  )
  } else {
    message("Report checkpoint matches the expected contract; render skipped.")
  }
}
# CHECKPOINT CONTRACT REGION: report END

message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
