#!/usr/bin/env Rscript

# Checkpointed R reproduction of the CellCharter CODEX mouse-spleen workflow.
# Large matrices stay HDF5-backed and the analysis state is orchestrated by a
# SingleCellExperiment. Set CODEX_STAGES to a comma-separated subset to resume
# particular stages (representation,graph,aggregation,autok,downstream,report).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(HDF5Array)
  library(Matrix)
})
source("R/sce_workflow.R")
source("R/aggregation.R")
source("R/torch_gmm.R")
source("R/autok.R")
source("R/enrichment.R")

input_dir <- normalizePath(
  Sys.getenv("CODEX_INPUT_DIR", file.path("results", "codex_mouse_spleen", "input")),
  winslash = "/", mustWork = FALSE
)
result_dir <- normalizePath(
  Sys.getenv("CODEX_RESULTS_DIR", file.path("results", "codex_mouse_spleen")),
  winslash = "/", mustWork = FALSE
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

requested_stages <- strsplit(
  Sys.getenv("CODEX_STAGES", "representation,graph,aggregation,autok,downstream,report"),
  ",", fixed = TRUE
)[[1L]]
requested_stages <- trimws(requested_stages)
known_stages <- c("representation", "graph", "aggregation", "autok", "downstream", "report")
if (any(!requested_stages %in% known_stages)) {
  stop("Unknown CODEX_STAGES value(s): ",
       paste(setdiff(requested_stages, known_stages), collapse = ", "))
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
batch_size <- as.integer(Sys.getenv("CODEX_BATCH_SIZE", "8192"))
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

representation_dir <- file.path(result_dir, "representation")
sce_path <- file.path(representation_dir, "codex_sce.rds")
if ("representation" %in% requested_stages && !file.exists(sce_path)) {
  started <- Sys.time()
  sce <- read_codex_sce(raw_h5ad, clustered_h5ad)
  metadata(sce)$rCellCharter$parameters <- list(
    seed = seed, autok_bounds = autok_bounds, max_runs = max_runs,
    convergence_tol = convergence_tol, batch_size = batch_size,
    num_threads = num_threads, n_perms = n_perms,
    gmm_max_iter = gmm_max_iter, gmm_tolerance = gmm_tolerance,
    covariance_regularization = covariance_regularization,
    graph_source = graph_source
  )
  .atomic_save_rds(sce, sce_path)
  .record_stage("representation", started, list(
    checkpoint = sce_path,
    cells = ncol(sce), markers = nrow(sce), trvae_dimensions = ncol(reducedDim(sce, "X_trVAE"))
  ))
}
if (!file.exists(sce_path)) stop("Representation checkpoint is absent: ", sce_path)
sce <- readRDS(sce_path)
# Older representation checkpoints are upgraded only when fields are absent.
# Never replace creation-time parameters with values from a later resume or
# report-only invocation; stage-specific provenance records the latter.
requested_parameters <- list(
  seed = seed, autok_bounds = autok_bounds, max_runs = max_runs,
  convergence_tol = convergence_tol, batch_size = batch_size,
  num_threads = num_threads, n_perms = n_perms,
  gmm_max_iter = gmm_max_iter, gmm_tolerance = gmm_tolerance,
  covariance_regularization = covariance_regularization,
  graph_source = graph_source
)
reference_updated <- is.null(metadata(sce)$rCellCharter$reference$scipy_version)
if (reference_updated) {
  metadata(sce)$rCellCharter$reference$scipy_version <-
    "not declared by cellcharter_analyses environment"
}
autok_provenance_updated <- file.exists(official_autok_zip) &&
  is.null(metadata(sce)$rCellCharter$inputs$official_autok_sha256)
if (autok_provenance_updated) {
  metadata(sce)$rCellCharter$inputs$official_autok_zip <-
    normalizePath(official_autok_zip)
  metadata(sce)$rCellCharter$inputs$official_autok_sha256 <-
    cellcharter_file_sha256(official_autok_zip)
}
labels_updated <- !"reference_cluster_cellcharter" %in% colnames(colData(sce))
if (labels_updated) {
  clustered_obs <- .cc_h5ad_dataframe(clustered_h5ad, "obs")
  if (nrow(clustered_obs) != ncol(sce) ||
      !"cluster_cellcharter" %in% colnames(clustered_obs)) {
    stop("Cannot attach the published K=11 sensitivity labels.")
  }
  colData(sce)$reference_cluster_cellcharter <-
    clustered_obs$cluster_cellcharter
}
parameters_updated <- is.null(metadata(sce)$rCellCharter$parameters)
if (parameters_updated) {
  metadata(sce)$rCellCharter$parameters <- requested_parameters
  metadata(sce)$rCellCharter$parameters_provenance <-
    "Backfilled from the first invocation that upgraded this legacy checkpoint."
}
if (reference_updated || autok_provenance_updated || labels_updated || parameters_updated) {
  .atomic_save_rds(sce, sce_path)
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
if (!any(requested_stages %in% c("graph", "aggregation", "autok", "downstream", "report"))) {
  message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
  quit(save = "no", status = 0L)
}

sample_ids <- levels(droplevels(colData(sce)$sample))
rows_by_sample <- stats::setNames(
  lapply(sample_ids, function(id) which(colData(sce)$sample == id)), sample_ids
)
rows_path <- file.path(result_dir, "graph", "rows_by_sample.rds")

graph_dir <- file.path(result_dir, "graph")
untrimmed_dir <- file.path(graph_dir, "untrimmed")
trimmed_dir <- file.path(graph_dir, "trimmed")
graph_manifest_path <- file.path(graph_dir, "graph_manifest.rds")

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
    if (file.exists(path)) {
      graph <- readRDS(path)
    } else {
      # Qhull construction is strictly sample-local, matching Squidpy's
      # scipy.spatial.Delaunay algorithm without cross-section edges.
      graph <- delaunay_graph(coordinates)
      .atomic_save_rds(graph, path)
    }
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
    if (!file.exists(global_path)) {
      trimmed <- remove_long_links(
        graphs[[sample_id]], coordinates_by_sample[[sample_id]],
        percentile = 99, cutoff = cutoff
      )
      .atomic_save_rds(trimmed, global_path)
    }
    r_global_paths[[sample_id]] <- global_path
    published_path <- .sample_path(published_trimmed_dir, sample_id)
    if (!file.exists(published_path)) {
      # The released H5AD empirically used a separate 99th-percentile cutoff
      # within each library. Reconstruct it in R and verify it edge-for-edge
      # against the published adjacency below.
      published_graph <- remove_long_links(
        graphs[[sample_id]], coordinates_by_sample[[sample_id]], percentile = 99
      )
      .atomic_save_rds(published_graph, published_path)
    } else {
      published_graph <- readRDS(published_path)
    }
    r_published_paths[[sample_id]] <- published_path
    official_path <- .sample_path(official_dir, sample_id)
    if (!file.exists(official_path)) {
      official <- .read_h5ad_csr_block(
        clustered_h5ad, "obsp/spatial_connectivities", rows_by_sample[[sample_id]]
      )
      .atomic_save_rds(official, official_path)
    } else {
      official <- readRDS(official_path)
    }
    official_paths[[sample_id]] <- official_path
    r_global <- readRDS(global_path)
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
    r_published_paths = r_published_paths,
    r_pinned_global_paths = r_global_paths,
    official_paths = official_paths,
    trimmed_paths = selected_paths
  )
  .atomic_save_rds(rows_by_sample, rows_path)
  .atomic_save_rds(manifest, graph_manifest_path)
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
if (!file.exists(rows_path)) .atomic_save_rds(rows_by_sample, rows_path)
if (!any(requested_stages %in% c("aggregation", "autok", "downstream", "report"))) {
  message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
  quit(save = "no", status = 0L)
}

aggregation_dir <- file.path(
  result_dir,
  if (identical(graph_source, "r_pinned_global")) {
    "aggregation_pinned_global"
  } else {
    "aggregation"
  }
)
aggregation_manifest_path <- file.path(aggregation_dir, "aggregation_manifest.rds")
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
    file.exists(existing_aggregation_manifest$published_graph_parity_checkpoint)
}
if ("aggregation" %in% requested_stages && !aggregation_checkpoint_complete) {
  started <- Sys.time()
  dir.create(aggregation_dir, recursive = TRUE, showWarnings = FALSE)
  sample_files <- stats::setNames(character(length(sample_ids)), sample_ids)
  for (sample_id in sample_ids) {
    rows <- rows_by_sample[[sample_id]]
    path <- .sample_path(aggregation_dir, sample_id, ".h5")
    sample_files[[sample_id]] <- path
    complete_path <- paste0(path, ".complete.rds")
    valid_checkpoint <- file.exists(path) && file.exists(complete_path)
    if (valid_checkpoint) {
      marker <- readRDS(complete_path)
      existing <- try(HDF5Array::HDF5Array(path, "X_cellcharter"), silent = TRUE)
      valid_checkpoint <- !inherits(existing, "try-error") &&
        identical(as.integer(dim(existing)), c(length(rows), 40L)) &&
        identical(marker$sha256, cellcharter_file_sha256(path))
    }
    if (!valid_checkpoint) {
      x_sample <- reducedDim(sce, "X_trVAE")[rows, , drop = FALSE]
      temporary <- tempfile(
        paste0(".", basename(path), "-"), dirname(path), fileext = ".h5"
      )
      cellcharter_aggregate_hdf5(
        x_sample, stats::setNames(list(graphs[[sample_id]]), sample_id),
        stats::setNames(list(seq_along(rows)), sample_id), temporary,
        name = "X_cellcharter", n_layers = 3L
      )
      produced <- HDF5Array::HDF5Array(temporary, "X_cellcharter")
      if (!identical(as.integer(dim(produced)), c(length(rows), 40L))) {
        stop("Incomplete aggregation output for sample ", sample_id, ".")
      }
      if (!file.rename(temporary, path)) {
        stop("Could not atomically publish aggregation checkpoint: ", path)
      }
      .atomic_save_rds(list(
        dimensions = dim(produced), sha256 = cellcharter_file_sha256(path),
        finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
      ), complete_path)
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
      published_aggregate <- cellcharter_aggregate(
        reducedDim(sce, "X_trVAE")[sample_rows, , drop = FALSE],
        published_graph, n_layers = 3L
      )
      published_differences[[i]] <-
        as.matrix(published_aggregate[selected_local_rows, , drop = FALSE]) -
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
  data_id <- digest::digest(list(
    sample_sha256 = sample_sha256,
    graph_manifest_sha256 = cellcharter_file_sha256(graph_manifest_path),
    dataset = "X_cellcharter", n_layers = 3L,
    analysis_graph_source = graph_manifest$analysis_graph_source
  ), algo = "sha256")
  manifest <- list(
    sample_files = sample_files, dataset = "X_cellcharter", n_layers = 3L,
    dimensions = dim(x_cellcharter), sce_checkpoint = aggregated_sce_path,
    parity_checkpoint = file.path(aggregation_dir, "python_r_parity.rds"),
    published_graph_parity_checkpoint = published_graph_parity_path,
    sample_sha256 = sample_sha256, data_id = data_id,
    analysis_graph_source = graph_manifest$analysis_graph_source
  )
  .atomic_save_rds(manifest, aggregation_manifest_path)
  .record_stage("aggregation", started, list(
    checkpoint = aggregation_manifest_path, parity = parity
  ))
  aggregation_checkpoint_complete <- TRUE
}
if (!aggregation_checkpoint_complete) {
  stop("Aggregation checkpoint bundle is incomplete; run with CODEX_STAGES=aggregation (or the default workflow) to regenerate missing artifacts.")
}
if (!file.exists(aggregation_manifest_path)) stop("Aggregation checkpoint is absent: ", aggregation_manifest_path)
aggregation_manifest <- readRDS(aggregation_manifest_path)
if (!identical(
      aggregation_manifest$analysis_graph_source,
      graph_manifest$analysis_graph_source
    )) {
  exact_published_reconstruction <-
    identical(graph_manifest$analysis_graph_source,
              "R-native Qhull with per-sample published-artifact cutoffs") &&
    all(graph_manifest$graph_parity$published_r_only_edges == 0L) &&
    all(graph_manifest$graph_parity$published_official_only_edges == 0L) &&
    identical(
      aggregation_manifest$analysis_graph_source,
      "official published H5AD adjacency (fixed tutorial-parity input)"
    )
  if (!exact_published_reconstruction) {
    stop("Aggregation checkpoint was produced from a different graph source; rerun it in a new CODEX_RESULTS_DIR.")
  }
  # The R reconstruction is edge-for-edge identical to the graph that produced
  # these HDF5 aggregates, so the expensive aggregate and its AutoK cache ID
  # remain valid. Record the audited provenance migration without changing ID.
  aggregation_manifest$analysis_graph_source <-
    graph_manifest$analysis_graph_source
  aggregation_manifest$data_id_migration <- paste(
    "Retained legacy data_id after edge-exact migration from published H5AD",
    "adjacency to its R-native per-sample Qhull reconstruction."
  )
  .atomic_save_rds(aggregation_manifest, aggregation_manifest_path)
}
representation_sce <- sce
sce <- readRDS(aggregation_manifest$sce_checkpoint)
if (!identical(colnames(sce), colnames(representation_sce))) {
  stop("Aggregated and representation SCE checkpoints have different cell order.")
}
if ("reference_cluster_cellcharter" %in% colnames(colData(representation_sce))) {
  colData(sce)$reference_cluster_cellcharter <-
    colData(representation_sce)$reference_cluster_cellcharter
}
metadata(sce)$rCellCharter <- metadata(representation_sce)$rCellCharter
x_cellcharter <- reducedDim(sce, "X_cellcharter")
if (!any(requested_stages %in% c("autok", "downstream", "report"))) {
  message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
  quit(save = "no", status = 0L)
}

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
k11_labels_path <- file.path(selected_dir, "tutorial_sensitivity_k11_labels.rds")

if ("autok" %in% requested_stages) {
  started <- Sys.time()
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("R `torch` cannot be loaded. Repair libtorch/processx before the AutoK stage; graph and aggregation checkpoints are safe.")
  }
  runtime <- torch_gmm_runtime_info("auto", batch_size, seed, num_threads)
  fit_data_id <- digest::digest(list(
    aggregation_data_id = aggregation_manifest$data_id,
    implementation = "rCellCharter torch full-covariance EM v1",
    batch_size = batch_size, max_iter = gmm_max_iter,
    tolerance = gmm_tolerance,
    covariance_regularization = covariance_regularization,
    device = runtime$device, num_threads = num_threads
  ), algo = "sha256")
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
  .atomic_save_rds(best_artifact$model, selected_model_path)
  .atomic_save_rds(list(
    k = autok$best_k, labels = best_artifact$labels, model_path = best_row$path[[1L]],
    selected_model_path = selected_model_path,
    role = "primary R-selected result"
  ), labels_path)
  k11_row <- autok$best_index[autok$best_index$k == 11L, , drop = FALSE]
  if (nrow(k11_row) == 1L) {
    k11_artifact <- readRDS(k11_row$path[[1L]])
    .atomic_save_rds(list(
      k = 11L, labels = k11_artifact$labels, model_path = k11_row$path[[1L]],
      role = "tutorial-matched sensitivity only"
    ), k11_labels_path)
  }
  .record_stage("autok", started, list(
    checkpoint = autok_path, selected_k = autok$best_k, runtime = runtime,
    candidate_bounds = autok_bounds, fitted_grid = autok$fit_k
  ))
}

downstream_dir <- file.path(
  result_dir,
  if (identical(graph_source, "r_pinned_global")) "downstream" else "downstream_published"
)
downstream_manifest_path <- file.path(downstream_dir, "downstream_manifest.rds")
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
}
if ("downstream" %in% requested_stages && !downstream_checkpoint_complete) {
  if (!file.exists(labels_path)) stop("Selected labels are absent; complete AutoK first.")
  started <- Sys.time()
  dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- readRDS(labels_path)
  labels <- factor(selected$labels, levels = seq.int(0L, selected$k - 1L))
  colData(sce)$r_domain <- labels

  k11_agreement_path <- file.path(
    downstream_dir, "tutorial_sensitivity_k11_reference_agreement.rds"
  )
  if (file.exists(k11_labels_path) &&
      "reference_cluster_cellcharter" %in% colnames(colData(sce))) {
    k11_selected <- readRDS(k11_labels_path)
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
  .record_stage("downstream", started, list(
    checkpoint = downstream_manifest_path, selected_k = selected$k,
    graph_bytes_measured = graph_bytes, n_perms = n_perms
  ))
  downstream_checkpoint_complete <- TRUE
}

if ("report" %in% requested_stages) {
  if (!downstream_checkpoint_complete && !file.exists(downstream_manifest_path)) {
    stop("Downstream checkpoints are absent; run with CODEX_STAGES=downstream,report.")
  }
  if (!downstream_checkpoint_complete) {
    stop("Downstream manifest is incomplete; run with CODEX_STAGES=downstream,report to regenerate missing artifacts.")
  }
  started <- Sys.time()
  output_dir <- file.path(result_dir, "report")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render(
    "reports/codex_mouse_spleen_report.Rmd",
    output_file = "codex_mouse_spleen_report.html", output_dir = output_dir,
    params = list(result_dir = result_dir), envir = new.env(parent = globalenv()),
    quiet = FALSE
  )
  .record_stage("report", started, list(
    html = file.path(output_dir, "codex_mouse_spleen_report.html")
  ))
}

message("Requested CODEX stages completed: ", paste(requested_stages, collapse = ", "))
