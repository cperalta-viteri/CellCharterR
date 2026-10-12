#!/usr/bin/env Rscript

# Build report inputs only after the revised workflow has produced three
# explicitly distinct label roles. The frozen Phase-0 tree at
# results/CellCharterR_mouse_spleen is intentionally not read or modified.
suppressPackageStartupMessages({
  library(Matrix)
})
source("R/enrichment.R")

workflow_request <- Sys.getenv(
  "CELLCHARTERR_WORKFLOW_RESULTS",
  "results/CellCharterR_mouse_spleen_v2"
)
if (!grepl("^/", workflow_request)) {
  workflow_request <- file.path(getwd(), workflow_request)
}
workflow_dir <- normalizePath(workflow_request, winslash = "/", mustWork = TRUE)

output_request <- Sys.getenv(
  "CELLCHARTERR_REPORT_RESULTS",
  file.path(workflow_dir, "second_review_deliverables")
)
if (!grepl("^/", output_request)) output_request <- file.path(getwd(), output_request)
output_dir <- normalizePath(output_request, winslash = "/", mustWork = FALSE)
phase0_dir <- normalizePath(
  "results/CellCharterR_mouse_spleen", winslash = "/", mustWork = FALSE
)
if (identical(workflow_dir, phase0_dir) || identical(output_dir, phase0_dir) ||
    startsWith(output_dir, paste0(phase0_dir, "/"))) {
  stop(
    "The Phase-0 result root is frozen and provisional/stale; choose a new ",
    "CELLCHARTERR_REPORT_RESULTS path."
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), dirname(path), ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("Could not checkpoint: ", path)
  invisible(path)
}

required <- c(
  report_input = file.path(workflow_dir, "downstream", "report_input.rds"),
  primary_labels = file.path(workflow_dir, "selected", "selected_labels.rds"),
  tutorial_k11_labels = file.path(workflow_dir, "tutorial_k11_standalone", "labels.rds"),
  published_reference = file.path(
    workflow_dir, "reference", "published_h5ad_cluster_cellcharter_labels.rds"
  ),
  graph_manifest = file.path(workflow_dir, "graph", "graph_manifest.rds"),
  autok_diagnostic = file.path(
    workflow_dir, "autok_pinned_global", "diagnostics", "stability_k1_k20.png"
  )
)
if (any(!file.exists(required))) {
  stop(
    "Required revised-workflow checkpoint(s) are absent:\n  - ",
    paste(required[!file.exists(required)], collapse = "\n  - "),
    "\nRun the missing workflow stages under the new result root first.",
    call. = FALSE
  )
}

report_input <- readRDS(required[["report_input"]])
primary <- readRDS(required[["primary_labels"]])
tutorial_k11 <- readRDS(required[["tutorial_k11_labels"]])
published_reference <- readRDS(required[["published_reference"]])

validate_role <- function(x, expected_role, context, n_cells) {
  required_fields <- c(
    "role_id", "role", "seed", "model_path", "cell_ids",
    "cell_order_sha256", "aggregation_data_id", "k", "labels"
  )
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields)) {
    stop(context, " is missing: ", paste(missing_fields, collapse = ", "))
  }
  if (!identical(x$role_id, expected_role)) {
    stop(context, " has role `", x$role_id, "`, expected `", expected_role, "`.")
  }
  if (length(x$labels) != n_cells || length(x$cell_ids) != n_cells ||
      anyNA(x$labels) || anyNA(x$cell_ids)) {
    stop(context, " does not contain one identified label per report row.")
  }
  invisible(x)
}

n_cells <- nrow(report_input)
validate_role(primary, "r_autok_primary", "Primary label checkpoint", n_cells)
validate_role(
  tutorial_k11, "r_tutorial_k11_sensitivity",
  "Standalone tutorial K=11 label checkpoint", n_cells
)
validate_role(
  published_reference, "published_h5ad_reference",
  "Published H5AD reference checkpoint", n_cells
)
if (!identical(tutorial_k11$k, 11L) || !identical(tutorial_k11$seed, 12345L) ||
    !identical(tutorial_k11$fit_function, "torch_gmm_fit")) {
  stop("Tutorial K=11 must be a standalone torch_gmm_fit with seed 12345.")
}
if (!is.na(published_reference$seed) || !is.na(published_reference$model_path)) {
  stop("Published H5AD labels must not claim an R fit seed or model path.")
}
if (identical(primary$model_path, tutorial_k11$model_path)) {
  stop("Primary AutoK and standalone tutorial K=11 must not share a model path.")
}
if (!file.exists(primary$model_path) || !file.exists(tutorial_k11$model_path)) {
  stop("An R label role refers to an absent model checkpoint.")
}

role_objects <- list(
  Primary = primary,
  TutorialK11 = tutorial_k11,
  PublishedH5ADReference = published_reference
)
expected_ids <- as.character(report_input$cell_id)
if (anyNA(expected_ids) || any(!nzchar(trimws(expected_ids))) ||
    anyDuplicated(expected_ids)) {
  stop("Report-input cell IDs must be non-missing, non-blank, and unique.")
}
computed_order_sha256 <- digest::digest(expected_ids, algo = "sha256")
if (any(vapply(role_objects, function(x) !identical(x$cell_ids, expected_ids), logical(1)))) {
  stop("A role checkpoint cell ID/order differs from the report input.")
}
role_order_sha256 <- vapply(
  role_objects, `[[`, character(1), "cell_order_sha256"
)
if (length(unique(role_order_sha256)) != 1L ||
    any(role_order_sha256 != computed_order_sha256)) {
  stop("Role checkpoints disagree on the immutable cell-order hash.")
}
role_aggregation_ids <- vapply(
  role_objects, `[[`, character(1), "aggregation_data_id"
)
if (anyNA(role_aggregation_ids) || any(!nzchar(role_aggregation_ids)) ||
    length(unique(role_aggregation_ids)) != 1L) {
  stop("Role checkpoints were not produced from one aggregation identity.")
}
input_metadata <- attr(report_input, "artifact_metadata", exact = TRUE)
if (is.null(input_metadata) ||
    !identical(input_metadata$role_id, primary$role_id) ||
    !identical(input_metadata$seed, primary$seed) ||
    !identical(input_metadata$model_path, primary$model_path) ||
    !identical(input_metadata$cell_order_sha256, primary$cell_order_sha256) ||
    !identical(input_metadata$aggregation_data_id, primary$aggregation_data_id)) {
  stop("Primary report input is missing or mismatches its role metadata.")
}

# Reassemble only the sparse graph. Assays and the aggregate representation stay
# disk-backed and are never materialized by report-input construction.
graph_manifest <- readRDS(required[["graph_manifest"]])
graph_paths <- graph_manifest$r_pinned_global_paths
sample_order <- unique(as.character(report_input$sample))
if (!identical(names(graph_paths), sample_order)) graph_paths <- graph_paths[sample_order]
if (anyNA(graph_paths) || any(!file.exists(graph_paths))) {
  stop("The sample-local pinned graph checkpoint bundle is incomplete.")
}
graphs <- lapply(sample_order, function(sample_id) {
  rows <- which(as.character(report_input$sample) == sample_id)
  graph <- readRDS(graph_paths[[sample_id]])
  if (nrow(graph) != length(rows) || ncol(graph) != length(rows)) {
    stop("Graph dimensions do not match sample ", sample_id, ".")
  }
  if (is.null(rownames(graph)) ||
      !identical(rownames(graph), expected_ids[rows]) ||
      !identical(colnames(graph), expected_ids[rows])) {
    stop("Graph cell order does not match report input for sample ", sample_id, ".")
  }
  graph
})
names(graphs) <- sample_order
combined_graph <- do.call(Matrix::bdiag, graphs)
if (nrow(combined_graph) != n_cells) {
  stop("Combined graph and report input differ in cell count.")
}

condition <- factor(as.character(report_input$condition), levels = c("BALBc", "MRL"))
library_id <- factor(as.character(report_input$sample), levels = sample_order)
# The v0.2.0 tutorial uses the complete categorical cell_type annotation, not
# the reduced publication-analysis grouping. Preserve its stored category order
# exactly; alphabetical reconstruction changes the displayed/reference matrix.
if (!is.factor(report_input$cell_type) || is.null(levels(report_input$cell_type))) {
  stop("`report_input$cell_type` must preserve the official categorical levels.")
}
cell_type <- factor(
  as.character(report_input$cell_type),
  levels = levels(report_input$cell_type)
)
analysis_paths <- list()

analysis_roles <- list(Primary = primary, TutorialK11 = tutorial_k11)
for (analysis_name in names(analysis_roles)) {
  selected <- analysis_roles[[analysis_name]]
  k <- as.integer(selected$k)
  labels <- factor(selected$labels, levels = seq.int(0L, k - 1L))
  analysis_slug <- if (identical(selected$role_id, "r_autok_primary")) {
    paste0("primary_AutoK_K", k)
  } else {
    "tutorial_standalone_K11"
  }
  analysis_dir <- file.path(output_dir, analysis_slug)
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

  enrichment <- .cc_enrichment(
    group = cell_type, label = labels, log2_fold_change = TRUE
  )
  nhood <- lapply(c("BALBc", "MRL"), function(group) {
    keep <- which(condition == group)
    .cc_nhood_enrichment(
      combined_graph[keep, keep, drop = FALSE],
      droplevels(labels[keep]),
      log2_fold_change = FALSE,
      only_inter = TRUE
    )
  })
  names(nhood) <- c("BALBc", "MRL")

  # Infer once in categorical order. The reverse panel negates only the effect
  # and reuses the exact same permutation counts and p-values.
  differential_inference <- .cc_diff_nhood_enrichment(
    adjacency = combined_graph,
    cluster = labels,
    condition = condition,
    library = library_id,
    condition_groups = c("BALBc", "MRL"),
    n_perms = 100L,
    seed = 0L,
    log2_fold_change = FALSE,
    only_inter = TRUE,
    permutation_mode = "cellcharter_v0.2.0"
  )
  differential_balbc_mrl <- .cc_diff_nhood_display(
    differential_inference, c("BALBc", "MRL")
  )
  differential_mrl_balbc <- .cc_diff_nhood_display(
    differential_inference, c("MRL", "BALBc")
  )

  map_input <- report_input[, c(
    "cell_id", "sample", "condition", "cell_type", "x", "y"
  )]
  map_input$domain <- labels
  attr(map_input, "artifact_metadata") <- list(
    role_id = selected$role_id,
    role = selected$role,
    seed = selected$seed,
    model_path = selected$model_path,
    cell_id_column = "cell_id",
    cell_order_sha256 = selected$cell_order_sha256,
    aggregation_data_id = selected$aggregation_data_id
  )
  paths <- c(
    labels = file.path(analysis_dir, "labels_with_role.rds"),
    map_input = file.path(analysis_dir, "map_input_with_role.rds"),
    enrichment = file.path(analysis_dir, "full_cell_type_enrichment.rds"),
    nhood = file.path(analysis_dir, "condition_neighborhood_enrichment.rds"),
    differential_inference = file.path(analysis_dir, "differential_inference.rds"),
    differential_mrl_balbc = file.path(analysis_dir, "MRL_vs_BALBc_display.rds"),
    differential_balbc_mrl = file.path(analysis_dir, "BALBc_vs_MRL_display.rds")
  )
  atomic_save_rds(selected, paths[["labels"]])
  atomic_save_rds(map_input, paths[["map_input"]])
  atomic_save_rds(enrichment, paths[["enrichment"]])
  atomic_save_rds(nhood, paths[["nhood"]])
  atomic_save_rds(differential_inference, paths[["differential_inference"]])
  atomic_save_rds(differential_mrl_balbc, paths[["differential_mrl_balbc"]])
  atomic_save_rds(differential_balbc_mrl, paths[["differential_balbc_mrl"]])
  analysis_paths[[analysis_name]] <- paths
}

role_registry <- lapply(role_objects, function(x) list(
  role_id = x$role_id,
  role = x$role,
  k = x$k,
  seed = x$seed,
  model_path = x$model_path,
  labels_path = switch(
    x$role_id,
    r_autok_primary = required[["primary_labels"]],
    r_tutorial_k11_sensitivity = required[["tutorial_k11_labels"]],
    published_h5ad_reference = required[["published_reference"]]
  ),
  cell_order_sha256 = x$cell_order_sha256,
  aggregation_data_id = x$aggregation_data_id
))

manifest <- list(
  tool = "CellCharterR",
  status = "Phase-3 regenerated report inputs; Phase-0 outputs are provisional/stale",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  workflow_root = workflow_dir,
  output_root = output_dir,
  source_git_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  reference_version = "CellCharter 0.2.0",
  reference_commit = "2a8fdde8213e32e5a00f1f9571bffd3dbb825a79",
  graph_source = graph_manifest$analysis_graph_source,
  graph_size_bytes = as.numeric(object.size(combined_graph)),
  cells = n_cells,
  samples = sample_order,
  displayed_samples = c("BALBc-1", "MRL-4"),
  primary_k = as.integer(primary$k),
  sensitivity_k = 11L,
  autok_stability_plot = required[["autok_diagnostic"]],
  roles = role_registry,
  enrichment_plot_parameters = list(
    group_key = "cell_type",
    size_threshold = 0.58,
    color_threshold = 0.58,
    cluster_labels = TRUE,
    cluster_labels_source = "CellCharter v0.2.0 default",
    statistic = "log2(observed / expected)"
  ),
  proximity_plot_parameters = list(
    groups = c("BALBc", "MRL"), annotate = TRUE,
    neighborhood_limits = c(-1, 1),
    differential_limits = "data range (CellCharter v0.2.0 default)"
  ),
  differential_parameters = list(
    inferential_contrast = c("BALBc", "MRL"),
    categorical_level_order = c("BALBc", "MRL"),
    displays = list(c("MRL", "BALBc"), c("BALBc", "MRL")),
    significance = 0.05,
    library_counts = c(BALBc = 3L, MRL = 6L),
    library_unit = "sample",
    n_perms = 100L,
    seed = 0L,
    permutation_mode = "cellcharter_v0.2.0",
    empirical_pvalue_denominator = 100L,
    reverse_display = "effect-only sign reversal; counts and p-values reused"
  ),
  artifacts = analysis_paths,
  session_info = utils::sessionInfo()
)
manifest_path <- file.path(output_dir, "CellCharterR_primary_K11_manifest.rds")
atomic_save_rds(manifest, manifest_path)

message("Revised primary/K=11 report inputs written to: ", output_dir)
