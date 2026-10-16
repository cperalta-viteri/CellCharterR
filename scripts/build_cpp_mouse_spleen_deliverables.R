#!/usr/bin/env Rscript

# Build lightweight report/package artifacts from the completed compiled-C++
# mouse-spleen AutoK run. This script never refits a GMM or reruns AutoK.
suppressPackageStartupMessages(library(Matrix))
source("R/enrichment.R")

cpp_dir <- "results/CellCharterR_mouse_spleen_cpp_v1"
workflow_dir <- "results/CellCharterR_mouse_spleen_v2"
output_dir <- file.path(cpp_dir, "report_artifacts")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cpp_dir <- normalizePath(cpp_dir, winslash = "/", mustWork = TRUE)
workflow_dir <- normalizePath(workflow_dir, winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

required <- c(
  autok = file.path(cpp_dir, "autok_result.rds"),
  provenance = file.path(cpp_dir, "provenance.rds"),
  primary_labels = file.path(cpp_dir, "selected_labels.rds"),
  primary_model = file.path(cpp_dir, "selected_model.rds"),
  k11_labels = file.path(cpp_dir, "tutorial_k11_labels.rds"),
  k11_model = file.path(cpp_dir, "tutorial_k11_model.rds"),
  report_input = file.path(workflow_dir, "downstream", "report_input.rds"),
  graph_manifest = file.path(workflow_dir, "graph", "graph_manifest.rds")
)
if (any(!file.exists(required))) {
  stop("Missing required artifact(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

atomic_save_rds <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), dirname(path), ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (!file.rename(temporary, path)) stop("Could not publish `", path, "`.")
  invisible(path)
}

autok <- readRDS(required[["autok"]])
provenance <- readRDS(required[["provenance"]])
primary <- readRDS(required[["primary_labels"]])
k11 <- readRDS(required[["k11_labels"]])
report_input <- readRDS(required[["report_input"]])
graph_manifest <- readRDS(required[["graph_manifest"]])

n_spots <- nrow(report_input)
cell_ids <- as.character(report_input$cell_id)
cell_order_sha256 <- digest::digest(cell_ids, algo = "sha256")
stopifnot(
  identical(cell_order_sha256, graph_manifest$cell_order_sha256),
  identical(autok$backend, "cpp"),
  isTRUE(autok$converged),
  identical(autok$best_k, 10L),
  identical(provenance$selected_k, 10L),
  identical(primary$k, 10L),
  identical(primary$seed, 0L),
  identical(k11$k, 11L),
  identical(k11$seed, 12345L),
  length(primary$labels) == n_spots,
  length(k11$labels) == n_spots
)

# Reassemble the validated sample-local graph without touching assays or the
# 40-column aggregated representation.
sample_order <- unique(as.character(report_input$sample))
graph_paths <- graph_manifest$r_pinned_global_paths[sample_order]
graphs <- lapply(sample_order, function(sample_id) {
  rows <- which(as.character(report_input$sample) == sample_id)
  graph <- readRDS(graph_paths[[sample_id]])
  stopifnot(
    nrow(graph) == length(rows),
    identical(rownames(graph), cell_ids[rows]),
    identical(colnames(graph), cell_ids[rows])
  )
  graph
})
combined_graph <- do.call(Matrix::bdiag, graphs)
condition <- factor(as.character(report_input$condition), levels = c("BALBc", "MRL"))
library_id <- factor(as.character(report_input$sample), levels = sample_order)
cell_type <- factor(
  as.character(report_input$cell_type),
  levels = levels(report_input$cell_type)
)

build_role <- function(name, label_object, model_path) {
  k <- as.integer(label_object$k)
  labels <- factor(label_object$labels, levels = seq.int(0L, k - 1L))
  role_dir <- file.path(output_dir, name)
  dir.create(role_dir, recursive = TRUE, showWarnings = FALSE)

  enrichment <- .cc_enrichment(cell_type, labels, log2_fold_change = TRUE)
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
  differential <- .cc_diff_nhood_enrichment(
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

  display <- as.character(report_input$sample) %in% c("BALBc-1", "MRL-4")
  map_input <- report_input[display, c("cell_id", "sample", "condition", "cell_type", "x", "y")]
  map_input$domain <- labels[display]
  attr(map_input, "artifact_metadata") <- list(
    backend = provenance$backend,
    k = k,
    seed = label_object$seed,
    cell_order_sha256 = cell_order_sha256,
    graph_source = graph_manifest$analysis_graph_source
  )

  paths <- c(
    labels = file.path(role_dir, "labels.rds"),
    map_input = file.path(role_dir, "map_input.rds"),
    enrichment = file.path(role_dir, "full_cell_type_enrichment.rds"),
    nhood = file.path(role_dir, "condition_neighborhood_enrichment.rds"),
    differential_inference = file.path(role_dir, "differential_inference.rds"),
    differential_mrl_balbc = file.path(role_dir, "MRL_vs_BALBc_display.rds"),
    differential_balbc_mrl = file.path(role_dir, "BALBc_vs_MRL_display.rds")
  )
  label_artifact <- list(
    format_version = 1L,
    labels = as.integer(label_object$labels),
    k = k,
    cell_order_sha256 = cell_order_sha256,
    backend = provenance$backend,
    seed = as.integer(label_object$seed),
    graph_source = graph_manifest$analysis_graph_source,
    reference_version = provenance$reference_version,
    reference_commit = provenance$reference_commit
  )
  atomic_save_rds(label_artifact, paths[["labels"]], compress = "xz")
  atomic_save_rds(map_input, paths[["map_input"]])
  atomic_save_rds(enrichment, paths[["enrichment"]])
  atomic_save_rds(nhood, paths[["nhood"]])
  atomic_save_rds(differential, paths[["differential_inference"]])
  atomic_save_rds(
    .cc_diff_nhood_display(differential, c("MRL", "BALBc")),
    paths[["differential_mrl_balbc"]]
  )
  atomic_save_rds(
    .cc_diff_nhood_display(differential, c("BALBc", "MRL")),
    paths[["differential_balbc_mrl"]]
  )
  list(paths = paths, labels = label_artifact, model_path = model_path)
}

primary_role <- build_role("primary_AutoK_K10", primary, required[["primary_model"]])
k11_role <- build_role("tutorial_standalone_K11", k11, required[["k11_model"]])

# Render a stable diagnostic directly from the saved AutoK state.
diagnostic_path <- file.path(output_dir, "autok_stability_k2_k20.png")
grDevices::png(diagnostic_path, width = 1800, height = 1000, res = 150)
graphics::par(mar = c(5, 5, 3, 1))
k_values <- as.integer(names(autok$stability_mean))
graphics::plot(
  k_values, autok$stability_mean, type = "b", pch = 19, lwd = 2,
  xlab = "Candidate K", ylab = "Mean neighboring-K Fowlkes-Mallows stability",
  main = "Compiled C++ AutoK stability: full mouse-spleen representation",
  xaxt = "n", ylim = range(c(autok$stability, autok$stability_mean), na.rm = TRUE)
)
graphics::axis(1, at = k_values)
for (i in seq_along(k_values)) {
  values <- autok$stability[i, ]
  graphics::points(rep(k_values[[i]], sum(is.finite(values))), values[is.finite(values)],
                   pch = 16, col = grDevices::adjustcolor("grey30", 0.45), cex = 0.75)
}
graphics::abline(v = 10, col = "#0F766E", lwd = 2, lty = 2)
graphics::abline(v = 11, col = "#8B0000", lwd = 2, lty = 3)
graphics::legend(
  "bottomleft", c("Mean stability", "AutoK primary K=10", "Sensitivity K=11"),
  col = c("black", "#0F766E", "#8B0000"), lty = c(1, 2, 3),
  pch = c(19, NA, NA), bty = "n"
)
grDevices::dev.off()

role_registry <- list(
  Primary = list(
    role_id = "r_cpp_autok_primary", role = "Compiled-C++ AutoK-selected primary",
    k = 10L, seed = 0L, model_path = required[["primary_model"]],
    cell_order_sha256 = cell_order_sha256
  ),
  TutorialK11 = list(
    role_id = "r_cpp_tutorial_k11_sensitivity",
    role = "Standalone compiled-C++ tutorial-matched sensitivity fit",
    k = 11L, seed = 12345L, model_path = required[["k11_model"]],
    cell_order_sha256 = cell_order_sha256
  )
)

# Preserve an honest before/after timing comparison for the same first five
# AutoK rounds (21 K fits per round) on the same 707,474-by-40 representation.
# The legacy runner did not record a single wall-clock duration, so its value is
# a conservative lower bound spanning the first and last completed checkpoint.
legacy_fit_dir <- file.path(workflow_dir, "autok_pinned_global", "fits")
legacy_five_round_paths <- unlist(lapply(seq_len(5L), function(run) {
  file.path(
    legacy_fit_dir,
    sprintf("k_%03d", seq_len(21L)),
    sprintf("run_%03d.rds", run)
  )
}), use.names = FALSE)
if (any(!file.exists(legacy_five_round_paths))) {
  stop("The preserved five-round legacy benchmark checkpoints are incomplete.")
}
legacy_mtime <- file.info(legacy_five_round_paths)$mtime
legacy_elapsed_lower_bound <- as.numeric(
  difftime(max(legacy_mtime), min(legacy_mtime), units = "secs")
)
legacy_example <- readRDS(legacy_five_round_paths[[1L]])
benchmark <- list(
  spots = n_spots,
  features = 40L,
  fit_k = seq_len(21L),
  rounds = 5L,
  fits = 105L,
  legacy = list(
    backend = "R torch CUDA (sequential K fits)",
    device = legacy_example$model$runtime$gpu_name,
    batch_size = legacy_example$model$runtime$batch_size,
    elapsed_lower_bound_seconds = legacy_elapsed_lower_bound,
    timing_source = "span of completed checkpoint modification times"
  ),
  cpp = list(
    backend = "compiled C++ CPU (concurrent K fits)",
    workers = provenance$workers,
    threads_per_fit = provenance$threads_per_fit,
    elapsed_seconds = provenance$elapsed_seconds,
    timing_source = "instrumented end-to-end runner"
  ),
  speedup_lower_bound = legacy_elapsed_lower_bound / provenance$elapsed_seconds,
  limitation = paste(
    "This is an observed workflow comparison, not a controlled microbenchmark.",
    "The legacy duration omits time before its first completed fit;",
    "the backends also use different initializers and numerical engines."
  )
)
manifest <- list(
  tool = "CellCharterR",
  status = "Compiled-C++ full-data AutoK and downstream report inputs",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source_git_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  reference_version = provenance$reference_version,
  reference_commit = provenance$reference_commit,
  backend = provenance$backend,
  backend_algorithm = provenance$backend_algorithm,
  graph_source = graph_manifest$analysis_graph_source,
  spots = n_spots,
  samples = sample_order,
  displayed_samples = c("BALBc-1", "MRL-4"),
  primary_k = 10L,
  sensitivity_k = 11L,
  autok = list(
    requested_k = autok$requested_k,
    max_runs = autok$max_runs,
    runs_completed = autok$runs_completed,
    convergence_tol = autok$convergence_tol,
    final_stability_change = tail(autok$convergence_history$stability_change, 1L),
    elapsed_seconds = provenance$elapsed_seconds,
    workers = provenance$workers,
    threads_per_fit = provenance$threads_per_fit
  ),
  autok_stability_plot = diagnostic_path,
  benchmark = benchmark,
  roles = role_registry,
  enrichment_plot_parameters = list(
    group_key = "cell_type", size_threshold = 0.58,
    color_threshold = 0.58, cluster_labels = TRUE
  ),
  differential_parameters = list(
    significance = 0.05, library_unit = "sample", n_perms = 100L,
    seed = 0L, permutation_mode = "cellcharter_v0.2.0"
  ),
  artifacts = list(Primary = primary_role$paths, TutorialK11 = k11_role$paths),
  session_info = utils::sessionInfo()
)
atomic_save_rds(manifest, file.path(output_dir, "CellCharterR_cpp_manifest.rds"))

# Package only compact, independently licensed handoff/statistical artifacts.
file.copy(
  primary_role$paths[c("labels", "enrichment", "differential_inference")],
  c(
    "inst/extdata/mouse_spleen_cpp_k10_labels.rds",
    "inst/extdata/mouse_spleen_cpp_k10_enrichment.rds",
    "inst/extdata/mouse_spleen_cpp_k10_diff_nhood.rds"
  ),
  overwrite = TRUE
)

message("Compiled-C++ report artifacts written to: ", output_dir)
