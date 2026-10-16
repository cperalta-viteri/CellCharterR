#!/usr/bin/env Rscript

# Rebuild only the inexpensive 1,000-spot teaching workflow with the supported
# compiled C++ backend. This does not read or recompute the full mouse spleen.
suppressPackageStartupMessages(library(CellCharterR))

spe <- cellcharter_mouse_spleen_tiny(with_results = FALSE)
spe <- cellcharter_graph(
  spe, spatial = "spatial", sample_col = "sample",
  graph_name = "cellcharter_graph", percentile = 99,
  progress = FALSE
)
spe <- cellcharter_aggregate(
  spe, representation = "X_trVAE", graph_name = "cellcharter_graph",
  sample_col = "sample", output = "X_cellcharter_tiny", n_layers = 3L,
  progress = FALSE
)
spe <- cellcharter_cluster(
  spe,
  representation = "X_cellcharter_tiny",
  candidate_k = c(2L, 10L),
  checkpoint_dir = tempfile("tiny-cpp-no-checkpoints-"),
  domain_col = "tiny_domain",
  max_runs = 3L,
  convergence_tol = 0.01,
  seed = 1L,
  backend = "cpp",
  workers = 2L,
  threads_per_fit = 1L,
  checkpoint_mode = "none",
  fit_args = list(max_iter = 20L)
)
spe <- cellcharter_enrichment(
  spe, group_col = "cell_type_reduced", domain_col = "tiny_domain",
  progress = FALSE
)
spe <- cellcharter_nhood_enrichment(
  spe, graph_name = "cellcharter_graph", domain_col = "tiny_domain"
)
spe <- cellcharter_diff_nhood_enrichment(
  spe,
  graph_name = "cellcharter_graph",
  domain_col = "tiny_domain",
  condition_col = "condition",
  sample_col = "sample",
  condition_groups = c("BALBc", "MRL"),
  n_perms = 10L,
  seed = 1L,
  progress = FALSE
)

# Store constructor-neutral components so the packaged object remains robust
# across compatible Bioconductor releases.
col_data <- as.data.frame(SummarizedExperiment::colData(spe))
colnames(col_data)[colnames(col_data) == "i-niche"] <- "i.niche"
components <- list(
  assays = lapply(SummarizedExperiment::assays(spe), as.matrix),
  row_data = as.data.frame(SummarizedExperiment::rowData(spe)),
  col_data = col_data,
  row_names = rownames(spe),
  col_names = colnames(spe),
  spatial = as.matrix(SpatialExperiment::spatialCoords(spe)),
  reduced_dims = lapply(
    SingleCellExperiment::reducedDims(spe, withDimnames = FALSE), as.matrix
  ),
  metadata = S4Vectors::metadata(spe),
  col_pairs = lapply(SingleCellExperiment::colPairNames(spe), function(name) {
    pair <- SingleCellExperiment::colPair(spe, name)
    data.frame(
      from = S4Vectors::from(pair),
      to = S4Vectors::to(pair),
      as.data.frame(S4Vectors::mcols(pair)),
      check.names = FALSE
    )
  })
)
names(components$col_pairs) <- SingleCellExperiment::colPairNames(spe)

path <- "inst/extdata/codex_mouse_spleen_tiny_workflow_components.rds"
temporary <- paste0(path, ".partial")
saveRDS(components, temporary, compress = "xz", version = 3L)
if (!file.rename(temporary, path)) stop("Could not publish `", path, "`.")
message("Rebuilt compiled-C++ teaching artifact: ", path)
