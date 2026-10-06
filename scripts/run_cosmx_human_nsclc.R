#!/usr/bin/env Rscript

# Disk-backed R-native reproduction of the post-scVI CellCharter CosMx NSCLC
# workflow.  It never builds a cohort-wide dense matrix: each physical section
# is transformed and graph-aggregated independently before a streamed GMM fit.
suppressPackageStartupMessages({
  library(anndataR); library(HDF5Array); library(DelayedArray)
})
source("R/aggregation.R")
source("R/torch_gmm.R")

input_h5ad <- Sys.getenv("COSMX_H5AD", "../data/cosmx_human_nsclc/cosmx_human_nsclc_clustered.h5ad")
out_dir <- "results/cosmx_human_nsclc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(12345)

# anndataR provides the annotation data frames without treating them as an R
# assay.  HDF5Array accesses numerical H5AD fields lazily.
ad <- read_h5ad(input_h5ad, as = "HDF5AnnData")
obs <- ad$obs
coordinates <- ad$obsm$spatial

# The official processed tutorial file stores the 10-D scVI representation.
# H5AD stores obsm matrices as observation-by-feature, so no transpose is used.
X_scvi <- HDF5Array(input_h5ad, "obsm/X_scVI")
stopifnot(nrow(X_scvi) == nrow(obs), nrow(coordinates) == nrow(obs))

# `sample` is the physical-section identifier.  It is intentionally used for
# graph construction and aggregation; patient IDs must never join sections.
sample_ids <- unique(as.character(obs$sample))
rows_by_sample <- setNames(lapply(sample_ids, function(s) which(obs$sample == s)), sample_ids)
graphs <- list()
for (sample_id in sample_ids) {
  rows <- rows_by_sample[[sample_id]]
  raw_graph <- delaunay_graph(coordinates[rows, 1:2, drop = FALSE])
  graphs[[sample_id]] <- remove_long_links(raw_graph, coordinates[rows, 1:2, drop = FALSE], percentile = 99)
}
saveRDS(graphs, file.path(out_dir, "sample_graphs.rds"), compress = FALSE)

# Own scVI values plus 1/2/3-hop mean neighbourhood features form a 40-D
# CellCharter-like input.  The writer is serial because HDF5 does not permit
# safe concurrent writes to the same file.
niche_file <- file.path(out_dir, "cellcharter_features.h5")
if (file.exists(niche_file)) file.remove(niche_file)
X_niche <- cellcharter_aggregate_hdf5(X_scvi, graphs, rows_by_sample, niche_file,
                                      name = "X_cellcharter", n_layers = 3L)

# K=20 follows the published analysis.  EM uses tensor algebra and CUDA when
# R Torch can see a GPU; on CPU it streams fixed-size disk-to-device batches.
gmm <- torch_gmm_fit(X_niche, n_clusters = 20L, batch_size = 8192L,
                     max_iter = 100L, tolerance = 1e-3,
                     covariance_regularization = 1e-5, device = "auto",
                     seed = 12345L, num_threads = parallel::detectCores(logical = FALSE))
saveRDS(gmm, file.path(out_dir, "torch_gmm_k20.rds"))
labels <- torch_gmm_predict(gmm, X_niche, batch_size = 8192L, device = "auto")

# Keep only a compact metadata result table for report rendering.  The runner
# tests for cluster_20 because older tutorial files may omit official labels.
official <- if ("cluster_20" %in% names(obs)) obs$cluster_20 else NA
result <- data.frame(sample = obs$sample, patient = obs$patient,
                     cell_type = obs$cell_type, official_cluster_20 = official,
                     x = coordinates[, 1], y = coordinates[, 2],
                     r_domain = factor(labels))
saveRDS(result, file.path(out_dir, "domain_metadata.rds"))
write.csv(as.data.frame(table(result$r_domain, result$cell_type)),
          file.path(out_dir, "domain_by_cell_type.csv"), row.names = FALSE)
write.csv(as.data.frame(table(result$r_domain, result$sample)),
          file.path(out_dir, "domain_by_sample.csv"), row.names = FALSE)
message("CosMx workflow completed: ", out_dir)
