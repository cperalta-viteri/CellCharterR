#!/usr/bin/env Rscript

# R-native, PCA-based CellCharter reproduction of the CODEX mouse-spleen tutorial.
# The assay stays HDF5-backed; one physical spleen is read at a time.
suppressPackageStartupMessages({ library(anndataR); library(HDF5Array); library(DelayedArray) })
source("R/aggregation.R")
source("R/torch_gmm.R")

input_h5ad <- Sys.getenv("CODEX_H5AD", "../data/codex_mouse_spleen.h5ad")
out_dir <- "results/codex_mouse_spleen"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(12345)

# anndataR reads official metadata and coordinates. Current SCE conversion
# materialises dense X, so HDF5Array reads the H5AD dataset directly instead.
ad <- read_h5ad(input_h5ad, as = "HDF5AnnData")
obs <- ad$obs
coordinates <- ad$obsm$spatial
X <- t(HDF5Array(input_h5ad, "X")) # H5AD X is feature x cell; we need cell x feature.
stopifnot(nrow(X) == nrow(obs), nrow(coordinates) == nrow(obs))

sample_ids <- unique(as.character(obs$sample))
rows_by_sample <- setNames(lapply(sample_ids, function(s) which(obs$sample == s)), sample_ids)

# Learn PCA directions from balanced, sample-wise scaled blocks. PCA replaces
# the tutorial's supplied trVAE embedding as the explicit R-native baseline.
n_pcs <- 10L
train <- lapply(rows_by_sample, function(rows) {
  take <- sample(rows, min(length(rows), 5000L))
  scale(as.matrix(X[take, , drop = FALSE]))
})
pca_fit <- prcomp(do.call(rbind, train), center = FALSE, scale. = FALSE, rank. = n_pcs)
saveRDS(pca_fit, file.path(out_dir, "pca_model.rds"))
rm(train); gc()

# Create on-disk PCA scores. A single HDF5 writer prevents concurrent-write corruption.
pca_file <- file.path(out_dir, "pca_scores.h5")
if (file.exists(pca_file)) file.remove(pca_file)
pca_sink <- HDF5RealizationSink(c(nrow(X), n_pcs), filepath = pca_file, name = "X_pca", type = "double")
write_rows <- function(sink, rows, block, nr, nc) {
  runs <- split(sort(rows), cumsum(c(TRUE, diff(sort(rows)) != 1L)))
  for (run in runs) {
    local <- match(run, rows)
    viewport <- S4Arrays::ArrayViewport(c(nr, nc), IRanges::IRanges(start = c(run[[1]], 1L), width = c(length(run), nc)))
    DelayedArray::write_block(sink, viewport, block[local, , drop = FALSE])
  }
}
for (sample_id in sample_ids) {
  rows <- rows_by_sample[[sample_id]]
  block <- as.matrix(X[rows, , drop = FALSE]) # Only one spleen enters RAM.
  mu <- colMeans(block); sigma <- apply(block, 2, sd); sigma[sigma == 0 | !is.finite(sigma)] <- 1
  projected <- sweep(sweep(block, 2, mu, "-"), 2, sigma, "/") %*% pca_fit$rotation
  write_rows(pca_sink, rows, projected, nrow(X), n_pcs)
  rm(block, projected); gc(FALSE)
}
X_pca <- HDF5Array(pca_file, "X_pca")

# Delaunay graph and top-1% long-edge trimming are strictly sample-local.
graphs <- list()
for (sample_id in sample_ids) {
  rows <- rows_by_sample[[sample_id]]
  graph <- delaunay_graph(coordinates[rows, 1:2, drop = FALSE])
  graphs[[sample_id]] <- remove_long_links(graph, coordinates[rows, 1:2, drop = FALSE], percentile = 99)
}
saveRDS(graphs, file.path(out_dir, "sample_graphs.rds"), compress = FALSE)

# Own PCs plus exclusive 1/2/3-hop neighbour means make a 40-feature niche vector.
niche_file <- file.path(out_dir, "cellcharter_features.h5")
if (file.exists(niche_file)) file.remove(niche_file)
X_niche <- cellcharter_aggregate_hdf5(X_pca, graphs, rows_by_sample, niche_file, name = "X_cellcharter", n_layers = 3L)

# The official tutorial fixes K=11. The Torch EM is full-covariance and streams rows.
gmm <- torch_gmm_fit(X_niche, n_clusters = 11L, batch_size = 8192L, max_iter = 100L,
                     tolerance = 1e-3, covariance_regularization = 1e-5, device = "auto", seed = 12345L,
                     num_threads = parallel::detectCores(logical = FALSE))
saveRDS(gmm, file.path(out_dir, "torch_gmm_k11.rds"))
labels <- torch_gmm_predict(gmm, X_niche, batch_size = 8192L, device = "auto")

# Compact artifacts are enough to render the report without reopening the HDF5 matrix.
result <- data.frame(sample = obs$sample, condition = sub("-.*$", "", as.character(obs$sample)),
                     cell_type = obs$cell_type, i_niche = obs$`i-niche`,
                     x = coordinates[, 1], y = coordinates[, 2], r_domain = factor(labels))
saveRDS(result, file.path(out_dir, "domain_metadata.rds"))
write.csv(as.data.frame(table(result$r_domain, result$cell_type)), file.path(out_dir, "domain_by_cell_type.csv"), row.names = FALSE)
write.csv(as.data.frame(table(result$r_domain, result$sample)), file.path(out_dir, "domain_by_sample.csv"), row.names = FALSE)
message("CODEX workflow completed: ", out_dir)
