#!/usr/bin/env Rscript

# Reconstruct the two-file CODEX import layout from the already-packaged tiny
# SCE. This script never opens or depends on either full H5AD.

input <- "data/codex_mouse_spleen_tiny.rda"
raw_path <- "inst/extdata/codex_mouse_spleen_tiny_raw.h5ad"
clustered_path <- "inst/extdata/codex_mouse_spleen_tiny_clustered.h5ad"
if (!file.exists(input)) stop("Packaged tiny SCE is absent: ", input)
if (file.exists(raw_path) || file.exists(clustered_path)) {
  stop("Tiny H5AD output already exists; remove it explicitly before regeneration.")
}

environment <- new.env(parent = emptyenv())
load(input, envir = environment)
sce <- environment$codex_mouse_spleen_tiny
ids <- colnames(sce)
markers <- rownames(sce)
col_data <- as.data.frame(SummarizedExperiment::colData(sce))

write_obs <- function(path, data, ids) {
  rhdf5::h5createGroup(path, "obs")
  rhdf5::h5write(ids, path, "obs/_index")
  for (name in colnames(data)) {
    value <- data[[name]]
    if (is.factor(value)) value <- as.character(value)
    if (is.logical(value)) value <- as.integer(value)
    rhdf5::h5write(value, path, paste0("obs/", name))
  }
}

rhdf5::h5createFile(raw_path)
write_obs(raw_path, col_data, ids)
rhdf5::h5createGroup(raw_path, "var")
rhdf5::h5write(markers, raw_path, "var/_index")
rhdf5::h5write(
  as.matrix(SummarizedExperiment::assay(sce, "protein")), raw_path, "X"
)
rhdf5::h5createGroup(raw_path, "obsm")
rhdf5::h5write(
  t(as.matrix(SingleCellExperiment::reducedDim(sce, "spatial"))),
  raw_path, "obsm/spatial"
)

rhdf5::h5createFile(clustered_path)
clustered_obs <- data.frame(
  sample = as.character(col_data$sample),
  cluster_cellcharter = as.character(col_data$reference_cluster_cellcharter),
  stringsAsFactors = FALSE
)
write_obs(clustered_path, clustered_obs, ids)
rhdf5::h5createGroup(clustered_path, "obsm")
rhdf5::h5write(
  t(as.matrix(SingleCellExperiment::reducedDim(sce, "X_trVAE"))),
  clustered_path, "obsm/X_trVAE"
)

manifest <- data.frame(
  artifact = basename(c(raw_path, clustered_path)),
  size_bytes = file.info(c(raw_path, clustered_path))$size,
  sha256 = vapply(
    c(raw_path, clustered_path), digest::digest, character(1),
    algo = "sha256", serialize = FALSE, USE.NAMES = FALSE
  ),
  source_tiny_sha256 = digest::digest(
    file = input, algo = "sha256", serialize = FALSE
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest, "inst/extdata/codex_mouse_spleen_tiny_h5ad_manifest.csv",
  row.names = FALSE
)
print(manifest)
