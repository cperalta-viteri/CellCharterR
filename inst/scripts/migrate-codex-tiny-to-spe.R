#!/usr/bin/env Rscript

# Rebuild only the portable 1,000-cell teaching object from the packaged H5AD
# pair as a native SpatialExperiment. Existing full-data labels and aggregate
# features are copied from the already validated tiny object; no full-data
# graph, clustering, or downstream statistic is recomputed.

args <- commandArgs(trailingOnly = TRUE)
raw_h5ad <- if (length(args) >= 1L) args[[1L]] else
  "inst/extdata/codex_mouse_spleen_tiny_raw.h5ad"
clustered_h5ad <- if (length(args) >= 2L) args[[2L]] else
  "inst/extdata/codex_mouse_spleen_tiny_clustered.h5ad"
tiny_rds <- if (length(args) >= 3L) args[[3L]] else
  "results/CellCharterR_mouse_spleen_v2/portable/CellCharterR_mouse_spleen_tiny_sce.rds"

required <- c(
  "CellCharterR", "digest", "SingleCellExperiment", "SpatialExperiment",
  "SummarizedExperiment"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "))

legacy_env <- new.env(parent = emptyenv())
utils::data(
  "codex_mouse_spleen_tiny", package = "CellCharterR", envir = legacy_env
)
legacy <- legacy_env$codex_mouse_spleen_tiny
imported <- CellCharterR::read_codex_spe(raw_h5ad, clustered_h5ad)

if (!identical(colnames(imported), colnames(legacy)) ||
    !identical(rownames(imported), rownames(legacy))) {
  stop("Fresh H5AD conversion does not match the validated tiny object identity.")
}

protein <- as.matrix(SummarizedExperiment::assay(imported, "protein"))
trvae <- as.matrix(SingleCellExperiment::reducedDim(imported, "X_trVAE"))
spatial <- as.matrix(SpatialExperiment::spatialCoords(imported))
legacy_protein <- as.matrix(SummarizedExperiment::assay(legacy, "protein"))
legacy_trvae <- as.matrix(SingleCellExperiment::reducedDim(legacy, "X_trVAE"))
legacy_spatial <- if (inherits(legacy, "SpatialExperiment")) {
  as.matrix(SpatialExperiment::spatialCoords(legacy))
} else {
  as.matrix(SingleCellExperiment::reducedDim(legacy, "spatial"))
}
stopifnot(
  identical(protein, legacy_protein),
  identical(trvae, legacy_trvae),
  identical(spatial, legacy_spatial)
)

reduced <- S4Vectors::SimpleList(X_trVAE = trvae)
if ("X_cellcharter" %in% SingleCellExperiment::reducedDimNames(legacy)) {
  reduced$X_cellcharter <- as.matrix(
    SingleCellExperiment::reducedDim(legacy, "X_cellcharter")
  )
}
metadata <- S4Vectors::metadata(legacy)
metadata$CellCharterR$inputs <- S4Vectors::metadata(imported)$CellCharterR$inputs
metadata$CellCharterR$spatialexperiment_migration <- list(
  source = "packaged checksum-pinned tiny H5AD pair",
  coordinate_store = "SpatialExperiment::spatialCoords",
  analytical_recomputation = FALSE,
  migrated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

codex_mouse_spleen_tiny <- SpatialExperiment::SpatialExperiment(
  assays = list(protein = protein),
  rowData = SummarizedExperiment::rowData(imported),
  colData = SummarizedExperiment::colData(legacy),
  reducedDims = reduced,
  spatialCoords = spatial,
  sample_id = as.character(SummarizedExperiment::colData(legacy)$sample),
  metadata = metadata
)
stopifnot(
  inherits(codex_mouse_spleen_tiny, "SpatialExperiment"),
  !"spatial" %in% SingleCellExperiment::reducedDimNames(codex_mouse_spleen_tiny),
  identical(rownames(SpatialExperiment::spatialCoords(codex_mouse_spleen_tiny)),
            colnames(codex_mouse_spleen_tiny))
)

dir.create(dirname(tiny_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(codex_mouse_spleen_tiny, tiny_rds, compress = "xz", version = 3L)
save(codex_mouse_spleen_tiny, file = "data/codex_mouse_spleen_tiny.rda",
     compress = "xz", version = 3L)

cat("Rebuilt native tiny SpatialExperiment:", tiny_rds, "\n")
cat("Packaged data: data/codex_mouse_spleen_tiny.rda\n")
cat("Cells:", ncol(codex_mouse_spleen_tiny), "\n")
cat("RDS SHA-256:", digest::digest(
  file = tiny_rds, algo = "sha256", serialize = FALSE
), "\n")
