#!/usr/bin/env Rscript

# Create the two data artifacts without recomputing any scientific analysis:
# (1) a movable, self-contained full SpatialExperiment and (2) its
# deterministic tiny subset. This script never reruns the scientific analysis.

args <- commandArgs(trailingOnly = TRUE)
source_rds <- if (length(args) >= 1L) args[[1L]] else
  "results/CellCharterR_mouse_spleen_v2/downstream/codex_sce_with_domains.rds"
k11_rds <- if (length(args) >= 2L) args[[2L]] else
  "results/CellCharterR_mouse_spleen_v2/tutorial_k11_standalone/labels.rds"
output_dir <- if (length(args) >= 3L) args[[3L]] else
  "results/CellCharterR_mouse_spleen_v2/portable"
tiny_n <- if (length(args) >= 4L) as.integer(args[[4L]]) else 1000L
tiny_seed <- if (length(args) >= 5L) as.integer(args[[5L]]) else 20261012L

for (path in c(source_rds, k11_rds)) {
  if (!file.exists(path)) stop("Required existing checkpoint is absent: ", path)
}
if (is.na(tiny_n) || tiny_n < 1L || is.na(tiny_seed)) {
  stop("Tiny size and seed must be valid whole numbers.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

source_sha256 <- digest::digest(file = source_rds, algo = "sha256", serialize = FALSE)
k11_sha256 <- digest::digest(file = k11_rds, algo = "sha256", serialize = FALSE)
source <- readRDS(source_rds)
if (!inherits(source, "SingleCellExperiment")) {
  stop("The source checkpoint is not a SingleCellExperiment.")
}
ids <- colnames(source)
if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) {
  stop("The source SCE must have complete, unique cell identifiers.")
}
if (tiny_n > ncol(source)) stop("Tiny size exceeds the number of source spots.")

k11 <- readRDS(k11_rds)
if (!identical(k11$cell_ids, ids) || length(k11$labels) != length(ids) || anyNA(k11$labels)) {
  stop("Standalone K=11 labels do not exactly align with the source SCE.")
}

col_data <- SummarizedExperiment::colData(source)
if (!"r_domain" %in% colnames(col_data) || length(unique(col_data$r_domain)) != 4L) {
  stop("The historical source checkpoint does not contain its validated K=4 labels.")
}
col_data$r_domain_k4 <- col_data$r_domain
col_data$r_domain_k11_sensitivity <- factor(
  k11$labels, levels = sort(unique(k11$labels))
)

# These are the only large realization calls. They copy already-computed arrays
# into ordinary in-object matrices; no graph, clustering, or statistic is rerun.
protein <- as.matrix(SummarizedExperiment::assay(source, "protein"))
source_reduced_names <- setdiff(
  SingleCellExperiment::reducedDimNames(source), "spatial"
)
reduced <- lapply(source_reduced_names, function(name) {
  value <- as.matrix(SingleCellExperiment::reducedDim(source, name))
  rownames(value) <- ids
  value
})
names(reduced) <- source_reduced_names
spatial <- if (inherits(source, "SpatialExperiment")) {
  as.matrix(SpatialExperiment::spatialCoords(source))
} else {
  as.matrix(SingleCellExperiment::reducedDim(source, "spatial"))
}
rownames(spatial) <- ids
colnames(spatial) <- c("x", "y")

metadata <- S4Vectors::metadata(source)
cc <- metadata$CellCharterR
if (!is.null(cc$inputs)) {
  # Keep immutable provenance but remove machine-local absolute source paths.
  if (!is.null(cc$inputs$raw_h5ad)) {
    cc$inputs$raw_filename <- basename(cc$inputs$raw_h5ad)
    cc$inputs$raw_h5ad <- NULL
  }
  if (!is.null(cc$inputs$clustered_h5ad)) {
    cc$inputs$clustered_filename <- basename(cc$inputs$clustered_h5ad)
    cc$inputs$clustered_h5ad <- NULL
  }
  if (!is.null(cc$inputs$official_autok_zip)) {
    cc$inputs$official_autok_filename <- basename(cc$inputs$official_autok_zip)
    cc$inputs$official_autok_zip <- NULL
  }
}
cc$portable_artifact <- list(
  self_contained = TRUE,
  external_array_dependencies = FALSE,
  source_checkpoint_filename = basename(source_rds),
  source_checkpoint_sha256 = source_sha256,
  k11_checkpoint_filename = basename(k11_rds),
  k11_checkpoint_sha256 = k11_sha256,
  k11_role = k11$role,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
metadata$CellCharterR <- cc

portable <- SpatialExperiment::SpatialExperiment(
  assays = list(protein = protein),
  rowData = SummarizedExperiment::rowData(source),
  colData = col_data,
  reducedDims = S4Vectors::SimpleList(reduced),
  spatialCoords = spatial,
  sample_id = as.character(col_data$sample),
  metadata = metadata
)
stopifnot(
  identical(colnames(portable), ids),
  inherits(portable, "SpatialExperiment"),
  is.matrix(SummarizedExperiment::assay(portable, "protein")),
  is.matrix(SpatialExperiment::spatialCoords(portable)),
  !"spatial" %in% SingleCellExperiment::reducedDimNames(portable),
  all(vapply(
    SingleCellExperiment::reducedDimNames(portable),
    function(name) is.matrix(SingleCellExperiment::reducedDim(portable, name)),
    logical(1)
  ))
)

full_path <- file.path(output_dir, "CellCharterR_mouse_spleen_sce.rds")
full_tmp <- paste0(full_path, ".partial")
saveRDS(portable, full_tmp, compress = "gzip", version = 3L)
if (!file.rename(full_tmp, full_path)) stop("Could not publish the full portable SCE.")
full_sha256 <- digest::digest(file = full_path, algo = "sha256", serialize = FALSE)

# Select the teaching spots directly from the portable object, never from the
# H5AD or a separately reconstructed representation.
set.seed(tiny_seed)
tiny_index <- sort(sample.int(ncol(portable), tiny_n, replace = FALSE))
codex_mouse_spleen_tiny <- portable[, tiny_index]
tiny_ids <- colnames(codex_mouse_spleen_tiny)
tiny_metadata <- S4Vectors::metadata(codex_mouse_spleen_tiny)
tiny_metadata$CellCharterR$tiny_subset <- list(
  source_artifact = basename(full_path),
  source_sha256 = full_sha256,
  seed = tiny_seed,
  n = tiny_n,
  source_column_indices = tiny_index,
  cell_ids = tiny_ids,
  cell_order_sha256 = digest::digest(tiny_ids, algo = "sha256")
)
S4Vectors::metadata(codex_mouse_spleen_tiny) <- tiny_metadata

tiny_rds <- file.path(output_dir, "CellCharterR_mouse_spleen_tiny_sce.rds")
saveRDS(codex_mouse_spleen_tiny, tiny_rds, compress = "xz", version = 3L)
save(codex_mouse_spleen_tiny, file = "data/codex_mouse_spleen_tiny.rda",
     compress = "xz", version = 3L)

manifest <- list(
  format_version = 1L,
  full = list(
    path = normalizePath(full_path),
    sha256 = full_sha256,
    size_bytes = file.info(full_path)$size,
    spots = ncol(portable), features = nrow(portable)
  ),
  tiny = list(
    rds_path = normalizePath(tiny_rds),
    rda_path = normalizePath("data/codex_mouse_spleen_tiny.rda"),
    rds_sha256 = digest::digest(file = tiny_rds, algo = "sha256", serialize = FALSE),
    rda_sha256 = digest::digest(file = "data/codex_mouse_spleen_tiny.rda",
                                algo = "sha256", serialize = FALSE),
    spots = tiny_n, seed = tiny_seed, spot_ids = tiny_ids
  )
)
saveRDS(manifest, file.path(output_dir, "portable_data_manifest.rds"))

cat("Full portable SpatialExperiment:", full_path, "\n")
cat("Full SHA-256:", full_sha256, "\n")
cat("Tiny packaged SpatialExperiment: data/codex_mouse_spleen_tiny.rda\n")
cat("Tiny spots:", tiny_n, "seed:", tiny_seed, "\n")
