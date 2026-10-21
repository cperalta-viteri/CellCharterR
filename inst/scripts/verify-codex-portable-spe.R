#!/usr/bin/env Rscript

# Read-only relocation and integrity verifier for the portable full SPE. This
# script never rebuilds a graph, embedding, clustering, or statistic.

args <- commandArgs(trailingOnly = TRUE)
full_rds <- if (length(args) >= 1L) args[[1L]] else
  "results/CellCharterR_mouse_spleen_v2/portable/CellCharterR_mouse_spleen_sce.rds"
tiny_rds <- if (length(args) >= 2L) args[[2L]] else
  "results/CellCharterR_mouse_spleen_v2/portable/CellCharterR_mouse_spleen_tiny_sce.rds"
output_rds <- if (length(args) >= 3L) args[[3L]] else
  "results/CellCharterR_mouse_spleen_v2/portable/portable_validation_result.rds"

required_packages <- c(
  "digest", "S4Vectors", "SingleCellExperiment", "SpatialExperiment",
  "SummarizedExperiment"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing verifier package(s): ", paste(missing_packages, collapse = ", "))
}
for (path in c(full_rds, tiny_rds)) {
  if (!file.exists(path)) stop("Required portable artifact is absent: ", path)
}

# Copy only the serialized object to a fresh unrelated directory. Reading this
# copy proves that no neighboring H5AD/HDF5 file is required.
relocation_dir <- tempfile("CellCharterR-portability-")
dir.create(relocation_dir)
on.exit(unlink(relocation_dir, recursive = TRUE, force = TRUE), add = TRUE)
relocated_rds <- file.path(relocation_dir, basename(full_rds))
if (!file.copy(full_rds, relocated_rds, copy.mode = FALSE)) {
  stop("Could not create the isolated validation copy.")
}

source_sha256 <- digest::digest(
  file = full_rds, algo = "sha256", serialize = FALSE
)
relocated_sha256 <- digest::digest(
  file = relocated_rds, algo = "sha256", serialize = FALSE
)
if (!identical(source_sha256, relocated_sha256)) {
  stop("The relocated RDS checksum differs from its source.")
}

spe <- readRDS(relocated_rds)
tiny <- readRDS(tiny_rds)
if (!inherits(spe, "SpatialExperiment") ||
    !inherits(tiny, "SpatialExperiment")) {
  stop("Both artifacts must be SpatialExperiment objects.")
}

expected_assays <- "protein"
expected_reduced_dims <- c("X_trVAE", "X_cellcharter")
required_labels <- c(
  "reference_cluster_cellcharter", "r_domain_k4",
  "r_domain_k11_sensitivity"
)
if (!identical(SummarizedExperiment::assayNames(spe), expected_assays)) {
  stop("The full SPE does not contain exactly the required protein assay.")
}
if (!setequal(
  SingleCellExperiment::reducedDimNames(spe), expected_reduced_dims
)) {
  stop("The full SPE lacks a required reduced dimension.")
}
if (!all(required_labels %in% colnames(SummarizedExperiment::colData(spe)))) {
  stop("The full SPE lacks required published, K=4, or K=11 labels.")
}
if (!identical(dim(spe), c(29L, 707474L))) {
  stop("Unexpected full SPE dimensions: ", paste(dim(spe), collapse = " x "))
}

large_fields <- list(
  protein = SummarizedExperiment::assay(spe, "protein"),
  X_trVAE = SingleCellExperiment::reducedDim(spe, "X_trVAE"),
  spatial = SpatialExperiment::spatialCoords(spe),
  X_cellcharter = SingleCellExperiment::reducedDim(spe, "X_cellcharter")
)
if (!all(vapply(large_fields, is.matrix, logical(1)))) {
  stop("Every large portable field must be an ordinary in-object R matrix.")
}

# Reject disk-backed seeds explicitly and reject stale absolute/H5AD paths in
# provenance. Large numeric matrices are treated as terminal leaves.
forbidden_classes <- c(
  "HDF5ArraySeed", "HDF5Array", "DelayedMatrix", "DelayedArray"
)
find_forbidden <- function(x, where = "root") {
  matched <- intersect(class(x), forbidden_classes)
  if (length(matched)) return(paste0(where, "<", matched[[1L]], ">"))
  if (is.matrix(x) || is.data.frame(x) || inherits(x, "DataFrame")) {
    return(character())
  }
  if (is.character(x)) {
    bad <- grepl("(^/|^file:|^[A-Za-z]:[/\\\\])", x, ignore.case = TRUE)
    return(if (any(bad, na.rm = TRUE)) paste0(where, "=", x[which(bad)[[1L]]]) else character())
  }
  if (is.list(x)) {
    labels <- names(x)
    if (is.null(labels)) labels <- as.character(seq_along(x))
    return(unlist(Map(
      function(value, label) find_forbidden(value, paste0(where, "$", label)),
      x, labels
    ), use.names = FALSE))
  }
  if (isS4(x)) {
    slots <- methods::slotNames(x)
    return(unlist(lapply(
      slots,
      function(slot) find_forbidden(methods::slot(x, slot), paste0(where, "@", slot))
    ), use.names = FALSE))
  }
  character()
}
forbidden <- c(
  unlist(lapply(large_fields, find_forbidden), use.names = FALSE),
  find_forbidden(S4Vectors::metadata(spe), "metadata")
)
if (length(forbidden)) {
  stop("Portable SPE retains an external dependency: ", forbidden[[1L]])
}

tiny_info <- S4Vectors::metadata(tiny)$CellCharterR$tiny_subset
if (is.null(tiny_info) || length(tiny_info$source_column_indices) != ncol(tiny)) {
  stop("The tiny SPE lacks its recorded full-object source indices.")
}
indices <- tiny_info$source_column_indices
if (!identical(colnames(tiny), colnames(spe)[indices])) {
  stop("Tiny spot identities do not match the recorded full-object subset.")
}
if (!identical(
  digest::digest(colnames(tiny), algo = "sha256"),
  tiny_info$cell_order_sha256
)) {
  stop("Tiny ordered-spot digest does not match its metadata.")
}
if (!identical(
  SummarizedExperiment::assay(tiny, "protein"),
  large_fields$protein[, indices, drop = FALSE]
)) {
  stop("Tiny protein values differ from the recorded full-object subset.")
}
for (name in expected_reduced_dims) {
  if (!identical(
    SingleCellExperiment::reducedDim(tiny, name),
    large_fields[[name]][indices, , drop = FALSE]
  )) {
    stop("Tiny reduced dimension differs from full object: ", name)
  }
}
if (!identical(
  SpatialExperiment::spatialCoords(tiny),
  large_fields$spatial[indices, , drop = FALSE]
)) {
  stop("Tiny spatial coordinates differ from the full object.")
}
if (!identical(
  SummarizedExperiment::colData(tiny),
  SummarizedExperiment::colData(spe)[indices, , drop = FALSE]
)) {
  stop("Tiny spot metadata differ from the recorded full-object subset.")
}

result <- list(
  format_version = 1L,
  status = "PASS",
  validation_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  verifier = "inst/scripts/verify-codex-portable-spe.R",
  analytical_recomputation = FALSE,
  source_filename = basename(full_rds),
  source_sha256 = source_sha256,
  relocated_sha256 = relocated_sha256,
  serialized_size_bytes = unname(file.info(full_rds)$size),
  object_size_bytes = as.numeric(utils::object.size(spe)),
  dimensions = dim(spe),
  assays = SummarizedExperiment::assayNames(spe),
  reduced_dims = vapply(large_fields[-1L], ncol, integer(1)),
  labels = stats::setNames(
    lapply(required_labels, function(name) {
      value <- SummarizedExperiment::colData(spe)[[name]]
      list(levels = levels(factor(value)), sha256 = digest::digest(value, algo = "sha256"))
    }),
    required_labels
  ),
  field_sha256 = vapply(
    large_fields, digest::digest, character(1), algo = "sha256"
  ),
  external_seed_or_path_matches = forbidden,
  tiny = list(
    spots = ncol(tiny),
    cell_order_sha256 = tiny_info$cell_order_sha256,
    exact_assay_match = TRUE,
    exact_reduced_dim_match = TRUE,
    exact_coldata_match = TRUE
  )
)
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_rds, version = 3L)
cat("Portable SPE validation: PASS\n")
cat("Full SHA-256:", source_sha256, "\n")
cat("Machine-readable result:", output_rds, "\n")
