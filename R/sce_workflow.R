# Utilities in this file keep the full workflow centred on a
# SingleCellExperiment while its large numeric matrices remain in HDF5.

.cc_h5ad_dataframe <- function(filepath, group) {
  listing <- rhdf5::h5ls(filepath, recursive = TRUE)
  direct <- listing[listing$group == paste0("/", group), , drop = FALSE]
  if (!nrow(direct)) stop("H5AD group `", group, "` was not found in ", filepath, ".")

  index_path <- paste0(group, "/_index")
  index <- as.character(rhdf5::h5read(filepath, index_path))
  columns <- setdiff(direct$name, "_index")
  values <- stats::setNames(vector("list", length(columns)), columns)

  for (column in columns) {
    row <- direct[direct$name == column, , drop = FALSE]
    path <- paste0(group, "/", column)
    if (row$otype[[1L]] == "H5I_GROUP") {
      categories <- as.character(rhdf5::h5read(filepath, paste0(path, "/categories")))
      codes <- as.integer(rhdf5::h5read(filepath, paste0(path, "/codes")))
      decoded <- rep.int(NA_character_, length(codes))
      keep <- codes >= 0L
      decoded[keep] <- categories[codes[keep] + 1L]
      values[[column]] <- factor(decoded, levels = categories)
    } else {
      values[[column]] <- as.vector(rhdf5::h5read(filepath, path))
    }
  }

  frame <- S4Vectors::DataFrame(values, check.names = FALSE)
  rownames(frame) <- index
  frame
}

.cc_h5ad_index <- function(filepath, group = "obs") {
  as.character(rhdf5::h5read(filepath, paste0(group, "/_index")))
}

#' Compute a SHA-256 checksum without loading a file into R memory
#'
#' @param filepath Existing file.
#' @return Lower-case hexadecimal SHA-256 digest.
#' @export
cellcharter_file_sha256 <- function(filepath) {
  if (!file.exists(filepath)) stop("File does not exist: ", filepath)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Install `digest` to record SHA-256 provenance.")
  }
  digest::digest(file = filepath, algo = "sha256", serialize = FALSE)
}

#' Read the fixed CODEX inputs as a disk-backed SingleCellExperiment
#'
#' The marker assay is read from the historical tutorial H5AD and the immutable
#' trVAE representation from the official clustered H5AD. Cell identifiers must
#' match exactly before either object is attached. The resulting SCE contains
#' HDF5-backed assay and reduced-dimension matrices rather than dense copies.
#'
#' @param raw_h5ad Official CODEX input H5AD (Figshare file 38538101).
#' @param clustered_h5ad Official clustered H5AD containing `X_trVAE`
#'   (Figshare file 46832722).
#' @param representation_dataset HDF5 dataset holding the fixed representation.
#' @param spatial_dataset HDF5 dataset holding x/y coordinates.
#' @return A `SingleCellExperiment` with disk-backed `protein`, `X_trVAE`, and
#'   `spatial` matrices plus input provenance in `metadata()`.
#' @export
read_codex_sce <- function(
    raw_h5ad, clustered_h5ad,
    representation_dataset = "obsm/X_trVAE",
    spatial_dataset = "obsm/spatial") {
  for (path in c(raw_h5ad, clustered_h5ad)) {
    if (!file.exists(path)) stop("Required CODEX input is absent: ", path)
  }
  raw_h5ad <- normalizePath(raw_h5ad)
  clustered_h5ad <- normalizePath(clustered_h5ad)
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Install `HDF5Array` to read CODEX without materialising the assay.")
  }

  raw_ids <- .cc_h5ad_index(raw_h5ad)
  clustered_ids <- .cc_h5ad_index(clustered_h5ad)
  ids_exact <- identical(raw_ids, clustered_ids)
  # The later official artifact resets each sample's local numeric index,
  # whereas the historical H5AD made repeated indices unique with `-<batch>`.
  # The order is accepted only when removing that suffix yields an exact
  # one-for-one match; the SCE itself keeps the historical globally unique IDs.
  ids_local_match <- identical(sub("-[0-9]+$", "", raw_ids), clustered_ids)
  if (!ids_exact && !ids_local_match) {
    stop("Raw and clustered CODEX H5AD cell identifiers/order do not match, even after documented sample suffix normalization.")
  }

  marker_assay <- HDF5Array::HDF5Array(raw_h5ad, "X")
  representation <- DelayedArray::t(HDF5Array::HDF5Array(
    clustered_h5ad, representation_dataset
  ))
  spatial <- DelayedArray::t(HDF5Array::HDF5Array(raw_h5ad, spatial_dataset))
  if (ncol(marker_assay) != length(raw_ids) ||
      nrow(representation) != length(raw_ids) || nrow(spatial) != length(raw_ids)) {
    stop("CODEX HDF5 matrix dimensions are inconsistent with cell metadata.")
  }

  col_data <- .cc_h5ad_dataframe(raw_h5ad, "obs")
  clustered_col_data <- .cc_h5ad_dataframe(clustered_h5ad, "obs")
  if (nrow(clustered_col_data) != nrow(col_data) ||
      !"cluster_cellcharter" %in% colnames(clustered_col_data)) {
    stop("Clustered CODEX metadata lacks the published CellCharter labels.")
  }
  # Preserve the published K=11 labels only as a sensitivity/reference column;
  # they are never used to select K or fit the primary R model.
  col_data$reference_cluster_cellcharter <- clustered_col_data$cluster_cellcharter
  if (!"condition" %in% colnames(col_data)) {
    col_data$condition <- factor(sub("-.*$", "", as.character(col_data$sample)))
  }
  marker_names <- .cc_h5ad_index(raw_h5ad, "var")
  row_data <- S4Vectors::DataFrame(marker = marker_names, row.names = marker_names)
  rownames(marker_assay) <- marker_names
  colnames(marker_assay) <- raw_ids
  rownames(representation) <- raw_ids
  rownames(spatial) <- raw_ids
  colnames(representation) <- paste0("trVAE_", seq_len(ncol(representation)))
  colnames(spatial) <- c("x", "y")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(protein = marker_assay), rowData = row_data,
    colData = col_data,
    reducedDims = S4Vectors::SimpleList(X_trVAE = representation, spatial = spatial)
  )
  S4Vectors::metadata(sce)$rCellCharter <- list(
    reference = list(
      cellcharter_version = "0.2.0",
      cellcharter_commit = "2a8fdde8213e32e5a00f1f9571bffd3dbb825a79",
      squidpy_version = "1.3.0",
      squidpy_commit = "e2ffc4d03d02e1b02c9bf15ae87fb7091a64df4b",
      python_version = "3.10",
      pytorch_version = "1.12.1",
      pycave_version = "3.2.1",
      scikit_learn_version = "1.3.0",
      scipy_version = "not declared by cellcharter_analyses environment",
      numpy_version = "1.23.4",
      scarches_version = "0.5.9"
    ),
    inputs = list(
      raw_h5ad = normalizePath(raw_h5ad),
      raw_sha256 = cellcharter_file_sha256(raw_h5ad),
      clustered_h5ad = normalizePath(clustered_h5ad),
      clustered_sha256 = cellcharter_file_sha256(clustered_h5ad),
      representation_dataset = representation_dataset,
      cell_order_sha256 = digest::digest(raw_ids, algo = "sha256"),
      clustered_cell_order_sha256 = digest::digest(clustered_ids, algo = "sha256"),
      cell_id_alignment = if (ids_exact) "exact" else "raw suffix normalization"
    )
  )
  sce
}

#' Attach an HDF5-backed CellCharter representation to an SCE
#'
#' @param sce A `SingleCellExperiment`.
#' @param filepath HDF5 aggregation output.
#' @param name Dataset name in `filepath`.
#' @param reduced_dim_name Name to use in `reducedDims(sce)`.
#' @return The updated SCE.
#' @export
cellcharter_attach_hdf5_reduced_dim <- function(
    sce, filepath, name = "X_cellcharter", reduced_dim_name = name) {
  if (!inherits(sce, "SingleCellExperiment")) stop("`sce` must be a SingleCellExperiment.")
  value <- HDF5Array::HDF5Array(filepath, name)
  if (nrow(value) != ncol(sce)) stop("Reduced dimension rows must match SCE cells.")
  rownames(value) <- colnames(sce)
  SingleCellExperiment::reducedDim(sce, reduced_dim_name) <- value
  sce
}
