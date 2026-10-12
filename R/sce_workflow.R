# Utilities in this file keep the full workflow centred on a
# SingleCellExperiment-compatible container. Spatial inputs use the
# SpatialExperiment subclass while large numeric matrices remain in HDF5.

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

.cc_validate_observation_ids <- function(
    ids, context = "Observation identifiers", allow_duplicates = FALSE) {
  if (is.null(ids)) {
    stop(context, " are missing; positional alignment is not accepted.", call. = FALSE)
  }
  ids <- as.character(ids)
  if (!length(ids)) stop(context, " are empty.", call. = FALSE)
  if (anyNA(ids)) stop(context, " contain missing values.", call. = FALSE)
  blank <- !nzchar(trimws(ids))
  if (any(blank)) stop(context, " contain blank values.", call. = FALSE)
  if (!allow_duplicates && anyDuplicated(ids)) {
    duplicate <- ids[duplicated(ids)][[1L]]
    stop(context, " contain duplicated identifier `", duplicate, "`.", call. = FALSE)
  }
  ids
}

.cc_observation_order_sha256 <- function(ids, context = "Observation identifiers") {
  ids <- .cc_validate_observation_ids(ids, context)
  digest::digest(ids, algo = "sha256")
}

.cc_align_observation_ids <- function(
    source_ids, target_ids, source_context = "Source observation identifiers",
    target_context = "Target observation identifiers", allow_reorder = FALSE) {
  source_ids <- .cc_validate_observation_ids(source_ids, source_context)
  target_ids <- .cc_validate_observation_ids(target_ids, target_context)
  if (length(source_ids) != length(target_ids) ||
      !setequal(source_ids, target_ids)) {
    missing <- setdiff(target_ids, source_ids)
    unexpected <- setdiff(source_ids, target_ids)
    detail <- c(
      if (length(missing)) paste0("missing `", missing[[1L]], "`"),
      if (length(unexpected)) paste0("unexpected `", unexpected[[1L]], "`")
    )
    stop(
      source_context, " and ", target_context, " have an identifier set mismatch",
      if (length(detail)) paste0(" (", paste(detail, collapse = "; "), ")") else "",
      ".", call. = FALSE
    )
  }
  if (identical(source_ids, target_ids)) {
    return(list(index = seq_along(target_ids), reordered = FALSE))
  }
  if (!isTRUE(allow_reorder)) {
    first <- which(source_ids != target_ids)[[1L]]
    stop(
      source_context, " and ", target_context,
      " have an order mismatch at position ", first, " (`", source_ids[[first]],
      "` versus `", target_ids[[first]], "`).", call. = FALSE
    )
  }
  list(index = match(target_ids, source_ids), reordered = TRUE)
}

.cc_validate_sce_identity <- function(
    sce, context = "SingleCellExperiment checkpoint", reduced_dims = NULL,
    recorded_order_sha256 = NULL) {
  if (!inherits(sce, "SingleCellExperiment")) {
    stop(context, " is not a SingleCellExperiment.", call. = FALSE)
  }
  ids <- .cc_validate_observation_ids(colnames(sce), paste0(context, " cell identifiers"))
  .cc_align_observation_ids(
    rownames(SummarizedExperiment::colData(sce)), ids,
    paste0(context, " colData row names"), paste0(context, " cell identifiers")
  )
  if (is.null(reduced_dims)) reduced_dims <- SingleCellExperiment::reducedDimNames(sce)
  for (reduced_dim_name in reduced_dims) {
    if (!reduced_dim_name %in% SingleCellExperiment::reducedDimNames(sce)) {
      stop(context, " lacks reduced dimension `", reduced_dim_name, "`.", call. = FALSE)
    }
    .cc_align_observation_ids(
      rownames(SingleCellExperiment::reducedDim(sce, reduced_dim_name)), ids,
      paste0(context, " `", reduced_dim_name, "` row names"),
      paste0(context, " cell identifiers")
    )
  }
  if (inherits(sce, "SpatialExperiment")) {
    spatial <- SpatialExperiment::spatialCoords(sce)
    if (!is.numeric(spatial) || length(dim(spatial)) != 2L ||
        nrow(spatial) != ncol(sce) || ncol(spatial) < 2L) {
      stop(context, " has invalid SpatialExperiment coordinates.", call. = FALSE)
    }
    .cc_align_observation_ids(
      rownames(spatial), ids,
      paste0(context, " spatialCoords row names"),
      paste0(context, " cell identifiers")
    )
  }
  order_sha256 <- .cc_observation_order_sha256(ids, paste0(context, " cell identifiers"))
  if (!is.null(recorded_order_sha256) &&
      (!is.character(recorded_order_sha256) || length(recorded_order_sha256) != 1L ||
       is.na(recorded_order_sha256) || !identical(recorded_order_sha256, order_sha256))) {
    stop(context, " cell order does not match its recorded SHA-256 identity.", call. = FALSE)
  }
  invisible(list(ids = ids, cell_order_sha256 = order_sha256))
}

.cc_codex_h5ad_identity <- function(
    raw_h5ad, clustered_h5ad, raw_col_data = NULL, clustered_col_data = NULL) {
  raw_ids <- .cc_validate_observation_ids(
    .cc_h5ad_index(raw_h5ad), "Raw CODEX H5AD observation identifiers"
  )
  clustered_ids <- .cc_validate_observation_ids(
    .cc_h5ad_index(clustered_h5ad),
    "Clustered CODEX H5AD local observation identifiers", allow_duplicates = TRUE
  )
  if (is.null(raw_col_data)) raw_col_data <- .cc_h5ad_dataframe(raw_h5ad, "obs")
  if (is.null(clustered_col_data)) {
    clustered_col_data <- .cc_h5ad_dataframe(clustered_h5ad, "obs")
  }
  if (!"sample" %in% colnames(raw_col_data) ||
      !"sample" %in% colnames(clustered_col_data)) {
    stop("Both CODEX H5AD inputs must provide `obs/sample` for immutable identity alignment.",
         call. = FALSE)
  }
  raw_sample <- .cc_validate_observation_ids(
    as.character(raw_col_data$sample), "Raw CODEX H5AD sample identifiers",
    allow_duplicates = TRUE
  )
  clustered_sample <- .cc_validate_observation_ids(
    as.character(clustered_col_data$sample), "Clustered CODEX H5AD sample identifiers",
    allow_duplicates = TRUE
  )
  if (length(raw_ids) != length(raw_sample) ||
      length(clustered_ids) != length(clustered_sample)) {
    stop("CODEX H5AD observation identifiers and sample metadata have different lengths.",
         call. = FALSE)
  }

  ids_exact <- !anyDuplicated(clustered_ids) && identical(raw_ids, clustered_ids)
  if (ids_exact) {
    if (!identical(raw_sample, clustered_sample)) {
      mismatch <- which(raw_sample != clustered_sample)
      position <- if (length(mismatch)) mismatch[[1L]] else NA_integer_
      if (is.na(position)) {
        stop("Raw and clustered CODEX H5AD sample metadata have an order mismatch.",
             call. = FALSE)
      }
      stop("Raw and clustered CODEX H5AD sample metadata have an order mismatch at position ",
           position, ".", call. = FALSE)
    }
    canonical_raw <- raw_ids
    canonical_clustered <- clustered_ids
    alignment <- "exact global observation identifiers"
  } else {
    # The official clustered artifact uses sample-local indices. Qualifying
    # them by `sample` makes that documented identity scheme unambiguous.
    canonical_raw <- paste(raw_sample, sub("-[0-9]+$", "", raw_ids), sep = "\037")
    canonical_clustered <- paste(clustered_sample, clustered_ids, sep = "\037")
    canonical_raw <- .cc_validate_observation_ids(
      canonical_raw, "Raw CODEX H5AD sample-scoped observation identifiers"
    )
    canonical_clustered <- .cc_validate_observation_ids(
      canonical_clustered, "Clustered CODEX H5AD sample-scoped observation identifiers"
    )
    .cc_align_observation_ids(
      canonical_clustered, canonical_raw,
      "Clustered CODEX H5AD sample-scoped observation identifiers",
      "raw CODEX H5AD sample-scoped observation identifiers"
    )
    alignment <- "exact sample-scoped local observation identifiers"
  }
  list(
    raw_ids = raw_ids, clustered_ids = clustered_ids,
    canonical_raw_ids = canonical_raw,
    canonical_clustered_ids = canonical_clustered,
    raw_sample = raw_sample, clustered_sample = clustered_sample,
    alignment = alignment,
    cell_order_sha256 = .cc_observation_order_sha256(raw_ids),
    clustered_cell_order_sha256 = digest::digest(clustered_ids, algo = "sha256"),
    canonical_cell_order_sha256 = .cc_observation_order_sha256(canonical_raw)
  )
}

#' Compute a SHA-256 checksum without loading a file into R memory
#'
#' @param filepath Existing file.
#' @return Lower-case hexadecimal SHA-256 digest.
#' @examples
#' example_file <- tempfile(fileext = ".txt")
#' writeLines("CellCharterR", example_file)
#' cellcharter_file_sha256(example_file)
#' unlink(example_file)
#' @export
cellcharter_file_sha256 <- function(filepath) {
  if (!file.exists(filepath)) stop("File does not exist: ", filepath)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Install `digest` to record SHA-256 provenance.")
  }
  digest::digest(file = filepath, algo = "sha256", serialize = FALSE)
}

#' Read the fixed CODEX inputs as a disk-backed SpatialExperiment
#'
#' The marker assay is read from the historical tutorial H5AD and the immutable
#' trVAE representation from the official clustered H5AD. Cell identifiers must
#' match exactly before either object is attached. The resulting
#' `SpatialExperiment` contains an HDF5-backed protein assay and trVAE reduced
#' dimension. Only the two-column coordinate matrix is realized, because
#' `spatialCoords()` is the canonical SpatialExperiment coordinate store.
#'
#' @param raw_h5ad Official CODEX input H5AD (Figshare file 38538101).
#' @param clustered_h5ad Official clustered H5AD containing `X_trVAE`
#'   (Figshare file 46832722).
#' @param representation_dataset HDF5 dataset holding the fixed representation.
#' @param spatial_dataset HDF5 dataset holding x/y coordinates.
#' @return A `SpatialExperiment` with disk-backed `protein` and `X_trVAE`,
#'   canonical `spatialCoords()`, and input provenance in `metadata()`.
#' @examples
#' raw <- system.file("extdata", "codex_mouse_spleen_tiny_raw.h5ad",
#'                    package = "CellCharterR")
#' clustered <- system.file("extdata", "codex_mouse_spleen_tiny_clustered.h5ad",
#'                          package = "CellCharterR")
#' example_spe <- read_codex_spe(raw, clustered)
#' example_spe
#' @export
read_codex_spe <- function(
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

  col_data <- .cc_h5ad_dataframe(raw_h5ad, "obs")
  clustered_col_data <- .cc_h5ad_dataframe(clustered_h5ad, "obs")
  identity <- .cc_codex_h5ad_identity(
    raw_h5ad, clustered_h5ad, col_data, clustered_col_data
  )
  raw_ids <- identity$raw_ids
  clustered_ids <- identity$clustered_ids

  marker_assay <- HDF5Array::HDF5Array(raw_h5ad, "X")
  representation <- DelayedArray::t(HDF5Array::HDF5Array(
    clustered_h5ad, representation_dataset
  ))
  spatial <- DelayedArray::t(HDF5Array::HDF5Array(raw_h5ad, spatial_dataset))
  if (ncol(marker_assay) != length(raw_ids) ||
      nrow(representation) != length(raw_ids) || nrow(spatial) != length(raw_ids)) {
    stop("CODEX HDF5 matrix dimensions are inconsistent with cell metadata.")
  }

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
  marker_names <- .cc_validate_observation_ids(
    .cc_h5ad_index(raw_h5ad, "var"), "Raw CODEX H5AD variable identifiers"
  )
  row_data <- S4Vectors::DataFrame(marker = marker_names, row.names = marker_names)
  rownames(marker_assay) <- marker_names
  colnames(marker_assay) <- raw_ids
  rownames(representation) <- raw_ids
  rownames(spatial) <- raw_ids
  colnames(representation) <- paste0("trVAE_", seq_len(ncol(representation)))
  colnames(spatial) <- c("x", "y")

  # SpatialExperiment requires an ordinary numeric coordinate matrix. This
  # realizes only two columns, never the protein assay or trVAE representation.
  spatial <- as.matrix(spatial)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(protein = marker_assay), rowData = row_data,
    colData = col_data,
    reducedDims = S4Vectors::SimpleList(X_trVAE = representation),
    spatialCoords = spatial,
    sample_id = as.character(col_data$sample)
  )
  S4Vectors::metadata(spe)$CellCharterR <- list(
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
      spatial_dataset = spatial_dataset,
      cell_order_sha256 = identity$cell_order_sha256,
      clustered_cell_order_sha256 = identity$clustered_cell_order_sha256,
      canonical_cell_order_sha256 = identity$canonical_cell_order_sha256,
      cell_id_alignment = identity$alignment
    )
  )
  spe
}

#' @rdname read_codex_spe
#' @description `read_codex_sce()` is retained as a backward-compatible alias;
#'   it returns the same `SpatialExperiment` as `read_codex_spe()`.
#' @export
read_codex_sce <- function(
    raw_h5ad, clustered_h5ad,
    representation_dataset = "obsm/X_trVAE",
    spatial_dataset = "obsm/spatial") {
  read_codex_spe(
    raw_h5ad = raw_h5ad,
    clustered_h5ad = clustered_h5ad,
    representation_dataset = representation_dataset,
    spatial_dataset = spatial_dataset
  )
}

#' Attach an HDF5-backed CellCharter representation to an SCE
#'
#' @param sce A `SingleCellExperiment`.
#' @param filepath HDF5 aggregation output.
#' @param name Dataset name in `filepath`.
#' @param reduced_dim_name Name to use in `reducedDims(sce)`.
#' @param cell_ids Observation identifiers for the HDF5 rows. When omitted,
#'   row names on the HDF5 dataset or a `<name>_cell_ids` HDF5 dataset are used.
#' @param allow_reorder Whether an explicitly identified HDF5 matrix may be
#'   reordered into SCE order. The reorder is recorded in SCE metadata.
#' @return The updated SCE.
#' @examples
#' ids <- paste0("cell", seq_len(3))
#' example_sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = matrix(0, 1, 3, dimnames = list("marker", ids)))
#' )
#' output <- tempfile(fileext = ".h5")
#' HDF5Array::writeHDF5Array(matrix(seq_len(6), nrow = 3),
#'                           filepath = output, name = "X_cellcharter")
#' example_sce <- cellcharter_attach_hdf5_reduced_dim(
#'   example_sce, output, cell_ids = ids
#' )
#' SingleCellExperiment::reducedDimNames(example_sce)
#' unlink(output)
#' @export
cellcharter_attach_hdf5_reduced_dim <- function(
    sce, filepath, name = "X_cellcharter", reduced_dim_name = name,
    cell_ids = NULL, allow_reorder = FALSE) {
  if (!inherits(sce, "SingleCellExperiment")) stop("`sce` must be a SingleCellExperiment.")
  target <- .cc_validate_observation_ids(colnames(sce), "SCE cell identifiers")
  value <- HDF5Array::HDF5Array(filepath, name)
  if (is.null(cell_ids)) cell_ids <- rownames(value)
  if (is.null(cell_ids)) {
    cell_ids <- try(
      as.character(rhdf5::h5read(filepath, paste0(name, "_cell_ids"))), silent = TRUE
    )
    if (inherits(cell_ids, "try-error")) cell_ids <- NULL
  }
  source <- .cc_validate_observation_ids(
    cell_ids, paste0("HDF5 dataset `", name, "` observation identifiers")
  )
  if (nrow(value) != length(source)) {
    stop("HDF5 reduced-dimension rows do not match its observation identifiers.",
         call. = FALSE)
  }
  alignment <- .cc_align_observation_ids(
    source, target, paste0("HDF5 dataset `", name, "` observation identifiers"),
    "SCE cell identifiers", allow_reorder = allow_reorder
  )
  if (alignment$reordered) value <- value[alignment$index, , drop = FALSE]
  rownames(value) <- target
  SingleCellExperiment::reducedDim(sce, reduced_dim_name) <- value
  if (is.null(S4Vectors::metadata(sce)$CellCharterR)) {
    S4Vectors::metadata(sce)$CellCharterR <- list()
  }
  S4Vectors::metadata(sce)$CellCharterR$identity_alignment[[reduced_dim_name]] <- list(
    source = normalizePath(filepath), dataset = name,
    source_cell_order_sha256 = .cc_observation_order_sha256(source),
    target_cell_order_sha256 = .cc_observation_order_sha256(target),
    reordered = alignment$reordered
  )
  sce
}
