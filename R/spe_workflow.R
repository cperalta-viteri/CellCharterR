# Utilities in this file keep the full workflow centred on a
# SpatialExperiment container while large numeric matrices remain in HDF5.

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

.cc_observation_order_sha256 <- function(
    ids, context = "Observation identifiers", allow_duplicates = FALSE) {
  ids <- .cc_validate_observation_ids(
    ids, context, allow_duplicates = allow_duplicates
  )
  digest::digest(ids, algo = "sha256")
}

.cc_align_observation_ids <- function(
    source_ids, target_ids, source_context = "Source observation identifiers",
    target_context = "Target observation identifiers", allow_reorder = FALSE,
    allow_duplicates = FALSE) {
  source_ids <- .cc_validate_observation_ids(
    source_ids, source_context, allow_duplicates = allow_duplicates
  )
  target_ids <- .cc_validate_observation_ids(
    target_ids, target_context, allow_duplicates = allow_duplicates
  )
  if (identical(source_ids, target_ids)) {
    return(list(index = seq_along(target_ids), reordered = FALSE))
  }
  if (allow_duplicates &&
      (anyDuplicated(source_ids) || anyDuplicated(target_ids))) {
    first <- if (length(source_ids) == length(target_ids)) {
      which(source_ids != target_ids)[[1L]]
    } else {
      min(length(source_ids), length(target_ids)) + 1L
    }
    stop(
      source_context, " and ", target_context,
      " have an order mismatch at position ", first,
      "; duplicated identifiers require exact positional agreement.",
      call. = FALSE
    )
  }
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

.cc_legacy_sce_as_spe <- function(x, spatial_reduced_dim = "spatial") {
  if (inherits(x, "SpatialExperiment")) return(x)
  if (!inherits(x, "SingleCellExperiment")) {
    stop("Legacy object is not a SingleCellExperiment.", call. = FALSE)
  }
  if (!spatial_reduced_dim %in%
      SingleCellExperiment::reducedDimNames(x)) {
    stop(
      "Legacy SingleCellExperiment lacks reducedDim `", spatial_reduced_dim,
      "` required to reconstruct spatialCoords.", call. = FALSE
    )
  }
  spatial <- as.matrix(SingleCellExperiment::reducedDim(
    x, spatial_reduced_dim, withDimnames = TRUE
  ))
  if (!is.numeric(spatial) || nrow(spatial) != ncol(x) || ncol(spatial) < 2L ||
      any(!is.finite(spatial[, seq_len(2L), drop = FALSE]))) {
    stop("Legacy spatial reduced dimension must contain finite x/y coordinates.",
         call. = FALSE)
  }
  spatial <- spatial[, seq_len(2L), drop = FALSE]
  rownames(spatial) <- colnames(x)
  if (is.null(colnames(spatial))) colnames(spatial) <- c("x", "y")

  # Rebuild only the container shell. HDF5-backed assays and reduced dimensions
  # remain lazy seeds; the small two-column coordinate matrix is the only value
  # intentionally materialised for SpatialExperiment's canonical store.
  SpatialExperiment::SpatialExperiment(
    assays = SummarizedExperiment::assays(x, withDimnames = FALSE),
    rowData = SummarizedExperiment::rowData(x),
    colData = SummarizedExperiment::colData(x),
    reducedDims = SingleCellExperiment::reducedDims(x),
    altExps = SingleCellExperiment::altExps(x),
    rowPairs = SingleCellExperiment::rowPairs(x),
    colPairs = SingleCellExperiment::colPairs(x),
    metadata = S4Vectors::metadata(x),
    spatialCoords = spatial,
    mainExpName = SingleCellExperiment::mainExpName(x)
  )
}

.cc_validate_spe_identity <- function(
    spe, context = "SpatialExperiment checkpoint", reduced_dims = NULL,
    recorded_order_sha256 = NULL) {
  if (!inherits(spe, "SpatialExperiment")) {
    stop(context, " is not a SpatialExperiment.", call. = FALSE)
  }
  n_spots <- ncol(spe)
  column_ids <- colnames(spe)
  if (!is.null(column_ids)) {
    column_ids <- .cc_validate_observation_ids(
      column_ids, paste0(context, " spot identifiers"),
      allow_duplicates = TRUE
    )
  }

  # SPE/SCE components are positionally aligned by the container contract.
  # Their optional dimnames are annotations, not independent foreign keys, so
  # a missing or legacy component rowname must never make valid current data
  # unusable. Prefer a genuinely unique identifier for provenance hashes when
  # one is available; otherwise use an explicitly positional fallback and let
  # content fingerprints protect checkpoint and result reuse.
  col_data <- SummarizedExperiment::colData(spe)
  candidates <- list(
    colnames = column_ids,
    colData_rownames = rownames(col_data)
  )
  if (!is.null(column_ids) && "sample_id" %in% colnames(col_data)) {
    sample_ids <- as.character(col_data$sample_id)
    if (length(sample_ids) == n_spots && !anyNA(sample_ids) &&
        all(nzchar(sample_ids))) {
      candidates$sample_scoped_colnames <- paste(
        sample_ids, column_ids, sep = "\037"
      )
    }
  }
  unique_candidate <- vapply(candidates, function(candidate) {
    !is.null(candidate) && length(candidate) == n_spots &&
      !anyNA(candidate) && all(nzchar(as.character(candidate))) &&
      !anyDuplicated(candidate)
  }, logical(1L))
  if (any(unique_candidate)) {
    identity_source <- names(candidates)[which(unique_candidate)[[1L]]]
    ids <- as.character(candidates[[identity_source]])
  } else {
    identity_source <- "position"
    labels <- if (is.null(column_ids)) rep.int("", n_spots) else column_ids
    ids <- paste0("position:", seq_len(n_spots), "\037", labels)
  }

  # NULL means that no reduced dimension is consumed by this operation. This
  # prevents an unrelated legacy reducedDim name from blocking a colData-only
  # analysis. Callers that consume a representation pass its name explicitly.
  if (is.null(reduced_dims)) reduced_dims <- character()
  for (reduced_dim_name in reduced_dims) {
    if (!reduced_dim_name %in% SingleCellExperiment::reducedDimNames(spe)) {
      stop(context, " lacks reduced dimension `", reduced_dim_name, "`.", call. = FALSE)
    }
    reduced <- SingleCellExperiment::reducedDim(
      spe, reduced_dim_name, withDimnames = FALSE
    )
    if (length(dim(reduced)) != 2L || nrow(reduced) != n_spots) {
      stop(context, " reduced dimension `", reduced_dim_name,
           "` is not positionally aligned with its spots.", call. = FALSE)
    }
  }
  if (inherits(spe, "SpatialExperiment")) {
    spatial <- SpatialExperiment::spatialCoords(spe)
    if (!is.numeric(spatial) || length(dim(spatial)) != 2L ||
        nrow(spatial) != ncol(spe) || ncol(spatial) < 2L) {
      stop(context, " has invalid SpatialExperiment coordinates.", call. = FALSE)
    }
    # SpatialExperiment keeps coordinates positionally aligned with columns
    # when the object is constructed and subset. Coordinate row names are
    # optional annotations, not the container's canonical column identifiers;
    # imported objects may therefore retain unique H5AD barcodes here while
    # using repeated slide-level labels in colnames(spe). Requiring the two
    # naming systems to match rejects a valid, internally aligned SPE.
  }
  order_sha256 <- .cc_observation_order_sha256(
    ids, paste0(context, " spot identifiers"), allow_duplicates = TRUE
  )
  if (!is.null(recorded_order_sha256) &&
      (!is.character(recorded_order_sha256) || length(recorded_order_sha256) != 1L ||
       is.na(recorded_order_sha256) || !identical(recorded_order_sha256, order_sha256))) {
    stop(context, " spot order does not match its recorded SHA-256 identity.", call. = FALSE)
  }
  invisible(list(
    ids = ids, cell_order_sha256 = order_sha256,
    identity_source = identity_source,
    identity_is_unique = !identical(identity_source, "position")
  ))
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
#' @param filepath One existing readable file path. The file is streamed by
#'   [digest::digest()] and is not loaded into R memory.
#' @return One 64-character lower-case hexadecimal SHA-256 digest.
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
#' trVAE representation from the official clustered H5AD. Spot identifiers must
#' match exactly before either object is attached. The resulting
#' `SpatialExperiment` contains an HDF5-backed protein assay and trVAE reduced
#' dimension. Only the two-column coordinate matrix is realized, because
#' `spatialCoords()` is the canonical SpatialExperiment coordinate store.
#'
#' @param raw_h5ad Path to the official CODEX input H5AD (Figshare file
#'   38538101). Its dense `X` dataset becomes the HDF5-backed `protein` assay;
#'   its `obs` and `var` tables populate `colData` and `rowData`.
#' @param clustered_h5ad Path to the official clustered H5AD (Figshare file
#'   46832722). Its published cluster labels are stored in
#'   `colData(spe)$reference_cluster_cellcharter` after immutable spot-identity
#'   alignment against `raw_h5ad`.
#' @param representation_dataset Non-empty HDF5 path in `clustered_h5ad` to the
#'   fixed trVAE spot-by-feature representation. It is attached disk-backed as
#'   `reducedDim(spe, "X_trVAE")` and is not materialised.
#' @param spatial_dataset Non-empty HDF5 path in `raw_h5ad` to the two-column
#'   spatial coordinate matrix. Only these coordinates are materialised and
#'   stored in `spatialCoords(spe)` with columns `x` and `y`.
#' @return A `SpatialExperiment` with disk-backed `protein` and `X_trVAE`,
#'   canonical `spatialCoords()`, published labels in `colData`, and file paths,
#'   checksums, pinned reference versions, and identity-alignment provenance in
#'   `metadata(spe)$CellCharterR`.
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
    stop("CODEX HDF5 matrix dimensions are inconsistent with spot metadata.")
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

#' Attach an HDF5-backed representation to a SpatialExperiment
#'
#' The HDF5 dataset is attached without materialising its values. Its rows must
#' identify the same spots as `spe`; order differences are rejected unless
#' `allow_reorder = TRUE`. The returned object therefore remains dependent on
#' `filepath` and the file must be moved together with the object.
#'
#' @param spe A `SpatialExperiment` whose columns are spatial spots.
#' @param filepath Path to an existing HDF5 file containing the representation.
#' @param name Non-empty name of the two-dimensional HDF5 dataset in
#'   `filepath`.
#' @param reduced_dim_name Non-empty name under which the disk-backed matrix is
#'   stored in `reducedDim(spe, reduced_dim_name)`. It defaults to `name` and
#'   replaces a reduced dimension with the same name.
#' @param cell_ids Character vector identifying the HDF5 rows. When omitted,
#'   row names on the HDF5 dataset or a `<name>_cell_ids` dataset in the same
#'   file are used. Positional alignment without identifiers is rejected.
#' @param allow_reorder A single logical value. If `TRUE`, rows with the same
#'   unique identifiers may be reordered into `colnames(spe)` order. If
#'   `FALSE`, their order must already match. The alignment decision and both
#'   order checksums are stored in
#'   `metadata(spe)$CellCharterR$identity_alignment[[reduced_dim_name]]`.
#' @return The updated `SpatialExperiment`, with an `HDF5Array` stored in
#'   `reducedDim(spe, reduced_dim_name)` and alignment provenance stored in
#'   `metadata(spe)$CellCharterR$identity_alignment`.
#' @examples
#' ids <- paste0("spot", seq_len(3))
#' example_spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = matrix(0, 1, 3, dimnames = list("marker", ids))),
#'   spatialCoords = matrix(seq_len(6), ncol = 2,
#'                          dimnames = list(ids, c("x", "y")))
#' )
#' output <- tempfile(fileext = ".h5")
#' HDF5Array::writeHDF5Array(matrix(seq_len(6), nrow = 3),
#'                           filepath = output, name = "X_cellcharter")
#' example_spe <- cellcharter_attach_hdf5_reduced_dim(
#'   spe = example_spe, filepath = output, cell_ids = ids
#' )
#' SingleCellExperiment::reducedDimNames(example_spe)
#' unlink(output)
#' @export
cellcharter_attach_hdf5_reduced_dim <- function(
    spe, filepath, name = "X_cellcharter", reduced_dim_name = name,
    cell_ids = NULL, allow_reorder = FALSE) {
  if (!inherits(spe, "SpatialExperiment")) {
    stop("`spe` must be a SpatialExperiment.", call. = FALSE)
  }
  if (!is.character(filepath) || length(filepath) != 1L || is.na(filepath) ||
      !nzchar(filepath) || !file.exists(filepath)) {
    stop("`filepath` must name one existing HDF5 file.", call. = FALSE)
  }
  names_to_validate <- list(name = name, reduced_dim_name = reduced_dim_name)
  for (argument in names(names_to_validate)) {
    value_name <- names_to_validate[[argument]]
    if (!is.character(value_name) || length(value_name) != 1L ||
        is.na(value_name) || !nzchar(value_name)) {
      stop("`", argument, "` must be one non-empty name.", call. = FALSE)
    }
  }
  if (!is.logical(allow_reorder) || length(allow_reorder) != 1L ||
      is.na(allow_reorder)) {
    stop("`allow_reorder` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  target <- .cc_validate_observation_ids(colnames(spe), "SPE spot identifiers")
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
    "SPE spot identifiers", allow_reorder = allow_reorder
  )
  if (alignment$reordered) value <- value[alignment$index, , drop = FALSE]
  rownames(value) <- target
  SingleCellExperiment::reducedDim(spe, reduced_dim_name) <- value
  if (is.null(S4Vectors::metadata(spe)$CellCharterR)) {
    S4Vectors::metadata(spe)$CellCharterR <- list()
  }
  S4Vectors::metadata(spe)$CellCharterR$identity_alignment[[reduced_dim_name]] <- list(
    source = normalizePath(filepath), dataset = name,
    source_cell_order_sha256 = .cc_observation_order_sha256(source),
    target_cell_order_sha256 = .cc_observation_order_sha256(target),
    reordered = alignment$reordered
  )
  spe
}
