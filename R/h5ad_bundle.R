# Fast, resumable, disk-backed H5AD -> SpatialExperiment conversion for R 4.4.
# The source H5AD is cloned byte-for-byte; no Python or anndataR is required.

.transition_require <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
    if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "),
                              ". Install them with BiocManager::install().", call. = FALSE)
}

.transition_h5ls_cache <- new.env(parent = emptyenv())
.transition_children <- function(file, group) {
    key <- normalizePath(file, mustWork = TRUE)
    if (!exists(key, envir = .transition_h5ls_cache, inherits = FALSE)) {
        assign(key, rhdf5::h5ls(file, recursive = TRUE, datasetinfo = FALSE),
               envir = .transition_h5ls_cache)
    }
    listing <- get(key, envir = .transition_h5ls_cache, inherits = FALSE)
    listing[listing$group == group, , drop = FALSE]
}

.transition_file_samples <- function(path, width = 4096L) {
    size <- unname(file.info(path)$size)
    offsets <- unique(pmax(0, c(0, floor(size / 2), size - width)))
    connection <- file(path, "rb")
    on.exit(close(connection), add = TRUE)
    lapply(offsets, function(offset) {
        seek(connection, where = offset, origin = "start")
        readBin(connection, "raw", n = min(width, size - offset))
    })
}

.transition_source_identity <- function(path) {
    info <- file.info(path)
    list(basename = basename(path), size = unname(info$size),
         mtime = as.numeric(info$mtime), samples = .transition_file_samples(path))
}

.transition_identity_equal <- function(x, y) {
    identical(x$size, y$size) && identical(x$mtime, y$mtime) && identical(x$samples, y$samples)
}

.transition_atomic_save_rds <- function(object, path) {
    temporary <- paste0(path, ".tmp-", Sys.getpid())
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(object, temporary, compress = FALSE)
    if (!file.rename(temporary, path)) stop("Could not publish ", path, call. = FALSE)
    invisible(path)
}

.transition_clone_file <- function(source, storage_dir, chunk_mib = 128L) {
    final <- file.path(storage_dir, "arrays.h5")
    partial <- paste0(final, ".partial")
    state_path <- file.path(storage_dir, "copy-state.rds")
    source_identity <- .transition_source_identity(source)
    if (file.exists(final)) {
        manifest_path <- file.path(storage_dir, "manifest.rds")
        if (!file.exists(manifest_path)) stop("arrays.h5 exists without a manifest; use a new output directory.", call. = FALSE)
        manifest <- readRDS(manifest_path)
        if (!.transition_identity_equal(manifest$source, source_identity)) {
            stop("The bundle belongs to a different or changed source H5AD.", call. = FALSE)
        }
        return(normalizePath(final, mustWork = TRUE))
    }
    if (file.exists(partial)) {
        if (!file.exists(state_path)) stop("Partial clone has no copy-state.rds; quarantine it or use a new directory.", call. = FALSE)
        state <- readRDS(state_path)
        if (!.transition_identity_equal(state$source, source_identity)) stop("Partial clone belongs to a changed source.", call. = FALSE)
    } else {
        state <- list(source = source_identity, started = Sys.time(), bytes = 0)
        .transition_atomic_save_rds(state, state_path)
    }
    offset <- if (file.exists(partial)) unname(file.info(partial)$size) else 0
    if (offset > source_identity$size) stop("Partial clone is larger than its source.", call. = FALSE)
    message(sprintf("Cloning %.2f GiB; resuming at %.2f GiB...", source_identity$size / 1024^3, offset / 1024^3))
    input <- file(source, "rb")
    output <- file(partial, if (offset) "ab" else "wb")
    on.exit(try(close(input), silent = TRUE), add = TRUE)
    on.exit(try(close(output), silent = TRUE), add = TRUE)
    seek(input, where = offset, origin = "start")
    chunk <- as.integer(chunk_mib * 1024^2)
    started <- proc.time()[["elapsed"]]
    repeat {
        buffer <- readBin(input, "raw", n = chunk)
        if (!length(buffer)) break
        writeBin(buffer, output)
        offset <- offset + length(buffer)
        flush(output)
        state$bytes <- offset
        state$updated <- Sys.time()
        .transition_atomic_save_rds(state, state_path)
        message(sprintf("  %.2f / %.2f GiB (%.1f%%)", offset / 1024^3,
                        source_identity$size / 1024^3, 100 * offset / source_identity$size))
    }
    close(input); close(output); on.exit(NULL, add = FALSE)
    if (unname(file.info(partial)$size) != source_identity$size ||
        !identical(.transition_file_samples(partial), source_identity$samples)) {
        stop("Completed clone failed byte-sample validation.", call. = FALSE)
    }
    if (!.transition_identity_equal(.transition_source_identity(source), source_identity)) {
        stop("Source changed during conversion; refusing to publish.", call. = FALSE)
    }
    probe <- rhdf5::h5ls(partial, recursive = FALSE, datasetinfo = FALSE)
    if (!all(c("X", "obs", "var", "obsm") %in% probe$name)) stop("Clone is not a readable AnnData HDF5 file.", call. = FALSE)
    if (!file.rename(partial, final)) stop("Could not atomically publish arrays.h5.", call. = FALSE)
    manifest <- list(format_version = 1L, source = source_identity, completed = Sys.time(),
                     copy_seconds = proc.time()[["elapsed"]] - started, arrays_file = "arrays.h5")
    .transition_atomic_save_rds(manifest, file.path(storage_dir, "manifest.rds"))
    unlink(state_path)
    normalizePath(final, mustWork = TRUE)
}

.transition_read_categorical <- function(file, group, expected_length) {
    categories <- rhdf5::h5read(file, paste0(group, "/categories"))
    codes <- as.integer(rhdf5::h5read(file, paste0(group, "/codes")))
    if (length(codes) != expected_length) stop("Unexpected categorical length: ", group, call. = FALSE)
    answer <- rep.int(NA_integer_, length(codes)); present <- codes >= 0L
    if (any(codes[present] >= length(categories))) stop("Invalid categorical code: ", group, call. = FALSE)
    answer[present] <- codes[present] + 1L
    attrs <- tryCatch(rhdf5::h5readAttributes(file, group), error = function(e) list())
    factor(answer, levels = seq_along(categories), labels = as.character(categories), ordered = isTRUE(attrs$ordered))
}

.transition_read_frame <- function(file, group, n, index_path, columns = NULL) {
    entries <- .transition_children(file, group); available <- setdiff(entries$name, "_index")
    frame_attrs <- tryCatch(rhdf5::h5readAttributes(file, group),
                            error = function(e) list())
    declared_order <- frame_attrs[["column-order"]]
    if (!is.null(declared_order)) {
        declared_order <- as.character(declared_order)
        if (!setequal(declared_order, available)) {
            stop("AnnData dataframe `", group,
                 "` has an invalid column-order attribute.", call. = FALSE)
        }
        available <- declared_order
    }
    if (!is.null(columns)) {
        unknown <- setdiff(columns, available)
        if (length(unknown)) stop("Unknown ", group, " columns: ", paste(unknown, collapse = ", "), call. = FALSE)
        available <- columns
    }
    values <- lapply(available, function(name) {
        path <- paste0(group, "/", name); entry <- entries[match(name, entries$name), , drop = FALSE]
        if (identical(entry$otype, "H5I_GROUP")) {
            children <- .transition_children(file, path)$name
            if (all(c("categories", "codes") %in% children)) return(.transition_read_categorical(file, path, n))
            stop("Unsupported AnnData dataframe encoding at `", path,
                 "`; no column was silently omitted.", call. = FALSE)
        }
        tryCatch(HDF5Array::HDF5Array(file, sub("^/", "", path)), error = function(e) {
            if (n > 1e6) stop("Cannot lazily represent ", path, ": ", conditionMessage(e), call. = FALSE)
            rhdf5::h5read(file, sub("^/", "", path))
        })
    })
    names(values) <- available
    ids <- as.character(rhdf5::h5read(file, index_path))
    if (length(ids) != n || anyNA(ids) || anyDuplicated(ids)) stop(index_path, " must contain unique identifiers.", call. = FALSE)
    answer <- do.call(S4Vectors::DataFrame, c(values, list(check.names = FALSE))); rownames(answer) <- ids
    answer
}

.transition_dataframe_obsm <- function(file, group, n) {
    entries <- .transition_children(file, group)
    names <- setdiff(entries$name[entries$otype == "H5I_DATASET"], "_index")
    columns <- lapply(names, function(name) {
        value <- HDF5Array::HDF5Array(file, paste0(sub("^/", "", group), "/", name))
        dim(value) <- c(n, 1L); value
    })
    if (!length(columns)) return(NULL)
    answer <- do.call(S4Arrays::acbind, columns); colnames(answer) <- names; answer
}

.transition_read_obsp <- function(file, group, n_spots) {
    entries <- .transition_children(file, group)
    required <- c("data", "indices", "indptr")
    if (!all(required %in% entries$name)) {
        stop("AnnData sparse matrix `", group,
             "` lacks data, indices, or indptr.", call. = FALSE)
    }
    attrs <- tryCatch(rhdf5::h5readAttributes(file, group),
                      error = function(e) list())
    encoding <- attrs[["encoding-type"]]
    if (is.null(encoding)) encoding <- attrs[["encoding_type"]]
    if (is.null(encoding)) {
        stop("AnnData obsp matrix `", group,
             "` lacks an encoding-type attribute; refusing to guess CSR versus CSC.",
             call. = FALSE)
    }
    encoding <- as.character(encoding)[1L]
    if (!encoding %in% c("csr_matrix", "csc_matrix")) {
        stop("Unsupported AnnData obsp encoding `", encoding, "` at ",
             group, ".", call. = FALSE)
    }
    shape <- attrs$shape
    if (is.null(shape) && "shape" %in% entries$name) {
        shape <- rhdf5::h5read(file, paste0(group, "/shape"))
    }
    if (is.null(shape)) {
        stop("AnnData obsp matrix `", group,
             "` lacks its required shape.", call. = FALSE)
    }
    if (anyNA(as.integer(shape)) ||
        !identical(as.integer(shape), c(n_spots, n_spots))) {
        stop("AnnData obsp matrix `", group, "` is not N-by-N.", call. = FALSE)
    }
    weights <- as.numeric(rhdf5::h5read(file, paste0(group, "/data")))
    indices <- as.integer(rhdf5::h5read(file, paste0(group, "/indices")))
    indptr <- as.integer(rhdf5::h5read(file, paste0(group, "/indptr")))
    if (anyNA(indices) || anyNA(indptr) ||
        length(weights) != length(indices) || length(indptr) != n_spots + 1L ||
        indptr[1L] != 0L || tail(indptr, 1L) != length(indices) ||
        any(diff(indptr) < 0L) || any(indices < 0L | indices >= n_spots) ||
        any(!is.finite(weights))) {
        stop("Invalid AnnData sparse matrix arrays at `", group, "`.", call. = FALSE)
    }
    major <- rep.int(seq_len(n_spots), diff(indptr))
    minor <- indices + 1L
    if (identical(encoding, "csr_matrix")) {
        from <- major; to <- minor
    } else {
        from <- minor; to <- major
    }
    answer <- S4Vectors::SelfHits(from = from, to = to, nnode = n_spots,
                                  value = weights, sort.by.query = FALSE)
    S4Vectors::metadata(answer)$anndata_obsp <- list(
        encoding_type = encoding,
        encoding_version = as.character(attrs[["encoding-version"]])[1L],
        shape = as.integer(shape), nnz = length(weights),
        value_min = if (length(weights)) min(weights) else NA_real_,
        value_max = if (length(weights)) max(weights) else NA_real_,
        explicit_zeros = sum(weights == 0),
        diagonal_entries = sum(from == to)
    )
    answer
}

.transition_build_spe <- function(arrays_file, assay_path, assay_name, obs_columns, var_columns,
                                  spatial_path, reduced_dims, aucell_path, sample_column) {
    assay <- HDF5Array::HDF5Array(arrays_file, assay_path)
    n_features <- nrow(assay); n_spots <- ncol(assay)
    col_data <- .transition_read_frame(arrays_file, "/obs", n_spots, "obs/_index", obs_columns)
    row_data <- .transition_read_frame(arrays_file, "/var", n_features, "var/_index", var_columns)
    if (!is.null(sample_column) && !sample_column %in% colnames(col_data)) {
        stop("sample_column must be imported and present; use NULL only for known single-sample data.", call. = FALSE)
    }
    rownames(assay) <- rownames(row_data); colnames(assay) <- rownames(col_data)
    coords <- HDF5Array::HDF5Array(arrays_file, spatial_path)
    if (!identical(dim(coords), c(2L, n_spots))) stop("Expected a 2-by-N spatial matrix.", call. = FALSE)
    spatial <- t(as.matrix(coords)); colnames(spatial) <- c("x", "y"); rownames(spatial) <- rownames(col_data)
    sample_id <- if (is.null(sample_column)) rep.int("sample01", n_spots) else as.character(col_data[[sample_column]])
    spe <- SpatialExperiment::SpatialExperiment(
        assays = S4Vectors::SimpleList(stats::setNames(list(assay), assay_name)), rowData = row_data,
        colData = col_data, spatialCoords = spatial, sample_id = sample_id)
    for (name in names(reduced_dims)) {
        value <- DelayedArray::t(HDF5Array::HDF5Array(arrays_file, unname(reduced_dims[[name]])))
        rownames(value) <- colnames(spe); SingleCellExperiment::reducedDim(spe, name) <- value
    }
    if (!is.null(aucell_path)) {
        aucell <- .transition_dataframe_obsm(arrays_file, paste0("/", sub("^/", "", aucell_path)), n_spots)
        rownames(aucell) <- colnames(spe); SingleCellExperiment::reducedDim(spe, "AUCell") <- aucell
    }
    graphs <- .transition_children(arrays_file, "/obsp")
    unsupported_graphs <- graphs$name[graphs$otype != "H5I_GROUP"]
    if (length(unsupported_graphs)) {
        stop("Unsupported dense or non-sparse AnnData obsp entries: ",
             paste(unsupported_graphs, collapse = ", "),
             ". No obsp entry was silently omitted.", call. = FALSE)
    }
    graph_names <- graphs$name
    graph_report <- stats::setNames(vector("list", length(graph_names)), graph_names)
    sample_codes <- match(sample_id, unique(sample_id))
    for (name in graph_names) {
        pair <- .transition_read_obsp(
            arrays_file, paste0("/obsp/", name), n_spots
        )
        from <- S4Vectors::from(pair); to <- S4Vectors::to(pair)
        cross_sample <- FALSE
        if (length(pair)) {
            for (start in seq.int(1L, length(pair), by = 1000000L)) {
                index <- start:min(length(pair), start + 999999L)
                if (any(sample_codes[from[index]] != sample_codes[to[index]])) {
                    cross_sample <- TRUE
                    break
                }
            }
        }
        if (cross_sample) {
            stop("AnnData obsp matrix `", name,
                 "` contains edges between different samples in `",
                 sample_column, "`.", call. = FALSE)
        }
        graph_report[[name]] <- c(
            S4Vectors::metadata(pair)$anndata_obsp,
            list(cross_sample_edges = 0L)
        )
        SingleCellExperiment::colPair(spe, name) <- pair
        rm(pair, from, to)
    }
    uns <- .transition_children(arrays_file, "/uns")
    # Keep potentially huge images and other unstructured payloads on disk.
    # The complete /uns hierarchy remains in arrays.h5 and this locator is
    # rebuilt whenever the portable directory is moved.
    S4Vectors::metadata(spe)$anndata_uns <- list(
        filepath = arrays_file, group = "/uns", entries = uns$name
    )
    S4Vectors::metadata(spe)$transitionR <- list(
        arrays_file = arrays_file, independent_of_source = TRUE, complete_h5ad_clone = TRUE,
        anndata_mapping = list(
            obs = colnames(col_data), var = colnames(row_data),
            obsm = c(spatial = spatial_path, reduced_dims, AUCell = aucell_path),
            obsp = graph_names,
            uns = uns$name
        ),
        obsp_groups = graph_names, obsp_report = graph_report,
        uns_entries = uns$name)
    valid <- methods::validObject(spe, test = TRUE)
    if (!identical(valid, TRUE)) stop("Invalid SPE: ", paste(valid, collapse = "; "), call. = FALSE)
    spe
}

.transition_rebase <- function(x, arrays_file) {
    if (!inherits(x, "DelayedArray")) return(x)
    DelayedArray::modify_seeds(x, function(seed) {
        if (inherits(seed, "HDF5ArraySeed")) methods::slot(seed, "filepath") <- arrays_file
        seed
    })
}

# Convert a large H5AD to an independent disk-backed SpatialExperiment bundle.
read_transition_h5ad_spe <- function(
        h5ad_file, output_dir, assay_path = "X", assay_name = "expression",
        obs_columns = NULL, var_columns = NULL, spatial_path = "obsm/spatial",
        reduced_dims = c(PCA = "obsm/X_pca", cellcharter = "obsm/X_cellcharter"),
        aucell_path = "obsm/score_aucell", sample_column = "slide_id_concat",
        chunk_mib = 128L, overwrite = FALSE) {
    .transition_require(c("DelayedArray", "HDF5Array", "rhdf5", "S4Arrays", "S4Vectors",
                          "SingleCellExperiment", "SpatialExperiment", "SummarizedExperiment"))
    h5ad_file <- normalizePath(h5ad_file, mustWork = TRUE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE); output_dir <- normalizePath(output_dir)
    lock <- file.path(output_dir, ".conversion-lock")
    if (!dir.create(lock, showWarnings = FALSE)) stop("Another conversion may be using this bundle.", call. = FALSE)
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    if (file.exists(file.path(output_dir, "spe_backing.h5"))) {
        stop("Legacy/incomplete spe_backing.h5 found; quarantine it before conversion.", call. = FALSE)
    }
    arrays_file <- .transition_clone_file(h5ad_file, output_dir, chunk_mib)
    spe <- .transition_build_spe(arrays_file, assay_path, assay_name, obs_columns, var_columns,
                                 spatial_path, reduced_dims, aucell_path, sample_column)
    # Keep the reconstruction recipe base-R-only. Serialized Bioconductor S4
    # objects are not reliably portable across releases (for example, newer
    # rowRanges can embed the separately packaged Seqinfo class).
    manifest_path <- file.path(output_dir, "manifest.rds")
    manifest <- readRDS(manifest_path)
    manifest$format_version <- 2L
    manifest$import_spec <- list(
        assay_path = assay_path, assay_name = assay_name,
        obs_columns = obs_columns, var_columns = var_columns,
        spatial_path = spatial_path, reduced_dims = reduced_dims,
        aucell_path = aucell_path, sample_column = sample_column
    )
    .transition_atomic_save_rds(manifest, manifest_path)
    S4Vectors::metadata(spe)$transitionR$bundle_dir <- output_dir
    rds_path <- file.path(output_dir, "spe.rds")
    if (file.exists(rds_path) && !overwrite) stop("spe.rds exists; set overwrite = TRUE.", call. = FALSE)
    .transition_atomic_save_rds(spe, rds_path); spe
}

# Load a moved SPE bundle by rebuilding its S4 shell against local packages.
load_transition_spe <- function(storage_dir) {
    storage_dir <- normalizePath(storage_dir, mustWork = TRUE)
    arrays_file <- normalizePath(file.path(storage_dir, "arrays.h5"), mustWork = TRUE)
    manifest_path <- file.path(storage_dir, "manifest.rds")
    if (!file.exists(manifest_path)) {
        stop("Bundle lacks manifest.rds and cannot be reconstructed portably.", call. = FALSE)
    }
    manifest <- readRDS(manifest_path)
    spec <- manifest$import_spec
    if (is.null(spec)) {
        # Format-1 bundles were created before reconstruction parameters were
        # recorded. They used these public defaults; rebuilding avoids loading
        # their release-specific serialized S4 skeleton.
        spec <- list(
            assay_path = "X", assay_name = "expression",
            obs_columns = NULL, var_columns = NULL,
            spatial_path = "obsm/spatial",
            reduced_dims = c(PCA = "obsm/X_pca", cellcharter = "obsm/X_cellcharter"),
            aucell_path = "obsm/score_aucell", sample_column = "slide_id_concat"
        )
    }
    spe <- do.call(
        .transition_build_spe,
        c(list(arrays_file = arrays_file), spec)
    )
    S4Vectors::metadata(spe)$transitionR$arrays_file <- arrays_file
    S4Vectors::metadata(spe)$transitionR$bundle_dir <- storage_dir
    S4Vectors::metadata(spe)$transitionR$manifest_format_version <-
        manifest$format_version
    S4Vectors::metadata(spe)$transitionR$reconstructed_locally <- TRUE
    valid <- methods::validObject(spe, test = TRUE)
    if (!identical(valid, TRUE)) stop("Reconstructed SPE is invalid: ", paste(valid, collapse = "; "), call. = FALSE)
    spe
}

#' Import an H5AD as an independent disk-backed SpatialExperiment
#'
#' The complete H5AD is cloned byte-for-byte into a resumable directory before a
#' [SpatialExperiment::SpatialExperiment()] is constructed. Large arrays remain
#' on disk. Python, reticulate, anndataR, and zellkonverter are not required.
#'
#' @param h5ad_file Path to the source H5AD.
#' @param output_dir Destination directory for `arrays.h5`, `spe.rds`, and the
#'   conversion manifest. Move these files together.
#' @param ... Additional importer arguments such as `sample_column`,
#'   `reduced_dims`, and `chunk_mib`.
#' @return A disk-backed `SpatialExperiment`.
#' @export
#' @examples
#' \dontrun{
#' spe <- cellcharter_import_h5ad("large.h5ad", "large-spe")
#' }
cellcharter_import_h5ad <- function(h5ad_file, output_dir, ...) {
    read_transition_h5ad_spe(h5ad_file = h5ad_file, output_dir = output_dir, ...)
}

#' Load a movable disk-backed SpatialExperiment
#'
#' Reconstructs the S4 experiment shell with the locally installed
#' Bioconductor classes and attaches it to the saved `arrays.h5`. This
#' avoids cross-release serialization dependencies such as `Seqinfo`. Use this
#' function after moving the complete directory.
#'
#' @param path Directory created by [cellcharter_import_h5ad()].
#' @return A valid disk-backed `SpatialExperiment`.
#' @export
#' @examples
#' \dontrun{
#' spe <- cellcharter_load_spe("large-spe")
#' }
cellcharter_load_spe <- function(path) {
    load_transition_spe(path)
}
