# Portable completed-analysis directorys.  Large rectangular data are rebuilt from
# analysis-owned HDF5 datasets; the RDS files contain structure and small results
# only and therefore never serialize an HDF5 seed pointing outside the bundle.

.cc_analysis_version <- 1L

.cc_analysis_progress <- function(enabled, step, message) {
  if (isTRUE(enabled)) message(sprintf("CellCharterR analysis [%s] %s", step, message))
}

.cc_analysis_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  if (!nzchar(x)) "unnamed" else x
}

.cc_analysis_sha256 <- function(path) {
  unname(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

.cc_analysis_file_record <- function(path, root) {
  rel <- substring(normalizePath(path, mustWork = TRUE),
                   nchar(normalizePath(root, mustWork = TRUE)) + 2L)
  list(path = .cc_analysis_canonical_path(rel),
       size = unname(file.info(path)$size), sha256 = .cc_analysis_sha256(path))
}

.cc_analysis_resolve <- function(root, relative, mustWork = TRUE) {
  if (!is.character(relative) || length(relative) != 1L || is.na(relative) ||
      !nzchar(relative) || grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", relative) ||
      grepl("^[/\\\\]", relative) || grepl("^[A-Za-z]:[/\\\\]", relative)) {
    stop("Unsafe relative path in analysis directory: ", relative, call. = FALSE)
  }
  relative <- .cc_analysis_canonical_path(relative)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  candidate <- file.path(root, relative)
  resolved <- normalizePath(candidate, winslash = "/", mustWork = mustWork)
  if (!identical(resolved, root) && !startsWith(resolved, paste0(root, "/"))) {
    stop("Analysis-directory path escapes its analysis directory: ", relative, call. = FALSE)
  }
  resolved
}

.cc_analysis_write_matrix <- function(x, directory, id, progress) {
  if (length(dim(x)) != 2L)
    stop("Bundle matrix component `", id, "` is not two-dimensional.", call. = FALSE)
  file <- file.path(directory, paste0(.cc_analysis_safe_name(id), ".h5"))
  dataset <- "data"
  .cc_analysis_progress(progress, "write", sprintf("%s [%s x %s]", id, nrow(x), ncol(x)))
  HDF5Array::writeHDF5Array(x, filepath = file, name = dataset, level = 0L,
                            with.dimnames = FALSE, verbose = isTRUE(progress))
  list(file = basename(file), dataset = dataset, dim = as.integer(dim(x)),
       dimnames = dimnames(x), class = class(x))
}

.cc_analysis_is_absolute <- function(x) {
  grepl("^[/\\\\]", x) | grepl("^[A-Za-z]:[/\\\\]", x) |
    grepl("^~([/\\\\]|$)", x)
}

.cc_analysis_canonical_path <- function(x) gsub("\\\\", "/", x)

.cc_analysis_pathish_key <- function(key) {
  !is.null(key) && grepl(
    "(^|_)(path|file|filepath|directory|dir|checkpoint|artifact|backing)(_|$)",
    key, ignore.case = TRUE
  )
}

.cc_analysis_new_registry <- function() {
  registry <- new.env(parent = emptyenv())
  registry$entries <- list()
  registry
}

.cc_analysis_registry_add <- function(registry, source, relative) {
  registry$entries[[length(registry$entries) + 1L]] <- list(
    source = source, relative = relative, directory = dir.exists(source)
  )
  invisible(NULL)
}

.cc_analysis_registry_find <- function(registry, source) {
  if (is.null(registry) || !length(registry$entries)) return(NULL)
  candidates <- Filter(function(entry) {
    identical(source, entry$source) ||
      (isTRUE(entry$directory) &&
         startsWith(source, paste0(entry$source, "/")))
  }, registry$entries)
  if (!length(candidates)) return(NULL)
  # Prefer the most specific registered ancestor when directories are nested.
  candidates[[which.max(vapply(candidates, function(x) nchar(x$source), integer(1)))]]
}

.cc_analysis_pack_paths <- function(
    x, directory, key = NULL, breadcrumb = "metadata", owned = FALSE,
    registry = .cc_analysis_new_registry()) {
  if (inherits(x, "cellcharter_bundle_path") ||
      inherits(x, "cellcharter_bundle_paths")) return(x)
  if (is.environment(x) || typeof(x) == "externalptr" || inherits(x, "connection")) {
    stop("Unsupported opaque external reference at `", breadcrumb,
         "`; remove it or provide a portable representation before saving.", call. = FALSE)
  }
  if (inherits(x, "cellcharter_autok") && !is.null(x$checkpoint_dir) &&
      file.exists(x$checkpoint_dir)) {
    old_root <- normalizePath(x$checkpoint_dir, mustWork = TRUE)
    packed_root <- .cc_analysis_pack_paths(
      x$checkpoint_dir, directory, "checkpoint_dir",
      paste0(breadcrumb, "$checkpoint_dir"), owned = TRUE,
      registry = registry
    )
    for (field in intersect(c("fit_index", "best_index"), names(x))) {
      if (is.data.frame(x[[field]]) && "path" %in% names(x[[field]])) {
        old <- as.character(x[[field]]$path)
        inside <- startsWith(normalizePath(old, mustWork = FALSE), paste0(old_root, .Platform$file.sep))
        if (!all(inside))
          stop("AutoK `", field, "$path` contains a fit outside `checkpoint_dir`.",
               call. = FALSE)
        relative <- substring(normalizePath(old, mustWork = FALSE), nchar(old_root) + 2L)
        x[[field]]$path <- structure(.cc_analysis_canonical_path(
          file.path(as.character(packed_root), relative)
        ),
                                     class = c("cellcharter_bundle_paths", "character"))
      }
    }
    x$checkpoint_dir <- packed_root
  }
  if (methods::is(x, "DataFrame") || is.data.frame(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      child_key <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else key
      x[[i]] <- .cc_analysis_pack_paths(
        x[[i]], directory, child_key,
        paste0(breadcrumb, "$", if (is.null(child_key)) i else child_key),
        owned = owned, registry = registry
      )
    }
    return(x)
  }
  if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      child_key <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else key
      child_owned <- owned || identical(child_key, "CellCharterR") ||
        identical(child_key, "transitionR")
      x[[i]] <- .cc_analysis_pack_paths(
        x[[i]], directory, child_key,
        paste0(breadcrumb, "$", if (is.null(child_key)) i else child_key),
        owned = child_owned, registry = registry
      )
    }
    return(x)
  }
  if (base::isS4(x)) {
    for (slot_name in methods::slotNames(x)) {
      packed <- .cc_analysis_pack_paths(
        methods::slot(x, slot_name), directory, slot_name,
        paste0(breadcrumb, "@", slot_name), owned = owned,
        registry = registry
      )
      tryCatch(
        methods::slot(x, slot_name) <- packed,
        error = function(e) stop(
          "Unsupported S4 metadata at `", breadcrumb, "@", slot_name,
          "`: ", conditionMessage(e), call. = FALSE
        )
      )
    }
    return(x)
  }
  if (!is.character(x)) return(x)
  pathish <- .cc_analysis_pathish_key(key)
  absolute <- !is.na(x) & .cc_analysis_is_absolute(x)
  if (any(absolute) && !owned) {
    stop("Unsupported absolute metadata value at `", breadcrumb,
         "`; remove it or store portable provenance instead.", call. = FALSE)
  }
  if (any(absolute) && owned && !isTRUE(pathish)) {
    # Preserve informative provenance without retaining a machine-specific
    # source location or copying a large file already exported as a component.
    x[absolute] <- basename(.cc_analysis_canonical_path(path.expand(x[absolute])))
    return(x)
  }
  if (isTRUE(pathish) && !owned && any(!is.na(x) & nzchar(x))) {
    stop("Unsupported external metadata path at `", breadcrumb,
         "`; move the referenced data into a supported component or remove the field.",
         call. = FALSE)
  }
  if (length(x) != 1L || !isTRUE(pathish) || !file.exists(path.expand(x))) {
    if (owned && isTRUE(pathish) && any(!is.na(x) & nzchar(x))) {
      stop("CellCharterR metadata path at `", breadcrumb,
           "` does not exist and cannot be bundled: ", x[which(!is.na(x) & nzchar(x))[1L]],
           call. = FALSE)
    }
    return(x)
  }
  source <- normalizePath(path.expand(x), winslash = "/", mustWork = TRUE)
  registered <- .cc_analysis_registry_find(registry, source)
  if (!is.null(registered)) {
    suffix <- if (identical(source, registered$source)) "" else
      substring(source, nchar(registered$source) + 2L)
    rel <- if (nzchar(suffix)) file.path(registered$relative, suffix) else
      registered$relative
    return(structure(.cc_analysis_canonical_path(rel),
                     class = c("cellcharter_bundle_path", "character")))
  }
  token <- substr(digest::digest(source, algo = "sha256"), 1L, 16L)
  parent <- file.path(directory, "artifacts", token)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  rel <- .cc_analysis_canonical_path(file.path("artifacts", token, basename(source)))
  destination <- file.path(directory, rel)
  if (!file.exists(destination) &&
      !file.copy(source, parent, recursive = TRUE, copy.mode = TRUE,
                 copy.date = TRUE)) {
    stop("Could not copy CellCharterR artifact referenced by `", key, "`: ", source,
         call. = FALSE)
  }
  .cc_analysis_registry_add(registry, source, rel)
  structure(rel, class = c("cellcharter_bundle_path", "character"))
}

.cc_analysis_unpack_paths <- function(x, directory) {
  if (inherits(x, "cellcharter_bundle_paths")) {
    rel <- .cc_analysis_canonical_path(as.character(x))
    if (any(grepl("(^|/)\\.\\.(/|$)", rel)) || any(grepl("^/", rel)))
      stop("Unsafe bundled artifact path in AutoK fit index.", call. = FALSE)
    paths <- vapply(rel, function(one) .cc_analysis_resolve(directory, one, FALSE),
                    character(1))
    if (any(!file.exists(paths)))
      stop("One or more bundled AutoK fit artifacts are missing.", call. = FALSE)
    return(paths)
  }
  if (inherits(x, "cellcharter_bundle_path")) {
    rel <- .cc_analysis_canonical_path(as.character(x))
    if (grepl("(^|/)\\.\\.(/|$)", rel) || grepl("^/", rel))
      stop("Unsafe bundled artifact path: ", rel, call. = FALSE)
    path <- .cc_analysis_resolve(directory, rel, mustWork = FALSE)
    if (!file.exists(path)) stop("Bundled CellCharterR artifact is missing: ", rel,
                                 call. = FALSE)
    return(path)
  }
  if (methods::is(x, "DataFrame") || is.data.frame(x)) {
    for (i in seq_along(x)) x[[i]] <- .cc_analysis_unpack_paths(x[[i]], directory)
    return(x)
  }
  if (is.list(x)) {
    for (i in seq_along(x)) x[[i]] <- .cc_analysis_unpack_paths(x[[i]], directory)
    return(x)
  }
  if (base::isS4(x)) {
    for (slot_name in methods::slotNames(x)) {
      methods::slot(x, slot_name) <- .cc_analysis_unpack_paths(
        methods::slot(x, slot_name), directory
      )
    }
    return(x)
  }
  x
}

.cc_analysis_seqfun <- function(name) {
  # Bioconductor 3.22 split Seqinfo out of GenomeInfoDb. Build the namespace
  # names dynamically so R CMD check does not require both mutually
  # version-specific packages to be declared simultaneously.
  seqinfo_package <- paste0("Seq", "info")
  genome_package <- paste0("Genome", "InfoDb")
  package <- if (requireNamespace(seqinfo_package, quietly = TRUE)) {
    seqinfo_package
  } else if (requireNamespace(genome_package, quietly = TRUE)) {
    genome_package
  } else {
    stop("Install the Bioconductor Seqinfo/GenomeInfoDb dependency required ",
         "to reconstruct rowRanges.", call. = FALSE)
  }
  getExportedValue(package, name)
}

.cc_analysis_validate_embedded <- function(x, breadcrumb, key = NULL) {
  if (inherits(x, "DelayedArray")) {
    stop("Disk-backed value at `", breadcrumb,
         "` must be moved to an assay or reducedDim before bundling.", call. = FALSE)
  }
  if (is.environment(x) || typeof(x) == "externalptr" || inherits(x, "connection")) {
    stop("Unsupported opaque external reference at `", breadcrumb, "`.", call. = FALSE)
  }
  if (is.character(x) && any(!is.na(x) & .cc_analysis_is_absolute(x))) {
    stop("Unsupported external path at `", breadcrumb, "`.", call. = FALSE)
  }
  if (methods::is(x, "DataFrame") || is.data.frame(x) || is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      child <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
      .cc_analysis_validate_embedded(
        x[[i]], paste0(breadcrumb, "$", child), child
      )
    }
  }
  invisible(TRUE)
}

.cc_analysis_pack_ranges <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "GRangesList") && all(S4Vectors::elementNROWS(x) == 0L)) {
    # The default rowRanges of an ordinary SummarizedExperiment is an empty
    # GRangesList. Its per-row annotations are already represented by rowData;
    # allowing the local constructor to recreate it avoids serializing Seqinfo.
    return(NULL)
  }
  if (!inherits(x, "GRanges")) {
    stop("Portable analysis directorys currently support `GRanges` rowRanges, not `",
         paste(class(x), collapse = "/"), "`.", call. = FALSE)
  }
  .cc_analysis_validate_embedded(S4Vectors::mcols(x), "rowRanges@mcols")
  .cc_analysis_validate_embedded(S4Vectors::metadata(x), "rowRanges@metadata")
  info <- .cc_analysis_seqfun("seqinfo")(x)
  list(
    format = "GRanges-v1", names = names(x),
    seqnames = as.character(GenomicRanges::seqnames(x)),
    start = as.numeric(BiocGenerics::start(x)),
    width = as.numeric(BiocGenerics::width(x)),
    strand = as.character(BiocGenerics::strand(x)),
    mcols = S4Vectors::mcols(x), metadata = S4Vectors::metadata(x),
    seqinfo = list(
      seqnames = .cc_analysis_seqfun("seqlevels")(info),
      seqlengths = .cc_analysis_seqfun("seqlengths")(info),
      isCircular = .cc_analysis_seqfun("isCircular")(info),
      genome = .cc_analysis_seqfun("genome")(info)
    )
  )
}

.cc_analysis_unpack_ranges <- function(x) {
  if (is.null(x)) return(NULL)
  if (!identical(x$format, "GRanges-v1")) {
    stop("Unsupported portable rowRanges representation.", call. = FALSE)
  }
  info <- do.call(.cc_analysis_seqfun("Seqinfo"), x$seqinfo)
  answer <- GenomicRanges::GRanges(
    seqnames = x$seqnames,
    ranges = IRanges::IRanges(start = x$start, width = x$width),
    strand = x$strand, seqinfo = info
  )
  names(answer) <- x$names
  S4Vectors::mcols(answer) <- x$mcols
  S4Vectors::metadata(answer) <- x$metadata
  answer
}

.cc_analysis_pack_pair <- function(x) {
  .cc_analysis_validate_embedded(S4Vectors::mcols(x), "colPair@mcols")
  .cc_analysis_validate_embedded(S4Vectors::metadata(x), "colPair@metadata")
  list(
    format = "SelfHits-v1", from = S4Vectors::from(x), to = S4Vectors::to(x),
    nnode = S4Vectors::nLnode(x), mcols = S4Vectors::mcols(x),
    metadata = S4Vectors::metadata(x)
  )
}

.cc_analysis_unpack_pair <- function(x, n_cells, name) {
  if (!identical(x$format, "SelfHits-v1") || !identical(as.integer(x$nnode),
                                                          as.integer(n_cells)) ||
      any(x$from < 1L | x$from > n_cells | x$to < 1L | x$to > n_cells)) {
    stop("Invalid bundled colPair `", name, "`.", call. = FALSE)
  }
  answer <- S4Vectors::SelfHits(
    from = x$from, to = x$to, nnode = n_cells
  )
  S4Vectors::mcols(answer) <- x$mcols
  S4Vectors::metadata(answer) <- x$metadata
  answer
}

.cc_analysis_save_sce <- function(
    x, directory, prefix, progress, artifact_root = directory,
    registry = .cc_analysis_new_registry()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  .cc_analysis_validate_embedded(SummarizedExperiment::rowData(x),
                               paste0(prefix, "$rowData"))
  .cc_analysis_validate_embedded(SummarizedExperiment::colData(x),
                               paste0(prefix, "$colData"))
  assay_names <- SummarizedExperiment::assayNames(x)
  assays <- lapply(seq_along(assay_names), function(i) {
    .cc_analysis_write_matrix(SummarizedExperiment::assay(x, i, withDimnames = FALSE),
                            directory, paste0(prefix, "_assay_", i), progress)
  })
  names(assays) <- assay_names

  rd_names <- SingleCellExperiment::reducedDimNames(x)
  reduced <- lapply(seq_along(rd_names), function(i) {
    .cc_analysis_write_matrix(SingleCellExperiment::reducedDim(x, i, withDimnames = FALSE),
                            directory, paste0(prefix, "_reduced_", i), progress)
  })
  names(reduced) <- rd_names

  # Register this experiment's checkpoint roots before descending into
  # alternative experiments, allowing all children to reuse their artifacts.
  packed_metadata <- .cc_analysis_pack_paths(
    S4Vectors::metadata(x), artifact_root, registry = registry
  )

  alt_names <- SingleCellExperiment::altExpNames(x)
  alternatives <- lapply(seq_along(alt_names), function(i) {
    subdir <- file.path(directory, paste0("alt_", i))
    .cc_analysis_save_sce(SingleCellExperiment::altExp(x, i), subdir,
                        paste0(prefix, "_alt_", i), progress,
                        artifact_root = artifact_root, registry = registry)
  })
  names(alternatives) <- alt_names

  pairs <- lapply(SingleCellExperiment::colPairNames(x), function(nm) {
    .cc_analysis_pack_pair(SingleCellExperiment::colPair(x, nm, asSparse = FALSE))
  })
  names(pairs) <- SingleCellExperiment::colPairNames(x)

  structure <- list(
    class = class(x), dim = as.integer(dim(x)), dimnames = dimnames(x),
    rowData = SummarizedExperiment::rowData(x),
    rowRanges = if (methods::is(x, "RangedSummarizedExperiment"))
      .cc_analysis_pack_ranges(SummarizedExperiment::rowRanges(x)) else NULL,
    colData = SummarizedExperiment::colData(x),
    metadata = packed_metadata,
    spatialCoords = if (methods::is(x, "SpatialExperiment"))
      SpatialExperiment::spatialCoords(x) else NULL,
    imgData = if (methods::is(x, "SpatialExperiment")) {
      image_data <- SpatialExperiment::imgData(x)
      if (nrow(image_data)) {
        stop("Non-empty `imgData` may contain externally backed images and is not yet ",
             "portable; remove it or materialize images before saving.", call. = FALSE)
      }
      image_data
    } else NULL,
    mainExpName = SingleCellExperiment::mainExpName(x),
    assays = assays, reducedDims = reduced, colPairs = pairs, altExps = alternatives
  )
  saveRDS(structure, file.path(directory, "structure.rds"), version = 3,
          compress = FALSE)
  list(directory = basename(directory), dim = as.integer(dim(x)), class = class(x))
}

.cc_analysis_read_matrix <- function(spec, directory, id) {
  path <- .cc_analysis_resolve(directory, spec$file, mustWork = FALSE)
  if (!file.exists(path)) stop("Analysis component `", id, "` is missing: ", path, call. = FALSE)
  ans <- HDF5Array::HDF5Array(path, spec$dataset)
  if (!identical(as.integer(dim(ans)), as.integer(spec$dim)))
    stop("Analysis component `", id, "` has unexpected dimensions.", call. = FALSE)
  dimnames(ans) <- spec$dimnames
  ans
}

.cc_analysis_load_sce <- function(directory, id = "root", artifact_root = directory) {
  structure_file <- file.path(directory, "structure.rds")
  if (!file.exists(structure_file)) stop("Missing analysis structure: ", structure_file, call. = FALSE)
  z <- readRDS(structure_file)
  z$metadata <- .cc_analysis_unpack_paths(z$metadata, artifact_root)
  assays <- lapply(seq_along(z$assays), function(i)
    .cc_analysis_read_matrix(z$assays[[i]], directory, paste0(id, "/assay/", names(z$assays)[i])))
  names(assays) <- names(z$assays)

  # Constructor-based reconstruction avoids carrying serialized HDF5 seed paths.
  row_args <- if (is.null(z$rowRanges)) list(rowData = z$rowData) else
    list(rowRanges = .cc_analysis_unpack_ranges(z$rowRanges))
  if ("SpatialExperiment" %in% z$class) {
    ans <- do.call(SpatialExperiment::SpatialExperiment, c(
      list(assays = assays, colData = z$colData, metadata = z$metadata,
           spatialCoords = z$spatialCoords, imgData = z$imgData), row_args))
  } else {
    ans <- do.call(SingleCellExperiment::SingleCellExperiment, c(
      list(assays = assays, colData = z$colData, metadata = z$metadata), row_args))
  }
  dimnames(ans) <- z$dimnames

  for (i in seq_along(z$reducedDims)) {
    value <- .cc_analysis_read_matrix(z$reducedDims[[i]], directory,
                                    paste0(id, "/reducedDim/", names(z$reducedDims)[i]))
    if (nrow(value) != ncol(ans)) stop("Reduced dimension row count does not match cells.", call. = FALSE)
    SingleCellExperiment::reducedDim(ans, names(z$reducedDims)[i]) <- value
  }
  for (nm in names(z$colPairs)) {
    SingleCellExperiment::colPair(ans, nm) <-
      .cc_analysis_unpack_pair(z$colPairs[[nm]], ncol(ans), nm)
  }
  for (i in seq_along(z$altExps)) {
    subdir <- .cc_analysis_resolve(directory, z$altExps[[i]]$directory)
    value <- .cc_analysis_load_sce(
      subdir, paste0(id, "/altExp/", names(z$altExps)[i]),
      artifact_root = artifact_root
    )
    if (ncol(value) != ncol(ans)) stop("Alternative experiment cell count does not match.", call. = FALSE)
    SingleCellExperiment::altExp(ans, names(z$altExps)[i]) <- value
  }
  SingleCellExperiment::mainExpName(ans) <- z$mainExpName
  methods::validObject(ans)
  ans
}

.cc_analysis_assert_no_links <- function(root) {
  paths <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                      no.. = TRUE, include.dirs = TRUE)
  links <- nzchar(Sys.readlink(paths))
  if (any(links)) {
    stop("Analysis directorys cannot contain symbolic-link payloads: ",
         paths[which(links)[1L]], call. = FALSE)
  }
  invisible(TRUE)
}

.cc_analysis_inventory <- function(root) {
  .cc_analysis_assert_no_links(root)
  paths <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                      no.. = TRUE)
  paths <- paths[!dir.exists(paths)]
  records <- lapply(paths, .cc_analysis_file_record, root = root)
  names(records) <- vapply(records, `[[`, character(1), "path")
  unname(records)
}

.cc_analysis_validate_files <- function(manifest, root, verify_checksums = TRUE) {
  .cc_analysis_assert_no_links(root)
  expected <- sort(vapply(manifest$files, `[[`, character(1), "path"))
  if (anyDuplicated(expected)) stop("Analysis manifest contains duplicate paths.", call. = FALSE)
  actual <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                       include.dirs = FALSE)
  actual <- sort(setdiff(actual, c("manifest.rds", "manifest.sha256")))
  if (!identical(actual, expected)) {
    stop("Analysis payload inventory differs from its manifest.", call. = FALSE)
  }
  for (entry in manifest$files) {
    path <- .cc_analysis_resolve(root, entry$path, mustWork = FALSE)
    if (!file.exists(path)) stop("Analysis file is missing: ", entry$path, call. = FALSE)
    if (!identical(as.numeric(file.info(path)$size), as.numeric(entry$size)))
      stop("Analysis file size mismatch: ", entry$path, call. = FALSE)
    if (isTRUE(verify_checksums) &&
        !identical(.cc_analysis_sha256(path), entry$sha256))
      stop("Analysis checksum mismatch: ", entry$path, call. = FALSE)
  }
  invisible(TRUE)
}

#' Save a portable completed CellCharterR analysis
#'
#' Large assays and embeddings are copied blockwise into analysis-owned HDF5
#' files. The returned directory can be moved independently of its source H5AD
#' and original backing files.
#'
#' @param spe A `SingleCellExperiment` or `SpatialExperiment`.
#' @param output_dir Destination directory, which must not exist unless
#'   `overwrite = TRUE`.
#' @param overwrite Replace an existing valid destination transactionally.
#' @param progress Show component- and HDF5-block-level progress messages.
#' @return The normalized analysis directory, invisibly.
#' @export
cellcharter_save_analysis <- function(
    spe, output_dir, overwrite = FALSE, progress = interactive()) {
  if (!methods::is(spe, "SingleCellExperiment"))
    stop("`spe` must be a SingleCellExperiment or SpatialExperiment.", call. = FALSE)
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))
    stop("`output_dir` must be one non-empty path.", call. = FALSE)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite) ||
      !is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`overwrite` and `progress` must each be `TRUE` or `FALSE`.", call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  parent <- dirname(output_dir)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
  if (file.exists(output_dir) && !isTRUE(overwrite))
    stop("Analysis destination already exists; use `overwrite = TRUE`: ", output_dir, call. = FALSE)
  stage <- tempfile(paste0(".", basename(output_dir), "-stage-"), tmpdir = parent)
  backup <- paste0(output_dir, ".backup")
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  dir.create(stage, recursive = TRUE)
  .cc_analysis_progress(progress, "1/4", "writing portable components")
  root_spec <- .cc_analysis_save_sce(spe, file.path(stage, "object"), "root", progress)
  manifest <- list(
    format = "CellCharterR-analysis-directory", format_version = .cc_analysis_version,
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    producer = list(
      CellCharterR = tryCatch(as.character(utils::packageVersion("CellCharterR")),
                             error = function(e) NA_character_),
      R = as.character(getRversion()),
      SpatialExperiment = as.character(utils::packageVersion("SpatialExperiment")),
      HDF5Array = as.character(utils::packageVersion("HDF5Array"))),
    root = root_spec, files = list())
  .cc_analysis_progress(progress, "2/4", "computing SHA-256 checksums")
  manifest$files <- .cc_analysis_inventory(stage)
  saveRDS(manifest, file.path(stage, "manifest.rds"), version = 3)
  writeLines(
    .cc_analysis_sha256(file.path(stage, "manifest.rds")),
    file.path(stage, "manifest.sha256"), useBytes = TRUE
  )
  .cc_analysis_progress(progress, "3/4", "validating staged analysis")
  # The inventory pass has just hashed every staged payload. Re-reading all
  # large HDF5 files here would double save-time I/O without adding evidence;
  # full checksum verification remains the default when the bundle is loaded.
  .cc_analysis_validate_files(manifest, stage, verify_checksums = FALSE)
  test <- .cc_analysis_load_sce(file.path(stage, "object"))
  if (!identical(dim(test), dim(spe))) stop("Staged analysis failed dimension validation.", call. = FALSE)
  rm(test)
  .cc_analysis_progress(progress, "4/4", "publishing atomically")
  if (file.exists(output_dir)) {
    if (file.exists(backup)) stop("Cannot overwrite while backup path exists: ", backup, call. = FALSE)
    if (!file.rename(output_dir, backup)) stop("Could not move existing analysis to backup.", call. = FALSE)
    restored <- FALSE
    on.exit(if (!restored && dir.exists(backup) && !dir.exists(output_dir))
      file.rename(backup, output_dir), add = TRUE)
  }
  if (!file.rename(stage, output_dir)) stop("Could not publish staged analysis on this filesystem.", call. = FALSE)
  if (exists("restored", inherits = FALSE)) {
    restored <- TRUE
    unlink(backup, recursive = TRUE)
  }
  invisible(output_dir)
}

#' Load a portable completed CellCharterR analysis
#'
#' @param path Directory created by [cellcharter_save_analysis()].
#' @param validate Verify all payload sizes and SHA-256 checksums before loading.
#' @return A locally reconstructed `SpatialExperiment` or
#'   `SingleCellExperiment` with analysis-owned HDF5 backing.
#' @export
cellcharter_load_analysis <- function(path, validate = TRUE) {
  if (!is.character(path) || length(path) != 1L || !dir.exists(path))
    stop("`path` must name an existing analysis-directory directory.", call. = FALSE)
  if (!is.logical(validate) || length(validate) != 1L || is.na(validate)) {
    stop("`validate` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  root <- normalizePath(path, mustWork = TRUE)
  manifest_file <- file.path(root, "manifest.rds")
  if (!file.exists(manifest_file)) stop("Analysis manifest is missing: ", manifest_file, call. = FALSE)
  manifest_digest <- file.path(root, "manifest.sha256")
  if (!file.exists(manifest_digest)) {
    stop("Analysis manifest checksum is missing: ", manifest_digest, call. = FALSE)
  }
  expected_manifest_digest <- trimws(readLines(manifest_digest, warn = FALSE, n = 1L))
  if (!identical(.cc_analysis_sha256(manifest_file), expected_manifest_digest)) {
    stop("Analysis manifest checksum mismatch; the bundle is modified or corrupt.", call. = FALSE)
  }
  manifest <- readRDS(manifest_file)
  if (!identical(manifest$format, "CellCharterR-analysis-directory"))
    stop("Directory is not a CellCharterR analysis directory.", call. = FALSE)
  if (!identical(as.integer(manifest$format_version), .cc_analysis_version))
    stop("Unsupported analysis-directory format version: ", manifest$format_version, call. = FALSE)
  if (isTRUE(validate)) .cc_analysis_validate_files(manifest, root)
  ans <- .cc_analysis_load_sce(file.path(root, "object"))
  if (!identical(as.integer(dim(ans)), as.integer(manifest$root$dim)))
    stop("Loaded object dimensions disagree with the manifest.", call. = FALSE)
  ans
}
