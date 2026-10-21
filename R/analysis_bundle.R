# Portable completed-analysis directories. Assays are rebuilt from
# analysis-owned HDF5 datasets; native annotations, coordinates, embeddings,
# graphs, and results are stored as RDS payloads. No serialized HDF5 seed points
# outside the directory.

.cc_analysis_version <- 2L

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

.cc_analysis_realize_matrix <- function(x, id, progress, block_mib = 64) {
  if (!inherits(x, "DelayedArray")) return(x)
  if (length(dim(x)) != 2L) {
    stop("Delayed component `", id, "` is not two-dimensional.", call. = FALSE)
  }
  nr <- nrow(x); nc <- ncol(x)
  .cc_analysis_progress(
    progress, "materialize",
    sprintf("%s [%s x %s] as an ordinary matrix", id, nr, nc)
  )
  if (!nr || !nc) return(as.matrix(x))
  # Probe one row to preserve integer/logical storage where possible, then fill
  # the result in bounded row blocks. This avoids `as.matrix()` creating a
  # second cohort-wide temporary object for legacy delayed metadata/embeddings.
  probe <- as.matrix(x[1L, , drop = FALSE])
  answer <- matrix(vector(typeof(probe), as.double(nr) * as.double(nc)),
                   nrow = nr, ncol = nc)
  bytes_per_value <- max(1, utils::object.size(probe) / max(1, length(probe)))
  rows_per_block <- max(
    1L,
    min(nr, floor(block_mib * 1024^2 / max(1, nc * bytes_per_value)))
  )
  for (start in seq.int(1L, nr, by = rows_per_block)) {
    end <- min(nr, start + rows_per_block - 1L)
    answer[start:end, ] <- as.matrix(x[start:end, , drop = FALSE])
  }
  dimnames(answer) <- dimnames(x)
  answer
}

.cc_analysis_write_native_matrix <- function(x, directory, id, progress) {
  if (inherits(x, "DelayedArray")) {
    # Preserve bounded-memory saving for legacy/delayed embeddings. The loader
    # realizes this analysis-owned copy into the canonical native reducedDim;
    # it never follows the original seed path.
    spec <- .cc_analysis_write_matrix(x, directory, id, progress)
    spec$format <- "native-on-load-hdf5-v1"
    return(spec)
  }
  value <- .cc_analysis_realize_matrix(x, id, progress)
  if (!is.matrix(value) || !is.numeric(value)) {
    stop("Reduced dimension `", id,
         "` must materialize to an ordinary numeric matrix.", call. = FALSE)
  }
  file <- file.path(directory, paste0(.cc_analysis_safe_name(id), ".rds"))
  .cc_analysis_progress(
    progress, "write", sprintf("%s [%s x %s] native matrix", id,
                               nrow(value), ncol(value))
  )
  saveRDS(value, file, version = 3, compress = FALSE)
  list(format = "native-matrix-rds-v1", file = basename(file),
       dim = as.integer(dim(value)), dimnames = dimnames(value))
}

.cc_analysis_materialize_value <- function(x, breadcrumb, progress) {
  if (inherits(x, "DelayedArray")) {
    if (length(dim(x)) == 2L) {
      value <- .cc_analysis_realize_matrix(x, breadcrumb, progress)
      if (ncol(value) == 1L) return(drop(value))
      return(value)
    }
    .cc_analysis_progress(progress, "materialize", breadcrumb)
    return(as.vector(as.array(x)))
  }
  if (methods::is(x, "DataFrame") || is.data.frame(x)) {
    for (i in seq_along(x)) {
      child <- if (!is.null(names(x)) && nzchar(names(x)[i])) names(x)[i] else i
      x[[i]] <- .cc_analysis_materialize_value(
        x[[i]], paste0(breadcrumb, "$", child), progress
      )
    }
    return(x)
  }
  if (is.list(x) && !base::isS4(x)) {
    for (i in seq_along(x)) {
      child <- if (!is.null(names(x)) && nzchar(names(x)[i])) names(x)[i] else i
      x[[i]] <- .cc_analysis_materialize_value(
        x[[i]], paste0(breadcrumb, "$", child), progress
      )
    }
  }
  x
}

.cc_analysis_copy_uns <- function(locator, directory, progress) {
  if (!.cc_analysis_is_anndata_uns_locator(locator)) return(locator)
  source <- normalizePath(locator$filepath, mustWork = TRUE)
  token <- substr(digest::digest(
    paste(source, locator$group, sep = "\n"), algo = "sha256"
  ), 1L, 16L)
  relative <- .cc_analysis_canonical_path(
    file.path("artifacts", paste0("anndata-uns-", token, ".h5"))
  )
  destination <- file.path(directory, relative)
  if (!file.exists(destination)) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    .cc_analysis_progress(
      progress, "copy",
      paste0("AnnData ", locator$group, " hierarchy (without the full H5AD clone)")
    )
    rhdf5::h5createFile(destination)
    source_id <- rhdf5::H5Fopen(source)
    destination_id <- rhdf5::H5Fopen(destination, flags = "H5F_ACC_RDWR")
    copy_error <- NULL
    tryCatch(
      rhdf5::H5Ocopy(source_id, locator$group, destination_id, locator$group),
      error = function(e) copy_error <<- e,
      finally = {
        rhdf5::H5Fclose(destination_id)
        rhdf5::H5Fclose(source_id)
      }
    )
    if (!is.null(copy_error)) {
      unlink(destination)
      stop("Could not copy AnnData `", locator$group,
           "` into the portable analysis: ", conditionMessage(copy_error),
           call. = FALSE)
    }
  }
  locator$filepath <- structure(
    relative, class = c("cellcharter_bundle_path", "character")
  )
  locator
}

.cc_analysis_prepare_import_metadata <- function(x, directory, progress) {
  if (!is.list(x)) return(x)
  if (.cc_analysis_is_anndata_uns_locator(x$anndata_uns)) {
    x$anndata_uns <- .cc_analysis_copy_uns(x$anndata_uns, directory, progress)
  }
  if (is.list(x$transitionR)) {
    transition <- x$transitionR
    # These paths point at the import directory, not at scientific results.
    # Copying either would duplicate the complete H5AD (potentially hundreds
    # of GiB) after the assay and /uns have already been exported.
    if (is.character(transition$arrays_file) && length(transition$arrays_file) == 1L) {
      transition$source_arrays_basename <- basename(transition$arrays_file)
    }
    if (is.character(transition$bundle_dir) && length(transition$bundle_dir) == 1L) {
      transition$source_bundle_basename <- basename(transition$bundle_dir)
    }
    transition$arrays_file <- NULL
    transition$bundle_dir <- NULL
    transition$complete_h5ad_clone <- FALSE
    transition$portable_analysis <- TRUE
    if (.cc_analysis_is_anndata_uns_locator(x$anndata_uns)) {
      transition$uns_file <- x$anndata_uns$filepath
    }
    x$transitionR <- transition
  }
  x
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

.cc_analysis_is_anndata_uns_locator <- function(x) {
  is.list(x) && !is.data.frame(x) &&
    all(c("filepath", "group", "entries") %in% names(x)) &&
    is.character(x$filepath) && length(x$filepath) == 1L &&
    !is.na(x$filepath) && nzchar(x$filepath) &&
    is.character(x$group) && length(x$group) == 1L &&
    !is.na(x$group) && startsWith(x$group, "/") &&
    (is.null(x$entries) || is.character(x$entries))
}

.cc_analysis_is_anndata_mapping <- function(x) {
  is.list(x) && !is.data.frame(x) &&
    all(c("obs", "var", "obsm", "obsp", "uns") %in% names(x))
}

.cc_analysis_child_path_role <- function(parent, key) {
  # These roles describe package-owned schemas, not arbitrary similarly named
  # third-party metadata. In particular, `checkpoint_mode` is an AutoK option,
  # while the H5AD importer deliberately stores /uns in its cloned arrays file.
  if (inherits(parent, "cellcharter_autok") &&
      identical(key, "checkpoint_mode")) {
    return("literal")
  }
  if (.cc_analysis_is_anndata_uns_locator(parent) &&
      identical(key, "filepath")) {
    return("owned")
  }
  if (.cc_analysis_is_anndata_uns_locator(parent) &&
      identical(key, "group")) {
    return("literal")
  }
  if (.cc_analysis_is_anndata_mapping(parent) &&
      identical(key, "obsm")) {
    # Values such as /obsm/spatial are HDF5 dataset identifiers inside the
    # relocated arrays file, not host filesystem paths.
    return("literal")
  }
  "auto"
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
    registry = .cc_analysis_new_registry(), path_role = "auto") {
  if (inherits(x, "cellcharter_bundle_path") ||
      inherits(x, "cellcharter_bundle_paths")) return(x)
  if (is.environment(x) || typeof(x) == "externalptr" || inherits(x, "connection")) {
    stop("Unsupported opaque external reference at `", breadcrumb,
         "`; remove it or provide a portable representation before saving.", call. = FALSE)
  }
  if (inherits(x, "cellcharter_autok") &&
      identical(x$checkpoint_mode, "none")) {
    # A no-checkpoint fit is self-contained through embedded_best_artifacts.
    # Its checkpoint_dir was only a required call argument and is not a
    # portable dependency of the completed result.
    x$checkpoint_dir <- NULL
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
        owned = owned, registry = registry,
        path_role = .cc_analysis_child_path_role(x, child_key)
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
        owned = child_owned, registry = registry,
        path_role = .cc_analysis_child_path_role(x, child_key)
      )
    }
    return(x)
  }
  if (base::isS4(x)) {
    for (slot_name in methods::slotNames(x)) {
      packed <- .cc_analysis_pack_paths(
        methods::slot(x, slot_name), directory, slot_name,
        paste0(breadcrumb, "@", slot_name), owned = owned,
        registry = registry, path_role = "auto"
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
  if (identical(path_role, "literal")) return(x)
  if (identical(path_role, "owned")) owned <- TRUE
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
    stop("Portable analysis directories currently support `GRanges` rowRanges, not `",
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

.cc_analysis_unpack_pair <- function(x, n_spots, name) {
  if (!identical(x$format, "SelfHits-v1") || !identical(as.integer(x$nnode),
                                                          as.integer(n_spots)) ||
      any(x$from < 1L | x$from > n_spots | x$to < 1L | x$to > n_spots)) {
    stop("Invalid bundled colPair `", name, "`.", call. = FALSE)
  }
  answer <- S4Vectors::SelfHits(
    from = x$from, to = x$to, nnode = n_spots
  )
  S4Vectors::mcols(answer) <- x$mcols
  S4Vectors::metadata(answer) <- x$metadata
  answer
}

.cc_analysis_save_spe <- function(
    x, directory, prefix, progress, artifact_root = directory,
    registry = .cc_analysis_new_registry()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  row_data <- .cc_analysis_materialize_value(
    SummarizedExperiment::rowData(x), paste0(prefix, "$rowData"), progress
  )
  col_data <- .cc_analysis_materialize_value(
    SummarizedExperiment::colData(x), paste0(prefix, "$colData"), progress
  )
  .cc_analysis_validate_embedded(row_data, paste0(prefix, "$rowData"))
  .cc_analysis_validate_embedded(col_data, paste0(prefix, "$colData"))
  assay_names <- SummarizedExperiment::assayNames(x)
  assays <- lapply(seq_along(assay_names), function(i) {
    .cc_analysis_write_matrix(SummarizedExperiment::assay(x, i, withDimnames = FALSE),
                            directory, paste0(prefix, "_assay_", i), progress)
  })
  names(assays) <- assay_names

  rd_names <- SingleCellExperiment::reducedDimNames(x)
  reduced <- lapply(seq_along(rd_names), function(i) {
    .cc_analysis_write_native_matrix(
      SingleCellExperiment::reducedDim(x, i, withDimnames = FALSE),
      directory, paste0(prefix, "_reduced_", i), progress
    )
  })
  names(reduced) <- rd_names

  # Register this experiment's checkpoint roots before descending into
  # alternative experiments, allowing all children to reuse their artifacts.
  raw_metadata <- .cc_analysis_prepare_import_metadata(
    S4Vectors::metadata(x), artifact_root, progress
  )
  packed_metadata <- .cc_analysis_pack_paths(
    raw_metadata, artifact_root, registry = registry
  )

  alt_names <- SingleCellExperiment::altExpNames(x)
  alternatives <- lapply(seq_along(alt_names), function(i) {
    subdir <- file.path(directory, paste0("alt_", i))
    .cc_analysis_save_spe(SingleCellExperiment::altExp(x, i), subdir,
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
    rowData = row_data,
    rowRanges = if (methods::is(x, "RangedSummarizedExperiment"))
      .cc_analysis_pack_ranges(SummarizedExperiment::rowRanges(x)) else NULL,
    colData = col_data,
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

.cc_analysis_read_reduced <- function(spec, directory, id) {
  if (identical(spec$format, "native-on-load-hdf5-v1")) {
    return(.cc_analysis_realize_matrix(
      .cc_analysis_read_matrix(spec, directory, id), id, progress = FALSE
    ))
  }
  if (!identical(spec$format, "native-matrix-rds-v1")) {
    # Backward compatibility with version-1 analysis directories that stored
    # reduced dimensions as HDF5Array objects.
    return(.cc_analysis_read_matrix(spec, directory, id))
  }
  path <- .cc_analysis_resolve(directory, spec$file, mustWork = FALSE)
  if (!file.exists(path)) {
    stop("Analysis component `", id, "` is missing: ", path, call. = FALSE)
  }
  answer <- readRDS(path)
  if (!is.matrix(answer) || !is.numeric(answer) ||
      !identical(as.integer(dim(answer)), as.integer(spec$dim))) {
    stop("Native reduced dimension `", id, "` is invalid.", call. = FALSE)
  }
  dimnames(answer) <- spec$dimnames
  answer
}

.cc_analysis_load_spe <- function(directory, id = "root", artifact_root = directory) {
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
    value <- .cc_analysis_read_reduced(
      z$reducedDims[[i]], directory,
      paste0(id, "/reducedDim/", names(z$reducedDims)[i])
    )
    if (nrow(value) != ncol(ans)) stop("Reduced dimension row count does not match spots.", call. = FALSE)
    SingleCellExperiment::reducedDim(ans, names(z$reducedDims)[i]) <- value
  }
  for (nm in names(z$colPairs)) {
    SingleCellExperiment::colPair(ans, nm) <-
      .cc_analysis_unpack_pair(z$colPairs[[nm]], ncol(ans), nm)
  }
  for (i in seq_along(z$altExps)) {
    subdir <- .cc_analysis_resolve(directory, z$altExps[[i]]$directory)
    value <- .cc_analysis_load_spe(
      subdir, paste0(id, "/altExp/", names(z$altExps)[i]),
      artifact_root = artifact_root
    )
    if (ncol(value) != ncol(ans)) stop("Alternative experiment spot count does not match.", call. = FALSE)
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
    stop("Analysis directories cannot contain symbolic-link payloads: ",
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
#' Large assays are copied blockwise into analysis-owned HDF5 files. Reduced
#' dimensions, annotations, coordinates, and graphs are saved as native R/
#' Bioconductor objects, so the loaded object keeps only its assays disk-backed.
#' The returned directory can be moved independently of its source H5AD and
#' original backing files. AnnData `/uns` content is copied separately without
#' duplicating the complete source H5AD clone.
#'
#' @param spe A `SpatialExperiment` containing the assays, `rowData`,
#'   `colData`, `spatialCoords`, reduced dimensions, `colPairs`, alternative
#'   experiments, and metadata to preserve.
#' @param output_dir One non-empty destination path. The function creates a
#'   self-contained directory with an object manifest, analysis-owned HDF5
#'   assays, native R payloads, file sizes, and SHA-256 checksums. The directory
#'   must not already exist unless `overwrite = TRUE`.
#' @param overwrite A single logical value. If `TRUE`, replace an existing
#'   destination transactionally after the staged replacement validates. The
#'   existing directory is restored if publication fails.
#' @param progress A single logical value controlling component-, block-,
#'   checksum-, validation-, and publication-progress messages.
#' @return The normalized path to `output_dir`, invisibly. The input `spe` is
#'   not modified. A later [cellcharter_load_analysis()] call reconstructs the
#'   `SpatialExperiment` from this movable directory without the original H5AD
#'   or other source backing files.
#' @export
cellcharter_save_analysis <- function(
    spe, output_dir, overwrite = FALSE, progress = interactive()) {
  if (!methods::is(spe, "SpatialExperiment"))
    stop("`spe` must be a SpatialExperiment.", call. = FALSE)
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
  root_spec <- .cc_analysis_save_spe(spe, file.path(stage, "object"), "root", progress)
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
  test <- .cc_analysis_load_spe(file.path(stage, "object"))
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
#' @param path One existing directory created by
#'   [cellcharter_save_analysis()]. Legacy version-1 directories remain
#'   readable; symbolic-link payloads and paths escaping the directory are
#'   rejected.
#' @param validate A single logical value. If `TRUE`, verify every payload's
#'   manifest-recorded size and SHA-256 checksum before reconstruction. If
#'   `FALSE`, skip payload sizes and checksums; the manifest checksum, safe path
#'   resolution, component dimensions and types, final object dimensions, and
#'   reconstructed object validity are still checked.
#' @return A locally reconstructed `SpatialExperiment`. Assays use
#'   analysis-owned HDF5 files below `path`; `rowData`, `colData`,
#'   `spatialCoords`, reduced dimensions, graphs, alternative experiments, and
#'   metadata are restored in their corresponding Bioconductor compartments.
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
  supported_formats <- c(
    "CellCharterR-analysis-directory",
    # Public releases before the API rename used this format label.  The v1
    # payload schema is otherwise readable by the current loader.
    "CellCharterR-analysis-bundle"
  )
  if (!manifest$format %in% supported_formats)
    stop("Directory is not a CellCharterR analysis directory.", call. = FALSE)
  if (!(as.integer(manifest$format_version) %in% c(1L, .cc_analysis_version)))
    stop("Unsupported analysis-directory format version: ", manifest$format_version, call. = FALSE)
  if (isTRUE(validate)) .cc_analysis_validate_files(manifest, root)
  ans <- .cc_analysis_load_spe(file.path(root, "object"))
  if (!identical(as.integer(dim(ans)), as.integer(manifest$root$dim)))
    stop("Loaded object dimensions disagree with the manifest.", call. = FALSE)
  ans
}
