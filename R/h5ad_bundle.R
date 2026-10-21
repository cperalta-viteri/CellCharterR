# Fast, resumable, disk-backed H5AD -> SpatialExperiment conversion for R 4.4.
# The source H5AD is cloned byte-for-byte; no Python or anndataR is required.

.transition_require <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
    if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "),
                              ". Install them with BiocManager::install().", call. = FALSE)
}

.transition_h5ls_cache <- new.env(parent = emptyenv())
.transition_integrity_cache <- new.env(parent = emptyenv())
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

.transition_file_sha256 <- function(path) {
    unname(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

.transition_file_fingerprint <- function(path) {
    info <- file.info(path)
    if (is.na(info$size) || is.na(info$mtime)) {
        stop("Cannot read arrays.h5 size or modification time.", call. = FALSE)
    }
    # Size and mtime survive a same-filesystem rename, while a normal copy gets
    # a new mtime. They are a fast accidental-change detector, not an
    # adversarial integrity proof; `verify = \"full\"` provides the latter.
    list(size = as.numeric(info$size), mtime = as.numeric(info$mtime))
}

.transition_verify_arrays <- function(path, manifest, verify = c("auto", "full", "none"),
                                      manifest_path = NULL,
                                      hash_fun = .transition_file_sha256) {
    verify <- match.arg(verify)
    observed_fingerprint <- .transition_file_fingerprint(path)
    expected_size <- manifest$arrays_size %||% manifest$source$size
    if (!is.null(expected_size) &&
        (!is.finite(observed_fingerprint$size) ||
         !identical(observed_fingerprint$size, as.numeric(expected_size)))) {
        stop("arrays.h5 size differs from its conversion manifest; the bundle is incomplete or corrupted.",
             call. = FALSE)
    }
    expected <- manifest$arrays_sha256
    if (is.null(expected) || !is.character(expected) || length(expected) != 1L ||
        !grepl("^[0-9a-f]{64}$", expected)) {
        if (!identical(verify, "none")) {
            warning("This legacy bundle has no arrays.h5 SHA-256. Re-import it to enable cryptographic integrity verification.",
                    call. = FALSE)
        }
        return(invisible(FALSE))
    }
    if (identical(verify, "none")) return(invisible(NA))
    cache_key <- normalizePath(path, mustWork = TRUE)
    cached <- get0(cache_key, envir = .transition_integrity_cache,
                   inherits = FALSE)
    cache_fingerprint <- c(observed_fingerprint, list(expected = expected))
    if (identical(verify, "auto")) {
        trusted <- manifest$arrays_fingerprint
        if (!is.null(trusted) &&
            identical(as.numeric(trusted$size), observed_fingerprint$size) &&
            identical(as.numeric(trusted$mtime), observed_fingerprint$mtime)) {
            return(invisible(TRUE))
        }
        if (!is.null(cached) &&
            identical(cached$fingerprint, cache_fingerprint)) {
            return(invisible(TRUE))
        }
    }
    message(sprintf("Verifying arrays.h5 SHA-256 (%.2f GiB sequential disk read; assays remain disk-backed)...",
                    observed_fingerprint$size / 1024^3))
    actual <- hash_fun(path)
    if (!identical(actual, expected)) {
        stop("arrays.h5 SHA-256 differs from its conversion manifest; refusing to load a corrupted bundle.",
             call. = FALSE)
    }
    assign(cache_key, list(fingerprint = cache_fingerprint, verified = Sys.time()),
           envir = .transition_integrity_cache)
    if (!is.null(manifest_path)) {
        manifest$arrays_fingerprint <- observed_fingerprint
        manifest$arrays_verified <- Sys.time()
        persisted <- tryCatch({
            .transition_atomic_save_rds(manifest, manifest_path)
            TRUE
        }, error = function(e) {
            message("SHA-256 passed, but the verified file fingerprint could not be persisted: ",
                    conditionMessage(e))
            FALSE
        })
        invisible(persisted)
    }
    invisible(TRUE)
}

.transition_lock_owner_path <- function(lock) file.path(lock, "owner.rds")

.transition_lock_write <- function(lock, owner) {
    .transition_atomic_save_rds(owner, .transition_lock_owner_path(lock))
}

.transition_lock_touch <- function(lock) {
    path <- .transition_lock_owner_path(lock)
    if (!file.exists(path)) return(invisible(FALSE))
    owner <- tryCatch(readRDS(path), error = function(e) NULL)
    expected_id <- attr(lock, "owner_id", exact = TRUE)
    if (is.null(owner) || !identical(owner$pid, Sys.getpid()) ||
        (!is.null(expected_id) && !identical(owner$owner_id, expected_id))) {
        return(invisible(FALSE))
    }
    owner$heartbeat <- Sys.time()
    .transition_lock_write(lock, owner)
    invisible(TRUE)
}

.transition_release_lock <- function(lock) {
    path <- .transition_lock_owner_path(lock)
    owner <- tryCatch(readRDS(path), error = function(e) NULL)
    expected_id <- attr(lock, "owner_id", exact = TRUE)
    if (!is.null(owner) && !is.null(expected_id) &&
        identical(owner$owner_id, expected_id)) {
        unlink(unclass(lock), recursive = TRUE)
        return(invisible(TRUE))
    }
    invisible(FALSE)
}

.transition_pid_alive <- function(pid) {
    is.numeric(pid) && length(pid) == 1L && is.finite(pid) && pid > 0 &&
        isTRUE(suppressWarnings(tryCatch(tools::pskill(as.integer(pid), signal = 0L),
                                        error = function(e) FALSE)))
}

.transition_acquire_lock <- function(output_dir) {
    lock <- file.path(output_dir, ".conversion-lock")
    host <- unname(Sys.info()[["nodename"]])
    if (!dir.create(lock, showWarnings = FALSE)) {
        owner_path <- .transition_lock_owner_path(lock)
        owner <- tryCatch(readRDS(owner_path), error = function(e) NULL)
        heartbeat <- if (!is.null(owner$heartbeat)) owner$heartbeat else
            if (!is.null(owner$created)) owner$created else file.info(lock)$mtime
        age_hours <- as.numeric(difftime(Sys.time(), heartbeat, units = "hours"))
        same_host <- !is.null(owner$host) && identical(as.character(owner$host), host)
        proven_dead <- same_host && !is.null(owner$pid) &&
            !.transition_pid_alive(owner$pid)
        if (!proven_dead) {
            description <- if (is.null(owner)) {
                "owner metadata is missing or unreadable"
            } else {
                sprintf("owner pid=%s host=%s, last heartbeat %.1f hours ago",
                        owner$pid %||% "unknown", owner$host %||% "unknown", age_hours)
            }
            stop("Another conversion may be using this bundle (", description,
                 "). If no conversion is active, rename `.conversion-lock` and retry; it is not deleted without process evidence.",
                 call. = FALSE)
        }
        stale <- file.path(
            output_dir,
            sprintf(".conversion-lock.stale-%s-pid-%s",
                    format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), owner$pid)
        )
        if (!file.rename(lock, stale) || !dir.create(lock, showWarnings = FALSE)) {
            stop("A dead conversion lock was detected but could not be quarantined safely.",
                 call. = FALSE)
        }
        message("Quarantined stale conversion lock at `", stale, "`.")
    }
    created <- Sys.time()
    owner_id <- digest::digest(
        list(pid = Sys.getpid(), host = host, created = created,
             lock = normalizePath(lock, mustWork = FALSE)),
        algo = "sha256"
    )
    owner <- list(pid = Sys.getpid(), host = host, owner_id = owner_id,
                  created = created, heartbeat = created)
    .transition_lock_write(lock, owner)
    attr(lock, "owner_id") <- owner_id
    lock
}

.transition_atomic_save_rds <- function(object, path) {
    temporary <- paste0(path, ".tmp-", Sys.getpid())
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(object, temporary, compress = FALSE)
    if (!file.rename(temporary, path)) stop("Could not publish ", path, call. = FALSE)
    invisible(path)
}

.transition_clone_file <- function(source, storage_dir, chunk_mib = 128L,
                                   lock = NULL,
                                   source_verification = c("sampled", "full")) {
    source_verification <- match.arg(source_verification)
    final <- file.path(storage_dir, "arrays.h5")
    partial <- paste0(final, ".partial")
    state_path <- file.path(storage_dir, "copy-state.rds")
    source_identity <- .transition_source_identity(source)
    source_sha256 <- NULL
    if (identical(source_verification, "full")) {
        message(sprintf(
            "Computing source H5AD SHA-256 before copying (%.2f GiB sequential disk read)...",
            source_identity$size / 1024^3
        ))
        source_sha256 <- .transition_file_sha256(source)
    }
    if (file.exists(final)) {
        manifest_path <- file.path(storage_dir, "manifest.rds")
        if (!file.exists(manifest_path)) {
            if (!file.exists(state_path)) {
                stop("arrays.h5 exists without a manifest or resumable copy state; quarantine it or use a new output directory.",
                     call. = FALSE)
            }
            state <- readRDS(state_path)
            if (!.transition_identity_equal(state$source, source_identity) ||
                unname(file.info(final)$size) != source_identity$size ||
                !identical(.transition_file_samples(final), source_identity$samples)) {
                stop("Unpublished arrays.h5 cannot be proven to match the current source; quarantine it.",
                     call. = FALSE)
            }
            message("Recovering an interrupted manifest publication; computing arrays.h5 SHA-256...")
            manifest <- list(
                format_version = 3L, source = source_identity,
                completed = Sys.time(), copy_seconds = NA_real_,
                arrays_file = "arrays.h5", arrays_size = source_identity$size,
                arrays_sha256 = .transition_file_sha256(final),
                arrays_fingerprint = .transition_file_fingerprint(final)
            )
            .transition_atomic_save_rds(manifest, manifest_path)
            unlink(state_path)
        } else {
            manifest <- readRDS(manifest_path)
        }
        if (!.transition_identity_equal(manifest$source, source_identity)) {
            stop("The bundle belongs to a different or changed source H5AD.", call. = FALSE)
        }
        if (is.null(manifest$arrays_sha256)) {
            message("Upgrading legacy conversion manifest with arrays.h5 SHA-256...")
            manifest$format_version <- max(as.integer(manifest$format_version %||% 1L), 3L)
            manifest$arrays_size <- unname(file.info(final)$size)
            manifest$arrays_sha256 <- .transition_file_sha256(final)
            manifest$arrays_fingerprint <- .transition_file_fingerprint(final)
            .transition_atomic_save_rds(manifest, manifest_path)
        }
        if (identical(source_verification, "full") &&
            !identical(source_sha256, manifest$arrays_sha256)) {
            stop("The full source SHA-256 differs from the published arrays.h5; the source changed or this directory belongs to different input.",
                 call. = FALSE)
        }
        if (identical(source_verification, "full") &&
            !identical(manifest$source_sha256, source_sha256)) {
            manifest$source_sha256 <- source_sha256
            manifest$source_verification <- "full"
            .transition_atomic_save_rds(manifest, manifest_path)
        }
        .transition_verify_arrays(final, manifest, verify = "auto",
                                  manifest_path = manifest_path)
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
        if (!is.null(lock)) .transition_lock_touch(lock)
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
    message("Computing arrays.h5 SHA-256 before publication (one sequential disk read)...")
    arrays_sha256 <- .transition_file_sha256(partial)
    if (!is.null(lock)) .transition_lock_touch(lock)
    if (identical(source_verification, "full") &&
        !identical(arrays_sha256, source_sha256)) {
        quarantine <- file.path(
            storage_dir,
            paste0("arrays.h5.source-mismatch-",
                   format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"))
        )
        if (!file.rename(partial, quarantine)) {
            quarantine <- partial
        }
        stop("The source changed while it was being copied: source and cloned SHA-256 differ. The unpublished copy was retained at `",
             quarantine, "` for inspection.", call. = FALSE)
    }
    if (!file.rename(partial, final)) stop("Could not atomically publish arrays.h5.", call. = FALSE)
    manifest <- list(format_version = 3L, source = source_identity, completed = Sys.time(),
                     copy_seconds = proc.time()[["elapsed"]] - started,
                     arrays_file = "arrays.h5", arrays_size = source_identity$size,
                     arrays_sha256 = arrays_sha256,
                     source_sha256 = source_sha256,
                     source_verification = source_verification,
                     arrays_fingerprint = .transition_file_fingerprint(final))
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

.transition_progress_bar <- function(enabled, total) {
    if (isTRUE(enabled) && total > 0L) utils::txtProgressBar(0L, total, style = 3L) else NULL
}

.transition_close_progress <- function(bar) {
    if (!is.null(bar)) close(bar)
    invisible(NULL)
}

.transition_frame_index <- function(file, group, fallback_path) {
    frame_attrs <- tryCatch(rhdf5::h5readAttributes(file, group),
                            error = function(e) list())
    fallback_name <- basename(fallback_path)
    index_value <- frame_attrs[["_index"]] %||% fallback_name
    if (!is.character(index_value) || length(index_value) != 1L ||
        is.na(index_value)) {
        stop("AnnData dataframe `", group,
             "` has an invalid `_index` attribute.", call. = FALSE)
    }
    index_name <- as.character(index_value)[[1L]]
    if (!nzchar(index_name) || grepl("[/\\\\]", index_name)) {
        stop("AnnData dataframe `", group,
             "` has an invalid `_index` attribute.", call. = FALSE)
    }
    entries <- .transition_children(file, group)
    index_entry <- entries[
        entries$name == index_name & entries$otype == "H5I_DATASET",
        , drop = FALSE
    ]
    if (nrow(index_entry) != 1L) {
        stop("AnnData dataframe `", group, "` index dataset `",
             index_name, "` does not exist or is ambiguous.", call. = FALSE)
    }
    list(
        name = index_name,
        path = paste0(sub("/$", "", group), "/", index_name),
        attributes = frame_attrs,
        entries = entries
    )
}

.transition_read_frame <- function(file, group, n, index_path, columns = NULL,
                                   progress = FALSE) {
    index <- .transition_frame_index(file, group, index_path)
    entries <- index$entries
    available <- setdiff(entries$name, index$name)
    frame_attrs <- index$attributes
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
    if (progress) message("Materialising ", group, " columns into ordinary vectors...")
    bar <- .transition_progress_bar(progress, length(available))
    on.exit(.transition_close_progress(bar), add = TRUE)
    values <- lapply(seq_along(available), function(i) {
        name <- available[[i]]
        path <- paste0(group, "/", name); entry <- entries[match(name, entries$name), , drop = FALSE]
        value <- if (identical(entry$otype, "H5I_GROUP")) {
            children <- .transition_children(file, path)$name
            if (all(c("categories", "codes") %in% children)) {
                .transition_read_categorical(file, path, n)
            } else {
                stop("Unsupported AnnData dataframe encoding at `", path,
                     "`; no column was silently omitted.", call. = FALSE)
            }
        } else {
            value <- rhdf5::h5read(file, sub("^/", "", path), drop = TRUE)
            if (length(value) != n) {
                stop("AnnData dataframe column `", path, "` has ", length(value),
                     " values; expected ", n, ".", call. = FALSE)
            }
            if (!is.null(dim(value))) value <- as.vector(value)
            value
        }
        if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
        value
    })
    names(values) <- available
    ids <- as.character(rhdf5::h5read(file, sub("^/", "", index$path)))
    if (length(ids) != n || anyNA(ids) || anyDuplicated(ids)) {
        stop(index$path, " must contain unique identifiers.", call. = FALSE)
    }
    answer <- if (length(values)) {
        do.call(S4Vectors::DataFrame, c(values, list(check.names = FALSE)))
    } else {
        # DataFrame() otherwise defaults to zero rows, even when the AnnData
        # dataframe legitimately contains only its index.
        S4Vectors::DataFrame(row.names = ids)
    }
    rownames(answer) <- ids
    answer
}

.transition_type_bytes <- function(type) {
    switch(type, raw = 1, logical = 4, integer = 4, double = 8,
           complex = 16, character = 16,
           stop("Unsupported HDF5 value type `", type, "`.", call. = FALSE))
}

.transition_dataset_info <- function(file, path, role, vector = FALSE) {
    .transition_matrix_preflight(file, path, role, "dataset")
    value <- HDF5Array::HDF5Array(file, sub("^/", "", path))
    dimensions <- dim(value)
    if (vector) {
        dimensions <- length(value)
    } else if (length(dimensions) != 2L) {
        stop("`", path, "` must be a two-dimensional HDF5 dataset.", call. = FALSE)
    }
    type <- DelayedArray::type(value)
    if (!type %in% c("integer", "double")) {
        stop("`", path, "` must contain numeric or integer values, not `",
             type, "`.", call. = FALSE)
    }
    list(path = path, dim = as.integer(dimensions), type = type,
         bytes = prod(as.double(dimensions)) * .transition_type_bytes(type))
}

.transition_allocate_matrix <- function(nrow, ncol, type) {
    if (identical(type, "integer")) {
        matrix(NA_integer_, nrow = nrow, ncol = ncol)
    } else {
        matrix(NA_real_, nrow = nrow, ncol = ncol)
    }
}

.transition_read_transposed_matrix <- function(file, info, n_observations,
                                               block_mib = 64, progress = FALSE,
                                               label = info$path) {
    if (!identical(info$dim[[2L]], as.integer(n_observations))) {
        stop("`", info$path, "` has ", info$dim[[2L]],
             " observations; expected ", n_observations, ".", call. = FALSE)
    }
    n_variables <- info$dim[[1L]]
    answer <- .transition_allocate_matrix(n_observations, n_variables, info$type)
    block_rows <- max(1L, floor(block_mib * 1024^2 /
                               max(1, n_variables * .transition_type_bytes(info$type))))
    starts <- seq.int(1L, n_observations, by = block_rows)
    if (progress) {
        message(sprintf("Materialising %s [%s x %s, %.2f GiB] blockwise...",
                        label, n_observations, n_variables, info$bytes / 1024^3))
    }
    bar <- .transition_progress_bar(progress, length(starts))
    on.exit(.transition_close_progress(bar), add = TRUE)
    for (i in seq_along(starts)) {
        index <- starts[[i]]:min(n_observations, starts[[i]] + block_rows - 1L)
        block <- rhdf5::h5read(
            file, sub("^/", "", info$path),
            index = list(seq_len(n_variables), index), drop = FALSE
        )
        answer[index, ] <- t(block)
        rm(block)
        if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
    }
    answer
}

.transition_aucell_info <- function(file, group, n_observations) {
    entries <- .transition_children(file, paste0("/", sub("^/", "", group)))
    attrs <- tryCatch(rhdf5::h5readAttributes(file, sub("^/", "", group)),
                      error = function(e) list())
    index_name <- as.character(attrs[["_index"]] %||% "_index")[[1L]]
    index_entry <- entries[
        entries$name == index_name & entries$otype == "H5I_DATASET",
        , drop = FALSE
    ]
    if (nrow(index_entry) != 1L ||
        suppressWarnings(as.double(index_entry$dim[[1L]])) != n_observations) {
        stop("AUCell dataframe `", group, "` must contain an observation index `",
             index_name, "` of length ", n_observations, ".", call. = FALSE)
    }
    names <- setdiff(entries$name[entries$otype == "H5I_DATASET"], index_name)
    declared_order <- attrs[["column-order"]]
    if (!is.null(declared_order)) {
        declared_order <- as.character(declared_order)
        if (!setequal(declared_order, names)) {
            stop("AUCell dataframe `", group,
                 "` has an invalid column-order attribute.", call. = FALSE)
        }
        names <- declared_order
    }
    infos <- lapply(names, function(name) {
        .transition_dataset_info(file, paste0(group, "/", name),
                                 "AUCell program", vector = TRUE)
    })
    if (!length(infos)) stop("AUCell dataframe contains no program datasets.", call. = FALSE)
    lengths <- vapply(infos, function(x) x$dim, integer(1L))
    if (any(lengths != n_observations)) {
        bad <- names[lengths != n_observations]
        stop("AUCell program(s) have an unexpected observation count: ",
             paste(bad, collapse = ", "), ".", call. = FALSE)
    }
    type <- if (any(vapply(infos, `[[`, character(1L), "type") == "double"))
        "double" else "integer"
    list(names = names, infos = infos, type = type,
         index_path = paste0(group, "/", index_name),
         bytes = as.double(n_observations) * length(names) * .transition_type_bytes(type))
}

.transition_read_aucell <- function(file, info, n_observations,
                                    block_mib = 64, progress = FALSE,
                                    expected_ids = NULL) {
    answer <- .transition_allocate_matrix(n_observations, length(info$names), info$type)
    block_rows <- max(1L, floor(block_mib * 1024^2 /
                               .transition_type_bytes(info$type)))
    starts <- seq.int(1L, n_observations, by = block_rows)
    total <- length(starts) * length(info$names)
    if (progress) {
        message(sprintf("Materialising AUCell [%s x %s, %.2f GiB] blockwise...",
                        n_observations, length(info$names), info$bytes / 1024^3))
    }
    bar <- .transition_progress_bar(progress, total)
    on.exit(.transition_close_progress(bar), add = TRUE)
    if (!is.null(expected_ids)) {
        if (length(expected_ids) != n_observations) {
            stop("Internal AUCell index validation received the wrong number of observation identifiers.",
                 call. = FALSE)
        }
        for (start in starts) {
            index <- start:min(n_observations, start + block_rows - 1L)
            observed_ids <- as.character(rhdf5::h5read(
                file, sub("^/", "", info$index_path),
                index = list(index), drop = TRUE
            ))
            mismatch <- which(observed_ids != expected_ids[index] |
                              is.na(observed_ids) | is.na(expected_ids[index]))
            if (length(mismatch)) {
                position <- index[[mismatch[[1L]]]]
                stop("AUCell dataframe observation index differs from `obs/_index` at position ",
                     position, "; refusing to attach scores to the wrong spots.",
                     call. = FALSE)
            }
        }
    }
    completed <- 0L
    for (column in seq_along(info$infos)) {
        dataset <- info$infos[[column]]$path
        for (start in starts) {
            index <- start:min(n_observations, start + block_rows - 1L)
            answer[index, column] <- rhdf5::h5read(
                file, sub("^/", "", dataset), index = list(index), drop = TRUE
            )
            completed <- completed + 1L
            if (!is.null(bar)) utils::setTxtProgressBar(bar, completed)
        }
    }
    colnames(answer) <- info$names
    answer
}

.transition_h5_object <- function(file, path) {
    path <- paste0("/", sub("^/", "", path))
    parent <- dirname(path)
    if (identical(parent, ".")) parent <- "/"
    entry <- .transition_children(file, parent)
    entry <- entry[entry$name == basename(path), , drop = FALSE]
    if (nrow(entry) != 1L) {
        stop("AnnData object `", path, "` does not exist or is ambiguous.",
             call. = FALSE)
    }
    entry
}

.transition_matrix_preflight <- function(file, path, role,
                                         expected = c("dataset", "dataframe")) {
    expected <- match.arg(expected)
    entry <- .transition_h5_object(file, path)
    if (identical(expected, "dataset") && !identical(entry$otype, "H5I_DATASET")) {
        attrs <- tryCatch(rhdf5::h5readAttributes(file, sub("^/", "", path)),
                          error = function(e) list())
        encoding <- as.character(attrs[["encoding-type"]] %||%
                                     attrs[["encoding_type"]] %||% "group")[1L]
        stop("Unsupported ", role, " layout at `", path, "` (", encoding,
             "). This importer currently requires a dense HDF5 dataset; export a dense representation or select one with the corresponding path argument.",
             call. = FALSE)
    }
    if (identical(expected, "dataframe") && !identical(entry$otype, "H5I_GROUP")) {
        stop("Unsupported ", role, " layout at `", path,
             "`. This importer currently requires an AnnData dataframe group with one HDF5 dataset per program; dense or sparse matrix obsm layouts are not silently converted.",
             call. = FALSE)
    }
    if (identical(expected, "dataframe")) {
        attrs <- tryCatch(rhdf5::h5readAttributes(file, sub("^/", "", path)),
                          error = function(e) list())
        encoding <- as.character(attrs[["encoding-type"]] %||%
                                     attrs[["encoding_type"]] %||% "")[1L]
        if (encoding %in% c("csr_matrix", "csc_matrix")) {
            stop("Unsupported sparse ", role, " layout at `", path, "` (",
                 encoding, "). Export it as an AnnData dataframe group before import.",
                 call. = FALSE)
        }
        children <- .transition_children(file, paste0("/", sub("^/", "", path)))
        children <- children[children$name != "_index", , drop = FALSE]
        if (!nrow(children) || any(children$otype != "H5I_DATASET")) {
            stop("Unsupported ", role, " dataframe at `", path,
                 "`: expected at least one direct HDF5 dataset and no nested encodings.",
                 call. = FALSE)
        }
    }
    invisible(entry)
}

.transition_validate_sample_ids <- function(sample_id, sample_column) {
    invalid <- is.na(sample_id) | !nzchar(trimws(sample_id))
    if (any(invalid)) {
        examples <- utils::head(which(invalid), 5L)
        stop("`", sample_column,
             "` contains missing or blank sample identifiers at observation position(s) ",
             paste(examples, collapse = ", "),
             ". Complete physical-sample identifiers are required before auditing graph isolation.",
             call. = FALSE)
    }
    invisible(sample_id)
}

.transition_character_estimate <- function(
        file, path, n,
        bytes_per_value = getOption(
            "CellCharterR.import.character_bytes_per_value", 256
        )) {
    if (!n) return(0)
    if (!is.numeric(bytes_per_value) || length(bytes_per_value) != 1L ||
        is.na(bytes_per_value) || !is.finite(bytes_per_value) ||
        bytes_per_value < 1) {
        stop(
            "`CellCharterR.import.character_bytes_per_value` must be one positive finite number.",
            call. = FALSE
        )
    }
    # Do not sample variable-length HDF5 strings here. rhdf5/HDF5 can segfault
    # when the same large vlen-string dataset is read through successive
    # hyperslabs (observed with the 3.5-million-spot DeepSpot-M index). The
    # actual reader performs one complete read later. This configurable
    # allowance covers one CHARSXP/pointer plus 256 payload bytes per value by
    # default, which is deliberately above ordinary spot IDs and annotations.
    as.double(n) * (64 + as.double(bytes_per_value))
}

.transition_frame_estimate <- function(file, group, n, columns = NULL) {
    fallback <- paste0(sub("/$", "", group), "/_index")
    index <- .transition_frame_index(file, group, fallback)
    entries <- index$entries
    available <- setdiff(entries$name, index$name)
    if (!is.null(columns)) {
        unknown <- setdiff(columns, available)
        if (length(unknown)) {
            stop("Unknown ", group, " columns: ", paste(unknown, collapse = ", "),
                 call. = FALSE)
        }
        available <- columns
    }
    bytes <- vapply(available, function(name) {
        entry <- entries[match(name, entries$name), , drop = FALSE]
        path <- paste0(sub("/$", "", group), "/", name)
        if (identical(entry$otype, "H5I_GROUP")) {
            categories_path <- paste0(path, "/categories")
            categories <- tryCatch(
                as.character(rhdf5::h5read(file, sub("^/", "", categories_path))),
                error = function(e) character()
            )
            category_bytes <- sum(64 + nchar(categories, type = "bytes"))
            return(as.double(n) * 4 + category_bytes)
        }
        # HDF5Array probes a dataset's first value to infer its R type. That
        # probe fails for valid AnnData ENUM/boolean and some string datasets,
        # even though rhdf5 can read them correctly. The h5ls class is enough
        # for a conservative pre-allocation estimate and avoids making the RAM
        # safety check less compatible than the actual reader.
        if (identical(as.character(entry$dclass), "STRING")) {
            return(.transition_character_estimate(file, path, n))
        }
        bytes_per_value <- switch(
            as.character(entry$dclass),
            FLOAT = 8, INTEGER = 4, ENUM = 4, STRING = 16,
            BITFIELD = 4, OPAQUE = 16,
            stop("Unsupported AnnData dataframe HDF5 class `",
                 entry$dclass, "` at `", group, "/", name, "`.",
                 call. = FALSE)
        )
        as.double(n) * bytes_per_value
    }, numeric(1L))
    index_bytes <- .transition_character_estimate(file, index$path, n)
    sum(bytes) + index_bytes
}

.transition_graph_estimate <- function(file) {
    graphs <- .transition_children(file, "/obsp")
    groups <- graphs$name[graphs$otype == "H5I_GROUP"]
    nnz <- vapply(groups, function(name) {
        path <- paste0("obsp/", name, "/data")
        tryCatch(as.double(length(HDF5Array::HDF5Array(file, path))),
                 error = function(e) 0)
    }, numeric(1L))
    names(nnz) <- groups
    # Retained SelfHits: integer from/to vectors plus numeric edge values.
    list(bytes = sum(nnz) * 16, nnz = sum(nnz), per_graph_nnz = nnz,
         graphs = length(groups))
}

.transition_materialization_preflight <- function(
        file, n_features, n_spots, obs_columns, var_columns,
        spatial_info, reduced_info, aucell_info,
        block_mib, max_materialized_gib) {
    metadata_bytes <- .transition_frame_estimate(file, "/obs", n_spots, obs_columns) +
        .transition_frame_estimate(file, "/var", n_features, var_columns) +
        as.double(n_spots) * 8 # sample_id vector/column pointer estimate
    embedding_bytes <- sum(vapply(reduced_info, `[[`, numeric(1L), "bytes")) +
        if (is.null(aucell_info)) 0 else aucell_info$bytes
    graph <- .transition_graph_estimate(file)
    components <- c(
        metadata = metadata_bytes,
        spatial_coordinates = spatial_info$bytes,
        reduced_dimensions = embedding_bytes,
        graphs = graph$bytes
    )
    total <- sum(components)
    # R's vector headers, DataFrame/list scaffolding, character payloads and
    # allocator fragmentation are data dependent. A 2x envelope is deliberately
    # conservative (and covered the measured DeepSpot-M retained object), while
    # the bounded graph/transpose workspace is added explicitly.
    overhead_factor <- 2
    graph_workspace <- if (graph$graphs) {
        4 * (as.double(n_spots) + 1) + 2 * block_mib * 1024^2
    } else {
        0
    }
    transient <- max(2 * block_mib * 1024^2, graph_workspace)
    estimated_peak <- total * overhead_factor + transient
    message(sprintf(
        paste0("Canonical import preflight: approximately %.2f GiB raw native payload ",
               "and %.2f GiB conservative peak-RAM allowance for metadata/coordinates/",
               "embeddings/graphs; count assay remains HDF5-backed. Matrix and graph ",
               "imports use bounded blocks."),
        total / 1024^3, estimated_peak / 1024^3
    ))
    if (is.finite(max_materialized_gib) &&
        estimated_peak > max_materialized_gib * 1024^3) {
        stop(sprintf(
            paste0("Estimated native peak materialization (%.2f GiB) exceeds ",
                   "`max_materialized_gib = %.2f`. Increase the limit only after ",
                   "confirming available RAM."),
            estimated_peak / 1024^3, max_materialized_gib
        ), call. = FALSE)
    }
    list(
        retained_bytes = total, retained_gib = total / 1024^3,
        estimated_peak_bytes = estimated_peak,
        estimated_peak_gib = estimated_peak / 1024^3,
        components_bytes = components, graph_nnz = graph$nnz,
        graph_nnz_by_name = graph$per_graph_nnz,
        block_mib = block_mib, estimated_peak_block_mib = transient / 1024^2,
        r_overhead_safety_factor = overhead_factor,
        excludes = "memory held by unrelated objects in the calling R session"
    )
}

.transition_read_obsp <- function(file, group, n_spots, block_mib = 64,
                                  progress = FALSE) {
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
    indptr <- as.integer(rhdf5::h5read(file, paste0(group, "/indptr")))
    if (anyNA(indptr) || length(indptr) != n_spots + 1L ||
        indptr[1L] != 0L || any(diff(indptr) < 0L)) {
        stop("Invalid AnnData sparse matrix arrays at `", group, "`.", call. = FALSE)
    }
    nnz <- utils::tail(indptr, 1L)
    data_length <- tryCatch(
        as.double(length(HDF5Array::HDF5Array(
            file, paste0(sub("^/", "", group), "/data")
        ))), error = function(e) NA_real_
    )
    index_length <- tryCatch(
        as.double(length(HDF5Array::HDF5Array(
            file, paste0(sub("^/", "", group), "/indices")
        ))), error = function(e) NA_real_
    )
    if (!is.finite(nnz) || nnz < 0 || nnz != data_length ||
        nnz != index_length) {
        stop("Invalid AnnData sparse matrix arrays at `", group, "`.", call. = FALSE)
    }
    from <- integer(nnz)
    to <- integer(nnz)
    weights <- numeric(nnz)
    block_edges <- max(1L, floor(block_mib * 1024^2 / 16))
    starts <- if (nnz) seq.int(1, nnz, by = block_edges) else numeric()
    if (progress && nnz) {
        message(sprintf("Materialising graph `%s` [%s edges] blockwise...",
                        basename(group), format(nnz, big.mark = ",")))
    }
    bar <- .transition_progress_bar(progress, length(starts))
    on.exit(.transition_close_progress(bar), add = TRUE)
    value_min <- Inf
    value_max <- -Inf
    explicit_zeros <- 0
    diagonal_entries <- 0
    for (i in seq_along(starts)) {
        start <- starts[[i]]
        end <- min(nnz, start + block_edges - 1)
        edge_index <- seq.int(start, end)
        indices <- as.integer(rhdf5::h5read(
            file, paste0(group, "/indices"), index = list(edge_index),
            drop = TRUE
        ))
        block_weights <- as.numeric(rhdf5::h5read(
            file, paste0(group, "/data"), index = list(edge_index),
            drop = TRUE
        ))
        if (length(indices) != length(edge_index) ||
            length(block_weights) != length(edge_index) || anyNA(indices) ||
            any(indices < 0L | indices >= n_spots) ||
            any(!is.finite(block_weights))) {
            stop("Invalid AnnData sparse matrix arrays at `", group, "`.",
                 call. = FALSE)
        }
        # indptr is zero-based. findInterval maps every edge offset to its
        # CSR row / CSC column without ever expanding all major indices at once.
        major <- findInterval(edge_index - 1, indptr, all.inside = TRUE)
        minor <- indices + 1L
        block_from <- if (identical(encoding, "csr_matrix")) major else minor
        block_to <- if (identical(encoding, "csr_matrix")) minor else major
        from[edge_index] <- block_from
        to[edge_index] <- block_to
        weights[edge_index] <- block_weights
        value_min <- min(value_min, block_weights)
        value_max <- max(value_max, block_weights)
        explicit_zeros <- explicit_zeros + sum(block_weights == 0)
        diagonal_entries <- diagonal_entries + sum(block_from == block_to)
        if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
    }
    answer <- S4Vectors::SelfHits(from = from, to = to, nnode = n_spots,
                                  value = weights, sort.by.query = FALSE)
    S4Vectors::metadata(answer)$anndata_obsp <- list(
        encoding_type = encoding,
        encoding_version = as.character(attrs[["encoding-version"]])[1L],
        shape = as.integer(shape), nnz = length(weights),
        value_min = if (length(weights)) value_min else NA_real_,
        value_max = if (length(weights)) value_max else NA_real_,
        explicit_zeros = explicit_zeros,
        diagonal_entries = diagonal_entries
    )
    answer
}

.transition_build_spe <- function(arrays_file, assay_path, assay_name, obs_columns, var_columns,
                                  spatial_path, reduced_dims, aucell_path, sample_column,
                                  materialize_block_mib = 64,
                                  max_materialized_gib = Inf,
                                  progress = interactive()) {
    if (!is.numeric(materialize_block_mib) || length(materialize_block_mib) != 1L ||
        !is.finite(materialize_block_mib) || materialize_block_mib <= 0) {
        stop("`materialize_block_mib` must be one positive finite number.", call. = FALSE)
    }
    if (!is.numeric(max_materialized_gib) || length(max_materialized_gib) != 1L ||
        is.na(max_materialized_gib) || max_materialized_gib <= 0) {
        stop("`max_materialized_gib` must be one positive number or `Inf`.", call. = FALSE)
    }
    if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
        stop("`progress` must be TRUE or FALSE.", call. = FALSE)
    }
    .transition_matrix_preflight(arrays_file, assay_path, "assay", "dataset")
    .transition_matrix_preflight(arrays_file, spatial_path, "spatial obsm", "dataset")
    for (path in unname(reduced_dims)) {
        .transition_matrix_preflight(arrays_file, path, "reduced-dimension obsm", "dataset")
    }
    if (!is.null(aucell_path)) {
        .transition_matrix_preflight(arrays_file, aucell_path, "AUCell obsm", "dataframe")
    }
    assay <- HDF5Array::HDF5Array(arrays_file, assay_path)
    n_features <- nrow(assay); n_spots <- ncol(assay)
    spatial_info <- .transition_dataset_info(
        arrays_file, spatial_path, "spatial obsm"
    )
    if (!identical(spatial_info$dim, c(2L, n_spots))) {
        stop("Expected a 2-by-N spatial matrix.", call. = FALSE)
    }
    reduced_info <- lapply(unname(reduced_dims), function(path) {
        .transition_dataset_info(arrays_file, path, "reduced-dimension obsm")
    })
    names(reduced_info) <- names(reduced_dims)
    for (name in names(reduced_info)) {
        if (!identical(reduced_info[[name]]$dim[[2L]], as.integer(n_spots))) {
            stop("Reduced dimension `", name, "` has an unexpected observation count.",
                 call. = FALSE)
        }
    }
    aucell_info <- if (is.null(aucell_path)) NULL else
        .transition_aucell_info(arrays_file, aucell_path, n_spots)
    memory_plan <- .transition_materialization_preflight(
        arrays_file, n_features, n_spots, obs_columns, var_columns,
        spatial_info, reduced_info, aucell_info,
        materialize_block_mib, max_materialized_gib
    )
    col_data <- .transition_read_frame(
        arrays_file, "/obs", n_spots, "obs/_index", obs_columns,
        progress = progress
    )
    row_data <- .transition_read_frame(
        arrays_file, "/var", n_features, "var/_index", var_columns,
        progress = progress
    )
    if (!is.null(sample_column) && !sample_column %in% colnames(col_data)) {
        stop("sample_column must be imported and present; use NULL only for known single-sample data.", call. = FALSE)
    }
    rownames(assay) <- rownames(row_data); colnames(assay) <- rownames(col_data)
    spatial <- .transition_read_transposed_matrix(
        arrays_file, spatial_info, n_spots, materialize_block_mib,
        progress, "spatial coordinates"
    )
    colnames(spatial) <- c("x", "y"); rownames(spatial) <- rownames(col_data)
    sample_id <- if (is.null(sample_column)) rep.int("sample01", n_spots) else as.character(col_data[[sample_column]])
    .transition_validate_sample_ids(sample_id, sample_column %||% "sample_id")
    spe <- SpatialExperiment::SpatialExperiment(
        assays = S4Vectors::SimpleList(stats::setNames(list(assay), assay_name)), rowData = row_data,
        colData = col_data, spatialCoords = spatial, sample_id = sample_id)
    for (name in names(reduced_dims)) {
        value <- .transition_read_transposed_matrix(
            arrays_file, reduced_info[[name]], n_spots, materialize_block_mib,
            progress, paste0("reducedDim `", name, "`")
        )
        rownames(value) <- colnames(spe); SingleCellExperiment::reducedDim(spe, name) <- value
    }
    if (!is.null(aucell_path)) {
        aucell <- .transition_read_aucell(
            arrays_file, aucell_info, n_spots, materialize_block_mib, progress,
            expected_ids = rownames(col_data)
        )
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
            arrays_file, paste0("/obsp/", name), n_spots,
            block_mib = materialize_block_mib, progress = progress
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
        uns_entries = uns$name,
        canonical_storage = list(
            assay = "HDF5-backed",
            colData = "ordinary vectors", rowData = "ordinary vectors",
            reducedDims = "ordinary matrices", spatialCoords = "ordinary matrix",
            colPairs = "SelfHits", uns = "ordinary locator metadata"
        ),
        materialization_preflight = memory_plan)
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
        chunk_mib = 128L, overwrite = FALSE,
        source_verification = c("sampled", "full"),
        materialize_block_mib = 64,
        max_materialized_gib = getOption("CellCharterR.import.max_materialized_gib", Inf),
        progress = interactive()) {
    .transition_require(c("DelayedArray", "HDF5Array", "rhdf5", "S4Arrays", "S4Vectors",
                          "SingleCellExperiment", "SpatialExperiment", "SummarizedExperiment"))
    source_verification <- match.arg(source_verification)
    h5ad_file <- normalizePath(h5ad_file, mustWork = TRUE)
    # Inspect an existing destination before acquiring the conversion lock.
    # Lock acquisition writes owner metadata (and may quarantine a stale lock),
    # so an overwrite rejection must happen first to remain strictly read-only.
    if (dir.exists(output_dir)) {
        output_dir <- normalizePath(output_dir, mustWork = TRUE)
        rds_path <- file.path(output_dir, "spe.rds")
        if (file.exists(rds_path) && !overwrite) {
            stop("spe.rds exists; set overwrite = TRUE.", call. = FALSE)
        }
    } else {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        output_dir <- normalizePath(output_dir, mustWork = TRUE)
        rds_path <- file.path(output_dir, "spe.rds")
    }
    lock <- .transition_acquire_lock(output_dir)
    on.exit(.transition_release_lock(lock), add = TRUE)
    if (file.exists(file.path(output_dir, "spe_backing.h5"))) {
        stop("Legacy/incomplete spe_backing.h5 found; quarantine it before conversion.", call. = FALSE)
    }
    arrays_file <- .transition_clone_file(
        h5ad_file, output_dir, chunk_mib, lock,
        source_verification = source_verification
    )
    .transition_lock_touch(lock)
    spe <- .transition_build_spe(arrays_file, assay_path, assay_name, obs_columns, var_columns,
                                 spatial_path, reduced_dims, aucell_path, sample_column,
                                 materialize_block_mib, max_materialized_gib,
                                 progress)
    .transition_lock_touch(lock)
    # Keep the reconstruction recipe base-R-only. Serialized Bioconductor S4
    # objects are not reliably portable across releases (for example, newer
    # rowRanges can embed the separately packaged Seqinfo class).
    manifest_path <- file.path(output_dir, "manifest.rds")
    manifest <- readRDS(manifest_path)
    manifest$format_version <- 3L
    manifest$import_spec <- list(
        assay_path = assay_path, assay_name = assay_name,
        obs_columns = obs_columns, var_columns = var_columns,
        spatial_path = spatial_path, reduced_dims = reduced_dims,
        aucell_path = aucell_path, sample_column = sample_column
    )
    .transition_atomic_save_rds(manifest, manifest_path)
    S4Vectors::metadata(spe)$transitionR$bundle_dir <- output_dir
    .transition_atomic_save_rds(spe, rds_path); spe
}

# Load a moved SPE bundle by rebuilding its S4 shell against local packages.
.transition_validate_import_spec <- function(spec, context = "import specification") {
    required <- c(
        "assay_path", "assay_name", "obs_columns", "var_columns",
        "spatial_path", "reduced_dims", "aucell_path", "sample_column"
    )
    if (!is.list(spec) || is.null(names(spec)) ||
        !setequal(names(spec), required) || anyDuplicated(names(spec))) {
        stop("`", context, "` must contain exactly: ",
             paste(required, collapse = ", "), ".", call. = FALSE)
    }
    scalar_path <- function(x, name, nullable = FALSE) {
        if (nullable && is.null(x)) return(invisible(NULL))
        if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
            stop("`", context, "$", name,
                 "` must be one non-empty HDF5 path",
                 if (nullable) " or NULL" else "", ".", call. = FALSE)
        }
    }
    scalar_path(spec$assay_path, "assay_path")
    scalar_path(spec$spatial_path, "spatial_path")
    scalar_path(spec$aucell_path, "aucell_path", nullable = TRUE)
    if (!is.character(spec$assay_name) || length(spec$assay_name) != 1L ||
        is.na(spec$assay_name) || !nzchar(spec$assay_name)) {
        stop("`", context, "$assay_name` must be one non-empty string.",
             call. = FALSE)
    }
    for (name in c("obs_columns", "var_columns")) {
        value <- spec[[name]]
        if (!is.null(value) &&
            (!is.character(value) || anyNA(value) || any(!nzchar(value)) ||
             anyDuplicated(value))) {
            stop("`", context, "$", name,
                 "` must be NULL or unique non-empty column names.", call. = FALSE)
        }
    }
    if (!is.character(spec$reduced_dims) || is.null(names(spec$reduced_dims)) ||
        anyNA(spec$reduced_dims) || any(!nzchar(spec$reduced_dims)) ||
        anyNA(names(spec$reduced_dims)) || any(!nzchar(names(spec$reduced_dims))) ||
        anyDuplicated(names(spec$reduced_dims))) {
        stop("`", context,
             "$reduced_dims` must be a uniquely named character vector of HDF5 paths.",
             call. = FALSE)
    }
    if (!is.null(spec$sample_column) &&
        (!is.character(spec$sample_column) || length(spec$sample_column) != 1L ||
         is.na(spec$sample_column) || !nzchar(spec$sample_column))) {
        stop("`", context,
             "$sample_column` must be one non-empty column name or NULL.",
             call. = FALSE)
    }
    spec[required]
}

load_transition_spe <- function(
        storage_dir, verify = c("auto", "full", "none"),
        legacy_import_spec = NULL,
        materialize_block_mib = 64,
        max_materialized_gib = getOption("CellCharterR.import.max_materialized_gib", Inf),
        progress = interactive()) {
    verify <- match.arg(verify)
    storage_dir <- normalizePath(storage_dir, mustWork = TRUE)
    manifest_path <- file.path(storage_dir, "manifest.rds")
    if (!file.exists(manifest_path)) {
        stop("Bundle lacks manifest.rds and cannot be reconstructed portably.", call. = FALSE)
    }
    manifest <- readRDS(manifest_path)
    version <- suppressWarnings(as.integer(manifest$format_version %||% 1L))
    if (length(version) != 1L || is.na(version) || !version %in% 1:3) {
        stop("Unsupported SPE import-directory format version: ",
             manifest$format_version %||% "missing", ".", call. = FALSE)
    }
    arrays_name <- manifest$arrays_file %||% "arrays.h5"
    if (!is.character(arrays_name) || length(arrays_name) != 1L ||
        is.na(arrays_name) || !nzchar(arrays_name) ||
        grepl("[/\\\\]", arrays_name) || arrays_name %in% c(".", "..")) {
        stop("The SPE import manifest contains an unsafe arrays-file name.",
             call. = FALSE)
    }
    arrays_file <- normalizePath(
        file.path(storage_dir, arrays_name), mustWork = TRUE
    )
    .transition_verify_arrays(arrays_file, manifest, verify = verify,
                              manifest_path = manifest_path)
    spec <- manifest$import_spec
    if (is.null(spec) && version != 1L) {
        stop("SPE import-directory format ", version,
             " manifest lacks the required `import_spec`; refusing to guess reconstruction parameters.",
             call. = FALSE)
    }
    if (is.null(spec)) {
        # Format-1 bundles were created before reconstruction parameters were
        # recorded. They used these public defaults; rebuilding avoids loading
        # their release-specific serialized S4 skeleton.
        default_spec <- list(
            assay_path = "X", assay_name = "expression",
            obs_columns = NULL, var_columns = NULL,
            spatial_path = "obsm/spatial",
            reduced_dims = c(PCA = "obsm/X_pca", cellcharter = "obsm/X_cellcharter"),
            aucell_path = "obsm/score_aucell", sample_column = "slide_id_concat"
        )
        if (is.null(legacy_import_spec)) {
            message(
                "Legacy format-1 import did not record converter options; ",
                "reconstructing with the historical defaults. If the original ",
                "conversion used custom paths or column subsets, pass ",
                "`legacy_import_spec` explicitly."
            )
            spec <- default_spec
        } else {
            spec <- legacy_import_spec
        }
    } else if (!is.null(legacy_import_spec)) {
        stop("`legacy_import_spec` is only valid for a format-1 import that lacks `import_spec`.",
             call. = FALSE)
    }
    spec <- .transition_validate_import_spec(spec)
    spe <- do.call(
        .transition_build_spe,
        c(list(
            arrays_file = arrays_file,
            materialize_block_mib = materialize_block_mib,
            max_materialized_gib = max_materialized_gib,
            progress = progress
        ), spec)
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
#' on disk. A SHA-256 of the independent `arrays.h5` artifact is recorded before
#' publication. Python, reticulate, anndataR, and zellkonverter are not required.
#'
#' The disk-backed importer deliberately accepts dense HDF5 datasets for the
#' assay, spatial coordinates, and reduced dimensions. Sparse `X` and sparse or
#' matrix-form AUCell `obsm` entries fail with an actionable error instead of
#' being silently materialised or omitted.
#'
#' The count assay is the only HDF5-backed slot in the returned object.
#' `colData`, `rowData`, spatial coordinates, reduced dimensions (including
#' AUCell), imported `SelfHits` graphs, and `uns` locator metadata are native R
#' objects. The importer reports both raw-payload and conservative peak-memory
#' estimates before allocating these slots and reads matrices and graphs in
#' bounded blocks without creating a second full transposed matrix.
#'
#' @param h5ad_file One existing source H5AD path. The file is read directly
#'   through HDF5; Python and AnnData are not required.
#' @param output_dir One destination directory. The importer writes an
#'   independent `arrays.h5` clone, a reconstructed `spe.rds` convenience
#'   object, and `manifest.rds`. Move the complete directory when transporting
#'   the result. The returned object does not depend on `h5ad_file`.
#' @param assay_path Non-empty HDF5 path to the feature-by-spot assay, by
#'   default AnnData `X`. It is exposed without materialisation as
#'   `assay(spe, assay_name)` backed by `output_dir/arrays.h5`.
#' @param assay_name Non-empty name assigned to the imported assay in
#'   `assayNames(spe)`.
#' @param obs_columns `NULL` to import every `obs` column, or a unique character
#'   vector selecting columns stored in `colData(spe)`. The AnnData observation
#'   index becomes `colnames(spe)` and `rownames(colData(spe))`.
#' @param var_columns `NULL` to import every `var` column, or a unique character
#'   vector selecting columns stored in `rowData(spe)`. The AnnData variable
#'   index becomes `rownames(spe)` and `rownames(rowData(spe))`.
#' @param spatial_path Non-empty HDF5 path to a dense AnnData spatial matrix.
#'   Through R's HDF5 dimension convention it must be two-by-spots; the importer
#'   transposes it blockwise to spots-by-two and stores it in
#'   `spatialCoords(spe)`.
#' @param reduced_dims A uniquely named character vector mapping each desired
#'   `reducedDim` name to a dense AnnData `obsm` HDF5 path. Through R's HDF5
#'   dimension convention each dataset must be features-by-spots; it is
#'   transposed blockwise and stored as an ordinary spots-by-features matrix in
#'   `reducedDim(spe, name)`. Use `character()` to import none.
#' @param aucell_path Non-empty HDF5 path to an AnnData dataframe-encoded AUCell
#'   `obsm` entry. Its named program columns are materialised blockwise as the
#'   ordinary spots-by-program matrix `reducedDim(spe, "AUCell")`. Use `NULL`
#'   to omit it; matrix-encoded and sparse AUCell entries are rejected.
#' @param sample_column One imported `obs` column used as the
#'   `SpatialExperiment` sample identifier. Use `NULL` only for known
#'   single-sample data; all spots are then assigned `sample01`. When non-`NULL`,
#'   it must be selected by `obs_columns` if a subset is requested.
#' @param chunk_mib Positive copy-buffer target in MiB used while cloning the
#'   H5AD to `arrays.h5`. It controls I/O workspace, not the HDF5 assay's final
#'   size or the memory retained by the returned object.
#' @param overwrite A single logical value. If `FALSE`, an existing published
#'   `spe.rds` is rejected. If `TRUE`, the destination is replaced under an
#'   ownership lock after the independent copy and reconstructed object pass
#'   validation.
#' @param source_verification Source-copy verification mode. `"sampled"`
#'   validates sampled source regions while copying and assumes the source does
#'   not change concurrently. `"full"` additionally hashes the entire source
#'   and requires it to match `arrays.h5` before publication.
#' @param materialize_block_mib Positive target size in MiB for temporary blocks
#'   used to materialise metadata, coordinates, reduced dimensions, AUCell, and
#'   sparse graphs. It does not materialise the count assay.
#' @param max_materialized_gib Maximum allowed conservative GiB estimate for
#'   native metadata, coordinates, embeddings, graphs, object overhead, and
#'   bounded workspaces. The disk-backed assay is excluded. Defaults to option
#'   `CellCharterR.import.max_materialized_gib` or `Inf`; imports exceeding the
#'   limit stop before materialisation.
#' @param progress A single logical value controlling copy, materialisation,
#'   graph-import, validation, and publication progress messages.
#' @return A valid disk-backed `SpatialExperiment`. The assay is stored in
#'   `assay(spe, assay_name)`; `obs` and `var` annotations are in `colData` and
#'   `rowData`; coordinates are in `spatialCoords`; requested `obsm` matrices
#'   are in `reducedDims`; sparse `obsp` matrices are imported as weighted
#'   `SelfHits` in `colPairs`; and an `/uns` locator plus conversion provenance
#'   is stored in `metadata(spe)$transitionR`. The independent directory also
#'   contains `arrays.h5`, `spe.rds`, and `manifest.rds`.
#' @export
#' @examples
#' \dontrun{
#' spe <- cellcharter_import_h5ad("large.h5ad", "large-spe")
#' }
cellcharter_import_h5ad <- function(
        h5ad_file, output_dir, assay_path = "X", assay_name = "expression",
        obs_columns = NULL, var_columns = NULL, spatial_path = "obsm/spatial",
        reduced_dims = c(PCA = "obsm/X_pca", cellcharter = "obsm/X_cellcharter"),
        aucell_path = "obsm/score_aucell", sample_column = "slide_id_concat",
        chunk_mib = 128L, overwrite = FALSE,
        source_verification = c("sampled", "full"),
        materialize_block_mib = 64,
        max_materialized_gib = getOption("CellCharterR.import.max_materialized_gib", Inf),
        progress = interactive()) {
    read_transition_h5ad_spe(
        h5ad_file = h5ad_file, output_dir = output_dir,
        assay_path = assay_path, assay_name = assay_name,
        obs_columns = obs_columns, var_columns = var_columns,
        spatial_path = spatial_path, reduced_dims = reduced_dims,
        aucell_path = aucell_path, sample_column = sample_column,
        chunk_mib = chunk_mib, overwrite = overwrite,
        source_verification = source_verification,
        materialize_block_mib = materialize_block_mib,
        max_materialized_gib = max_materialized_gib, progress = progress
    )
}

#' Load a movable disk-backed SpatialExperiment
#'
#' Reconstructs the S4 experiment shell with the locally installed
#' Bioconductor classes and attaches it to the saved `arrays.h5`. This
#' avoids cross-release serialization dependencies such as `Seqinfo`. Use this
#' function after moving the complete directory.
#'
#' @param path One existing directory created by
#'   [cellcharter_import_h5ad()]. Move `arrays.h5`, `manifest.rds`, and
#'   `spe.rds` together; reconstruction uses the manifest and local
#'   `arrays.h5`, not the original H5AD or the serialized convenience object.
#' @param verify Integrity verification mode. `"auto"` (the default) trusts the
#'   size and modification-time fingerprint recorded when the SHA-256 was
#'   computed, and re-hashes after that fingerprint changes (for example after
#'   a normal file copy). `"full"` always performs a complete sequential read
#'   and is the mode for adversarial cryptographic verification. `"none"`
#'   skips cryptographic verification. Even in `"none"` mode, file size, HDF5
#'   structure, and object validity are checked during reconstruction.
#' @param legacy_import_spec Optional named list describing the original
#'   converter arguments for a format-1 directory, whose manifest predates
#'   recorded import recipes. It must contain `assay_path`, `assay_name`,
#'   `obs_columns`, `var_columns`, `spatial_path`, `reduced_dims`,
#'   `aucell_path`, and `sample_column`. Omit it for historical defaults.
#' @param materialize_block_mib Maximum target size of each temporary block used
#'   while reading coordinates and embeddings into native matrices.
#' @param max_materialized_gib Optional safety limit for the conservative peak
#'   estimate for native metadata, coordinates, embeddings, graphs, R-object
#'   overhead, and bounded import workspaces. The count assay is excluded
#'   because it remains HDF5-backed. Defaults to the
#'   `CellCharterR.import.max_materialized_gib` option or `Inf`.
#' @param progress A single logical value controlling column-, matrix-block-,
#'   graph-, validation-, and reconstruction-progress messages.
#' @return A valid disk-backed `SpatialExperiment`. The assay points to
#'   `path/arrays.h5`; native `colData`, `rowData`, `spatialCoords`,
#'   `reducedDims`, and weighted `colPairs` are rebuilt from the recorded import
#'   specification; and the `/uns` locator and reconstruction provenance are in
#'   `metadata(spe)$transitionR`.
#' @export
#' @examples
#' \dontrun{
#' spe <- cellcharter_load_spe("large-spe")
#' }
cellcharter_load_spe <- function(
        path, verify = c("auto", "full", "none"),
        legacy_import_spec = NULL,
        materialize_block_mib = 64,
        max_materialized_gib = getOption("CellCharterR.import.max_materialized_gib", Inf),
        progress = interactive()) {
    load_transition_spe(
        path, verify = verify, legacy_import_spec = legacy_import_spec,
        materialize_block_mib = materialize_block_mib,
        max_materialized_gib = max_materialized_gib, progress = progress
    )
}
