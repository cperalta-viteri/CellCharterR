#' Build content-addressed checkpoint contracts
#'
#' These helpers keep expensive workflow checkpoints tied to the exact input
#' artifacts, cell order, parameters, implementation, and runtime assumptions
#' under which they may be reused.  Implementation digests intentionally hash
#' selected function definitions and/or explicitly delimited source regions;
#' unrelated edits elsewhere in a source file do not invalidate a stage.
#'
#' @param functions A uniquely named list of function objects whose formals and
#'   deparsed bodies define material stage behavior. Names and function content
#'   enter the canonical digest; an unnamed or non-function element is rejected.
#' @param source_regions A uniquely named list whose elements contain scalar
#'   `path`, `begin`, and `end` fields. Each path must be readable and each exact
#'   marker line must occur once, with `begin` before `end`. Only source strictly
#'   between those marker lines is hashed, allowing unrelated edits elsewhere.
#' @return One 64-character SHA-256 hexadecimal implementation digest.
#' @examples
#' implementation <- cellcharter_implementation_digest(
#'   functions = list(enrichment = cellcharter_enrichment)
#' )
#' cells <- digest::digest(c("cell1", "cell2"), algo = "sha256")
#' contract <- cellcharter_stage_contract(
#'   stage = "example", implementation_digest = implementation,
#'   upstream = list(input = "example-v1"), cell_order_sha256 = cells,
#'   parameters = list(method = "documented"), seeds = list(seed = NA_integer_),
#'   versions = list(CellCharterR = "0.99.0")
#' )
#' artifact <- tempfile(fileext = ".txt")
#' writeLines("checkpoint", artifact)
#' record_path <- tempfile(fileext = ".rds")
#' cellcharter_record_checkpoint(
#'   contract, c(output = artifact), path = record_path
#' )
#' cellcharter_checkpoint_reusable(
#'   record_path, contract, c(output = artifact)
#' )
#' unlink(c(artifact, record_path))
#' @export
cellcharter_implementation_digest <- function(functions = list(), source_regions = list()) {
  if (!is.list(functions) || !is.list(source_regions)) {
    stop("`functions` and `source_regions` must be lists.", call. = FALSE)
  }
  if (length(functions) && (is.null(names(functions)) ||
      any(!nzchar(names(functions))) || anyDuplicated(names(functions)))) {
    stop("`functions` must be uniquely named.", call. = FALSE)
  }
  if (length(source_regions) &&
      (is.null(names(source_regions)) || any(!nzchar(names(source_regions))) ||
       anyDuplicated(names(source_regions)))) {
    stop("`source_regions` must be uniquely named.", call. = FALSE)
  }
  function_source <- lapply(functions, function(fun) {
    if (!is.function(fun)) stop("Every `functions` element must be a function.", call. = FALSE)
    list(
      formals = paste(deparse(formals(fun), width.cutoff = 500L), collapse = "\n"),
      body = paste(deparse(body(fun), width.cutoff = 500L), collapse = "\n")
    )
  })
  region_source <- lapply(source_regions, function(region) {
    required <- c("path", "begin", "end")
    if (!is.list(region) || !all(required %in% names(region))) {
      stop("Every source region must contain `path`, `begin`, and `end`.", call. = FALSE)
    }
    scalar_string <- vapply(region[required], function(value) {
      is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
    }, logical(1L))
    if (!all(scalar_string)) {
      stop("Every source-region `path`, `begin`, and `end` must be one non-empty string.",
           call. = FALSE)
    }
    lines <- readLines(region$path, warn = FALSE)
    begin <- which(lines == region$begin)
    end <- which(lines == region$end)
    if (length(begin) != 1L || length(end) != 1L || end <= begin) {
      stop("Source-region markers must each occur once and in order: ", region$path,
           call. = FALSE)
    }
    paste(lines[seq.int(begin + 1L, end - 1L)], collapse = "\n")
  })
  digest::digest(
    .cc_contract_canonical(list(functions = function_source, regions = region_source)),
    algo = "sha256", serialize = TRUE
  )
}

.cc_contract_canonical <- function(value) {
  if (is.list(value) && !inherits(value, "data.frame")) {
    if (!is.null(names(value))) value <- value[order(names(value), method = "radix")]
    return(lapply(value, .cc_contract_canonical))
  }
  if (inherits(value, "package_version")) return(as.character(value))
  if (inherits(value, "POSIXt")) {
    return(format(value, tz = "UTC", usetz = TRUE))
  }
  if (is.factor(value)) return(as.character(value))
  value
}

.cc_contract_named_list <- function(value, name) {
  if (!is.list(value) || is.null(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value))) {
    stop("`", name, "` must be a uniquely named list.", call. = FALSE)
  }
  value
}

#' @rdname cellcharter_implementation_digest
#' @param stage One non-empty workflow-stage name.
#' @param implementation_digest Digest from
#'   [cellcharter_implementation_digest()].
#' @param upstream A uniquely named list of exact upstream artifact hashes,
#'   contract IDs, and/or immutable data IDs. Empty names and duplicates are
#'   rejected; values are recursively canonicalized before hashing.
#' @param cell_order_sha256 One 64-character SHA-256 for the canonical ordered
#'   spot identifiers consumed by the stage.
#' @param parameters A uniquely named list of all scientifically material stage
#'   parameters; include values even when they equal defaults.
#' @param seeds A uniquely named list of all seeds (an explicitly named
#'   `NA_integer_` is
#'   permitted when a stage is deterministic and consumes no seed).
#' @param versions A uniquely named list of package, runtime, and pinned
#'   reference versions that constrain safe reuse.
#' @param runtime A uniquely named list of device/runtime fields relevant to
#'   reuse; use `list()` for a device-independent stage.
#' @param schema_version One positive integer contract-schema version.
#' @return A recursively canonicalized contract list containing all supplied
#'   fields plus its own 64-character `contract_id` SHA-256.
#' @export
cellcharter_stage_contract <- function(
    stage, implementation_digest, upstream, cell_order_sha256, parameters,
    seeds, versions, runtime = list(), schema_version = 1L) {
  if (!is.character(stage) || length(stage) != 1L || !nzchar(stage)) {
    stop("`stage` must be one non-empty string.", call. = FALSE)
  }
  if (!is.character(implementation_digest) || length(implementation_digest) != 1L ||
      !grepl("^[0-9a-f]{64}$", implementation_digest)) {
    stop("`implementation_digest` must be one SHA-256 hexadecimal digest.", call. = FALSE)
  }
  if (!is.character(cell_order_sha256) || length(cell_order_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", cell_order_sha256)) {
    stop("`cell_order_sha256` must be one SHA-256 hexadecimal digest.", call. = FALSE)
  }
  schema_version <- as.integer(schema_version)
  if (length(schema_version) != 1L || is.na(schema_version) || schema_version < 1L) {
    stop("`schema_version` must be one positive integer.", call. = FALSE)
  }
  payload <- .cc_contract_canonical(list(
    schema = "CellCharterR_stage_contract",
    schema_version = schema_version,
    stage = stage,
    implementation_digest = implementation_digest,
    upstream = .cc_contract_named_list(upstream, "upstream"),
    cell_order_sha256 = cell_order_sha256,
    parameters = .cc_contract_named_list(parameters, "parameters"),
    seeds = .cc_contract_named_list(seeds, "seeds"),
    versions = .cc_contract_named_list(versions, "versions"),
    runtime = if (!length(runtime)) list() else .cc_contract_named_list(runtime, "runtime")
  ))
  payload$contract_id <- digest::digest(payload, algo = "sha256", serialize = TRUE)
  payload
}

.cc_contract_integral <- function(contract) {
  if (!is.list(contract) || is.null(contract$contract_id)) return(FALSE)
  supplied <- contract$contract_id
  contract$contract_id <- NULL
  identical(
    supplied,
    digest::digest(.cc_contract_canonical(contract), algo = "sha256", serialize = TRUE)
  )
}

.cc_adopt_legacy_implementation_contract <- function(
    record, expected, artifact_paths, allowed_implementation_digests,
    allowed_version_values = list()) {
  if (is.character(record) && length(record) == 1L) {
    if (!file.exists(record)) return(NULL)
    record <- try(readRDS(record), silent = TRUE)
    if (inherits(record, "try-error")) return(NULL)
  }
  if (!is.character(allowed_implementation_digests) ||
      anyNA(allowed_implementation_digests) ||
      any(!grepl("^[0-9a-f]{64}$", allowed_implementation_digests))) {
    stop(
      "`allowed_implementation_digests` must contain SHA-256 hexadecimal digests.",
      call. = FALSE
    )
  }
  if (!is.list(allowed_version_values) ||
      (length(allowed_version_values) &&
       (is.null(names(allowed_version_values)) ||
        any(!nzchar(names(allowed_version_values))) ||
        anyDuplicated(names(allowed_version_values))))) {
    stop("`allowed_version_values` must be a uniquely named list.",
         call. = FALSE)
  }
  if (!is.list(record) || !.cc_contract_integral(expected) ||
      !is.list(record$expected_contract) ||
      !.cc_contract_integral(record$expected_contract) ||
      !record$expected_contract$implementation_digest %in%
        allowed_implementation_digests ||
      !cellcharter_checkpoint_reusable(
        record, record$expected_contract, artifact_paths
      )) {
    return(NULL)
  }

  # A public-symbol/argument rename changes a function digest without changing
  # an already-created artifact. Adoption is deliberately narrow: after
  # replacing only the allowlisted implementation digest, every remaining
  # contract field must be byte-for-byte equivalent to the current expectation.
  normalized <- record$expected_contract
  normalized$implementation_digest <- expected$implementation_digest
  for (name in names(allowed_version_values)) {
    legacy_value <- normalized$versions[[name]]
    if (is.null(legacy_value) ||
        !legacy_value %in% allowed_version_values[[name]] ||
        is.null(expected$versions[[name]])) {
      return(NULL)
    }
    normalized$versions[[name]] <- expected$versions[[name]]
  }
  normalized$contract_id <- NULL
  normalized$contract_id <- digest::digest(
    .cc_contract_canonical(normalized), algo = "sha256", serialize = TRUE
  )
  if (!identical(normalized, expected)) return(NULL)
  record$expected_contract
}

#' Compare a saved checkpoint contract with an expected contract
#'
#' Verification/migration records store the expected contract separately from
#' their truthful legacy-creation status.  Both created and verified records
#' are reusable only while the full contract and artifact hashes still match.
#'
#' @rdname cellcharter_implementation_digest
#' @param record A checkpoint-record list or one path to an existing RDS record.
#'   Missing, unreadable, malformed, or non-integral records return `FALSE`.
#' @param expected An integral contract returned by
#'   [cellcharter_stage_contract()].
#' @param artifact_paths A uniquely named character vector of existing files.
#'   During reuse its names must exactly match the recorded artifacts and every
#'   current SHA-256 must match. During recording, absolute paths, SHA-256
#'   digests, and byte sizes are saved. Files are streamed rather than loaded
#'   into memory and their contents are never modified.
#' @return `TRUE` only for an exact, internally valid match.
#' @export
cellcharter_checkpoint_reusable <- function(record, expected, artifact_paths) {
  if (is.character(record) && length(record) == 1L) {
    if (!file.exists(record)) return(FALSE)
    record <- try(readRDS(record), silent = TRUE)
    if (inherits(record, "try-error")) return(FALSE)
  }
  if (!.cc_contract_integral(expected) || !is.list(record) ||
      !record$status %in% c("created_under_contract", "verified_existing_checkpoint") ||
      !is.list(record$expected_contract) || !.cc_contract_integral(record$expected_contract) ||
      !identical(record$expected_contract$contract_id, expected$contract_id)) {
    return(FALSE)
  }
  if (!is.character(artifact_paths) || is.null(names(artifact_paths)) ||
      any(!nzchar(names(artifact_paths))) || !identical(names(artifact_paths), names(record$artifacts)) ||
      !all(file.exists(artifact_paths))) {
    return(FALSE)
  }
  current <- vapply(artifact_paths, cellcharter_file_sha256, character(1))
  identical(unname(current), unname(vapply(record$artifacts, `[[`, character(1), "sha256")))
}

#' Record creation or read-only verification of checkpoint artifacts
#'
#' @rdname cellcharter_implementation_digest
#' @param expected An integral stage contract returned by
#'   [cellcharter_stage_contract()].
#' @param status Either `"created_under_contract"` or
#'   `"verified_existing_checkpoint"`.
#' @param evidence A named validation-evidence list. It must be non-empty when
#'   `status = "verified_existing_checkpoint"` so provenance does not imply
#'   that a legacy artifact was created under a contract it predates.
#' @param path One destination RDS path. Parent directories are created and the
#'   record is published by an atomic same-directory rename.
#' @return The saved checkpoint-record list, invisibly. Artifact contents are
#'   never modified.
#' @export
cellcharter_record_checkpoint <- function(
    expected, artifact_paths, status = "created_under_contract", evidence = list(), path) {
  if (!.cc_contract_integral(expected)) stop("`expected` is not an integral contract.", call. = FALSE)
  status <- match.arg(status, c("created_under_contract", "verified_existing_checkpoint"))
  if (!is.character(artifact_paths) || is.null(names(artifact_paths)) ||
      any(!nzchar(names(artifact_paths))) || anyDuplicated(names(artifact_paths)) ||
      !all(file.exists(artifact_paths))) {
    stop("`artifact_paths` must be uniquely named existing files.", call. = FALSE)
  }
  if (!is.list(evidence) || (length(evidence) &&
      (is.null(names(evidence)) || any(!nzchar(names(evidence)))))) {
    stop("`evidence` must be a named list.", call. = FALSE)
  }
  if (identical(status, "verified_existing_checkpoint") && !length(evidence)) {
    stop("Legacy checkpoint verification requires explicit evidence.", call. = FALSE)
  }
  artifacts <- Map(function(artifact_path) list(
    path = normalizePath(artifact_path, winslash = "/", mustWork = TRUE),
    sha256 = cellcharter_file_sha256(artifact_path),
    bytes = unname(file.info(artifact_path)$size)
  ), artifact_paths)
  record <- list(
    record_schema = "CellCharterR_checkpoint_record",
    record_schema_version = 1L,
    status = status,
    original_creation_contract = if (identical(status, "created_under_contract")) {
      expected$contract_id
    } else {
      "unknown"
    },
    expected_contract = expected,
    artifacts = artifacts,
    recorded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    evidence = evidence
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), dirname(path), ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(record, temporary)
  if (!file.rename(temporary, path)) stop("Could not atomically publish checkpoint record: ", path)
  invisible(record)
}
