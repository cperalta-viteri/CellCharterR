#' Build content-addressed checkpoint contracts
#'
#' These helpers keep expensive workflow checkpoints tied to the exact input
#' artifacts, cell order, parameters, implementation, and runtime assumptions
#' under which they may be reused.  Implementation digests intentionally hash
#' selected function definitions and/or explicitly delimited source regions;
#' unrelated edits elsewhere in a source file do not invalidate a stage.
#'
#' @param functions A named list of functions relevant to one stage.
#' @param source_regions A named list whose elements contain `path`, `begin`,
#'   and `end`.  Only source strictly between the unique marker lines is hashed.
#' @return A SHA-256 hexadecimal digest.
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
  if (length(functions) && (is.null(names(functions)) || any(!nzchar(names(functions))))) {
    stop("`functions` must be named.", call. = FALSE)
  }
  if (length(source_regions) &&
      (is.null(names(source_regions)) || any(!nzchar(names(source_regions))))) {
    stop("`source_regions` must be named.", call. = FALSE)
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
#' @param stage One workflow-stage name.
#' @param implementation_digest Digest from
#'   [cellcharter_implementation_digest()].
#' @param upstream Exact upstream artifact hashes and/or immutable data IDs.
#' @param cell_order_sha256 Canonical ordered-cell SHA-256.
#' @param parameters Complete material stage parameters.
#' @param seeds Complete named seeds (an explicitly named `NA_integer_` is
#'   permitted when a stage is deterministic and consumes no seed).
#' @param versions Package and pinned-reference versions.
#' @param runtime Device/runtime fields relevant to reuse; use a named empty
#'   list for device-independent stages.
#' @param schema_version Contract schema version.
#' @return A canonical contract with its own `contract_id`.
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

#' Compare a saved checkpoint contract with an expected contract
#'
#' Verification/migration records store the expected contract separately from
#' their truthful legacy-creation status.  Both created and verified records
#' are reusable only while the full contract and artifact hashes still match.
#'
#' @rdname cellcharter_implementation_digest
#' @param record Saved record or path to an RDS record.
#' @param expected Expected contract.
#' @param artifact_paths Named paths whose current SHA-256 values must match.
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
#' @param expected Expected stage contract.
#' @param artifact_paths Uniquely named existing artifact paths.
#' @param status Either `"created_under_contract"` or
#'   `"verified_existing_checkpoint"`.
#' @param evidence Named validation evidence. Required for legacy verification.
#' @param path Destination RDS path, written atomically.
#' @return The saved record, invisibly.
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
