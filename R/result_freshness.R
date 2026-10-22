# Scientific-input fingerprints for stored downstream results.  These are
# deliberately content based: object dimensions and observation order cannot
# detect an edited annotation, graph weight, or AUCell score.

.cc_result_parameters_sha256 <- function(parameters) {
  digest::digest(
    .cc_contract_canonical(parameters), algo = "sha256", serialize = TRUE
  )
}

.cc_result_vector_sha256 <- function(x, target_block_length = 1000000L) {
  n <- length(x)
  starts <- if (n) seq.int(1L, n, by = target_block_length) else integer()
  block_sha256 <- vapply(starts, function(start) {
    index <- seq.int(start, min(n, start + target_block_length - 1L))
    digest::digest(
      list(index = index, values = as.character(x[index])),
      algo = "sha256", serialize = TRUE
    )
  }, character(1L))
  digest::digest(
    list(
      schema = "CellCharterR_result_vector_v1",
      length = n,
      levels = if (is.factor(x)) levels(x) else NULL,
      ordered = is.ordered(x),
      block_sha256 = unname(block_sha256)
    ),
    algo = "sha256", serialize = TRUE
  )
}

.cc_result_graph_sha256 <- function(spe, graph_name,
                                    target_block_edges = 1000000L) {
  if (!graph_name %in% SingleCellExperiment::colPairNames(spe)) {
    stop("SPE colPairs lacks graph `", graph_name, "`.", call. = FALSE)
  }
  pairs <- SingleCellExperiment::colPair(spe, graph_name)
  if (!identical(S4Vectors::nnode(pairs), ncol(spe))) {
    stop("SPE graph `", graph_name,
         "` has dimensions inconsistent with its observations.", call. = FALSE)
  }
  starts <- if (length(pairs)) {
    seq.int(1L, length(pairs), by = target_block_edges)
  } else {
    integer()
  }
  active_edges <- 0
  block_sha256 <- vapply(starts, function(start) {
    index <- seq.int(start, min(length(pairs), start + target_block_edges - 1L))
    block <- pairs[index]
    metadata <- S4Vectors::mcols(block)
    values <- if (ncol(metadata)) metadata[[1L]] else rep.int(1, length(block))
    if ((!is.numeric(values) && !is.logical(values)) || anyNA(values) ||
        any(!is.finite(values))) {
      stop("SPE graph `", graph_name,
           "` must store finite numeric edge values.", call. = FALSE)
    }
    keep <- values != 0
    active_edges <<- active_edges + sum(keep)
    digest::digest(
      list(
        start = start,
        from = as.integer(S4Vectors::from(block)[keep]),
        to = as.integer(S4Vectors::to(block)[keep])
      ),
      algo = "sha256", serialize = TRUE
    )
  }, character(1L))
  digest::digest(
    list(
      schema = "CellCharterR_binary_graph_v1",
      nnode = S4Vectors::nnode(pairs),
      stored_edges = length(pairs),
      active_edges = active_edges,
      block_sha256 = unname(block_sha256)
    ),
    algo = "sha256", serialize = TRUE
  )
}

# Reproduce the fingerprint written before graph consumers were made explicit
# about their binary-topology semantics. This is verification-only: new result
# records always use .cc_result_graph_sha256().
.cc_result_graph_sha256_legacy <- function(
    spe, graph_name, target_block_edges = 1000000L) {
  if (!graph_name %in% SingleCellExperiment::colPairNames(spe)) {
    stop("SPE colPairs lacks graph `", graph_name, "`.", call. = FALSE)
  }
  pairs <- SingleCellExperiment::colPair(spe, graph_name)
  if (!identical(S4Vectors::nnode(pairs), ncol(spe))) {
    stop("SPE graph `", graph_name,
         "` has dimensions inconsistent with its observations.", call. = FALSE)
  }
  starts <- if (length(pairs)) {
    seq.int(1L, length(pairs), by = target_block_edges)
  } else {
    integer()
  }
  active_edges <- 0
  block_sha256 <- vapply(starts, function(start) {
    index <- seq.int(start, min(length(pairs), start + target_block_edges - 1L))
    block <- pairs[index]
    metadata <- S4Vectors::mcols(block)
    values <- if (ncol(metadata)) metadata[[1L]] else rep.int(1, length(block))
    if ((!is.numeric(values) && !is.logical(values)) || anyNA(values) ||
        any(!is.finite(values))) {
      stop("SPE graph `", graph_name,
           "` must store finite numeric edge values.", call. = FALSE)
    }
    keep <- values != 0
    active_edges <<- active_edges + sum(keep)
    digest::digest(
      list(
        start = start,
        from = as.integer(S4Vectors::from(block)[keep]),
        to = as.integer(S4Vectors::to(block)[keep]),
        weight = as.numeric(values[keep])
      ),
      algo = "sha256", serialize = TRUE
    )
  }, character(1L))
  digest::digest(
    list(
      schema = "CellCharterR_weighted_graph_v1",
      nnode = S4Vectors::nnode(pairs),
      stored_edges = length(pairs),
      active_edges = active_edges,
      block_sha256 = unname(block_sha256)
    ),
    algo = "sha256", serialize = TRUE
  )
}

.cc_result_input_sha256 <- function(spe, spec) {
  if (identical(spec$kind, "colData")) {
    return(.cc_result_vector_sha256(
      .cc_spe_coldata(spe, spec$name, paste0("result input `", spec$name, "`"))
    ))
  }
  if (identical(spec$kind, "domain")) {
    return(.cc_result_vector_sha256(.cc_spe_raw_domain(spe, spec$name)))
  }
  if (identical(spec$kind, "weighted_graph")) {
    return(.cc_result_graph_sha256(spe, spec$name))
  }
  if (identical(spec$kind, "zonation_scores")) {
    return(.cc_zonation_score_sha256(spe, spec$scores, spec$assay_name))
  }
  stop("Unknown CellCharterR result-input kind `", spec$kind, "`.", call. = FALSE)
}

.cc_result_inputs <- function(spe, specs) {
  if (is.null(specs)) return(list())
  if (!is.list(specs) || is.null(names(specs)) || any(!nzchar(names(specs))) ||
      anyDuplicated(names(specs))) {
    stop("Result input specifications must be a uniquely named list.", call. = FALSE)
  }
  lapply(specs, function(spec) {
    fingerprint <- c(spec, list(sha256 = .cc_result_input_sha256(spe, spec)))
    if (identical(spec$kind, "weighted_graph")) {
      fingerprint$semantics <- "binary_topology_v1"
    }
    fingerprint
  })
}

.cc_spe_require_result <- function(spe, result_name, context = "Stored result") {
  identity <- .cc_validate_spe_identity(spe)
  record <- .cc_spe_metadata(spe)$results[[result_name]]
  if (is.null(record) || is.null(record$value)) {
    stop(context, " `", result_name, "` is not available; recompute it.",
         call. = FALSE)
  }
  has_n_cells <- !is.null(record$n_cells)
  has_cell_order <- !is.null(record$cell_order_sha256)
  if (xor(has_n_cells, has_cell_order)) {
    stop(context, " `", result_name,
         "` has incomplete identity metadata; recompute it.", call. = FALSE)
  }
  if (has_n_cells &&
      (!identical(record$n_cells, ncol(spe)) ||
       !identical(record$cell_order_sha256, identity$cell_order_sha256))) {
    stop(context, " `", result_name,
         "` is stale after observation reordering or subsetting; recompute it.",
         call. = FALSE)
  }

  # Records created before result-input fingerprints were introduced retain
  # whatever identity checks their schema recorded.  The earliest package
  # records predate both identity fields, so they remain readable after an
  # upgrade; their scientific freshness cannot be inferred retroactively.
  # Partially populated identity metadata is rejected above rather than being
  # treated as legacy. Newly computed records always enter the stronger branch.
  fingerprints <- record$input_fingerprints
  if (!has_n_cells || is.null(fingerprints) ||
      is.null(record$parameters_sha256)) {
    warning(
      context, " `", result_name,
      "` uses an unverified legacy freshness record; it remains readable for compatibility, but CellCharterR cannot prove that its scientific inputs are unchanged. Recompute it before drawing new conclusions.",
      call. = FALSE
    )
  }
  for (label in names(fingerprints %||% list())) {
    fingerprint <- fingerprints[[label]]
    current <- try(.cc_result_input_sha256(spe, fingerprint), silent = TRUE)
    matches <- !inherits(current, "try-error") &&
      identical(current, fingerprint$sha256)
    if (!matches && identical(fingerprint$kind, "weighted_graph") &&
        is.null(fingerprint$semantics)) {
      legacy <- try(
        .cc_result_graph_sha256_legacy(spe, fingerprint$name), silent = TRUE
      )
      matches <- !inherits(legacy, "try-error") &&
        identical(legacy, fingerprint$sha256)
    }
    if (!matches) {
      stop(context, " `", result_name, "` is stale because its `", label,
           "` scientific input changed; recompute it.", call. = FALSE)
    }
  }
  if (!is.null(record$parameters_sha256) &&
      !identical(record$parameters_sha256,
                 .cc_result_parameters_sha256(record$parameters))) {
    stop(context, " `", result_name,
         "` is stale because its stored analysis parameters changed; recompute it.",
         call. = FALSE)
  }
  record
}
