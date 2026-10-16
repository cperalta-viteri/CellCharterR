# High-level SpatialExperiment API. Numerical kernels use internal `.cc_*`
# names, while these functions keep identity-sensitive workflow state attached
# to one SpatialExperiment (or compatible SingleCellExperiment).

.cc_sce_coldata <- function(sce, name, context = "column") {
  .cc_validate_sce_identity(sce)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`", context, "` must name one colData column.", call. = FALSE)
  }
  if (!name %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("SCE colData lacks requested ", context, " `", name, "`.", call. = FALSE)
  }
  value <- SummarizedExperiment::colData(sce)[[name]]
  if (length(value) != ncol(sce) || anyNA(value)) {
    stop("SCE ", context, " `", name, "` must contain one non-missing value per cell.",
         call. = FALSE)
  }
  value
}

.cc_sce_domain <- function(sce, name) {
  .cc_validate_sce_identity(sce)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`domain_col` must name one cluster-label variable.", call. = FALSE)
  }

  column_data <- SummarizedExperiment::colData(sce)
  if (name %in% colnames(column_data)) {
    value <- column_data[[name]]
  } else {
    reduced <- SingleCellExperiment::reducedDims(sce, withDimnames = FALSE)
    reduced_names <- names(reduced)
    if (name %in% reduced_names) {
      candidate <- reduced[[name]]
      if (length(dim(candidate)) != 2L || ncol(candidate) != 1L) {
        stop("reducedDim `", name, "` must have exactly one cluster-label column.",
             call. = FALSE)
      }
      value <- candidate[, 1L]
    } else {
      hits <- which(vapply(reduced, function(x) {
        !is.null(colnames(x)) && name %in% colnames(x)
      }, logical(1L)))
      if (!length(hits)) {
        stop("Cluster labels `", name,
             "` were not found in `colData(sce)` or `reducedDims(sce)`.",
             call. = FALSE)
      }
      if (length(hits) > 1L) {
        stop("Cluster-label column `", name,
             "` occurs in multiple reducedDims; use an unambiguous name.", call. = FALSE)
      }
      candidate <- reduced[[hits]]
      value <- candidate[, match(name, colnames(candidate))]
    }
  }
  if (length(value) != ncol(sce) || anyNA(value)) {
    stop("Cluster labels `", name,
         "` must contain one non-missing value per observation.", call. = FALSE)
  }
  clustering_stage <- paste0("clustering:", name)
  if (!is.null(.cc_sce_metadata(sce)$stages[[clustering_stage]])) {
    .cc_sce_require_stage(sce, clustering_stage)
  }
  value
}

.cc_sce_reduced_dim <- function(sce, name, context = "representation") {
  .cc_validate_sce_identity(sce, reduced_dims = name)
  value <- SingleCellExperiment::reducedDim(sce, name)
  # `is.numeric()` is FALSE for an S4 DelayedMatrix even when its seed stores
  # doubles or integers. Inspect the delayed storage type without realizing the
  # representation, which is essential for multi-million-spot inputs.
  numeric_value <- is.numeric(value) ||
    (inherits(value, "DelayedArray") &&
       DelayedArray::type(value) %in% c("integer", "double"))
  if (!numeric_value || length(dim(value)) != 2L || nrow(value) != ncol(sce) ||
      ncol(value) == 0L) {
    stop("SCE ", context, " `", name, "` must be a non-empty numeric matrix-like object.",
         call. = FALSE)
  }
  if (is.matrix(value) && any(!is.finite(value))) {
    stop("SCE ", context, " `", name, "` contains non-finite values.", call. = FALSE)
  }
  value
}

.cc_sce_spatial_coordinates <- function(sce, name = "spatial") {
  .cc_validate_sce_identity(sce, reduced_dims = NULL)
  if (inherits(sce, "SpatialExperiment") && identical(name, "spatial")) {
    value <- SpatialExperiment::spatialCoords(sce)
    source <- "spatialCoords"
  } else {
    value <- .cc_sce_reduced_dim(sce, name, "spatial reduced dimension")
    source <- paste0("reducedDim:", name)
  }
  if (!is.numeric(value) || length(dim(value)) != 2L ||
      nrow(value) != ncol(sce) || ncol(value) < 2L) {
    stop("Spatial coordinates must be a numeric cell-by-coordinate matrix with at least two columns.",
         call. = FALSE)
  }
  if (any(!is.finite(as.matrix(value[, seq_len(2L), drop = FALSE])))) {
    stop("Spatial coordinates contain non-finite values.", call. = FALSE)
  }
  attr(value, "cellcharterr_coordinate_source") <- source
  value
}

.cc_sce_metadata <- function(sce) {
  value <- S4Vectors::metadata(sce)$CellCharterR
  if (is.null(value)) value <- list()
  if (is.null(value$stages)) value$stages <- list()
  if (is.null(value$results)) value$results <- list()
  value
}

.cc_sce_vector_sha256 <- function(x, target_block_length = 1000000L) {
  n <- length(x)
  starts <- if (n) seq.int(1L, n, by = target_block_length) else integer()
  block_sha256 <- vapply(starts, function(start) {
    index <- seq.int(start, min(n, start + target_block_length - 1L))
    values <- x[index]
    if (is.factor(values)) values <- as.character(values)
    digest::digest(
      list(index = index, values = unname(values)),
      algo = "sha256", serialize = TRUE
    )
  }, character(1L))
  digest::digest(
    list(schema = "CellCharterR_vector_v1", length = n,
         block_sha256 = unname(block_sha256)),
    algo = "sha256", serialize = TRUE
  )
}

.cc_sce_graph_sha256 <- function(sce, graph_name,
                                 target_block_edges = 1000000L) {
  if (!graph_name %in% SingleCellExperiment::colPairNames(sce)) {
    stop("SCE colPairs lacks graph `", graph_name, "`.", call. = FALSE)
  }
  pairs <- SingleCellExperiment::colPair(sce, graph_name)
  if (!identical(S4Vectors::nnode(pairs), ncol(sce))) {
    stop("SCE graph `", graph_name, "` has dimensions inconsistent with its cells.",
         call. = FALSE)
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
      stop("SCE graph `", graph_name,
           "` must store finite numeric edge values.", call. = FALSE)
    }
    keep <- values != 0
    active_edges <<- active_edges + sum(keep)
    digest::digest(
      list(start = start, from = S4Vectors::from(block)[keep],
           to = S4Vectors::to(block)[keep]),
      algo = "sha256", serialize = TRUE
    )
  }, character(1L))
  digest::digest(
    list(schema = "CellCharterR_binary_graph_v1", nnode = S4Vectors::nnode(pairs),
         stored_edges = length(pairs), active_edges = active_edges,
         block_sha256 = unname(block_sha256)),
    algo = "sha256", serialize = TRUE
  )
}

.cc_sce_stage_input_specs <- function(stage, parameters) {
  usable <- function(specs) {
    Filter(function(spec) {
      is.character(spec$name) && length(spec$name) == 1L &&
        !is.na(spec$name) && nzchar(spec$name)
    }, specs)
  }
  if (startsWith(stage, "graph:")) {
    coordinate_source <- parameters$coordinate_source %||% "reducedDim:spatial"
    return(usable(list(
      coordinates = list(
        kind = "coordinates", name = parameters$spatial,
        source = coordinate_source
      ),
      sample = list(kind = "colData", name = parameters$sample_col),
      graph = list(kind = "graph", name = parameters$graph_name)
    )))
  }
  if (startsWith(stage, "aggregation:")) {
    return(usable(list(
      representation = list(kind = "reducedDim", name = parameters$representation),
      graph = list(kind = "graph", name = parameters$graph_name),
      sample = list(kind = "colData", name = parameters$sample_col),
      output = list(kind = "reducedDim", name = parameters$output)
    )))
  }
  if (startsWith(stage, "clustering:")) {
    return(usable(list(
      representation = list(kind = "reducedDim", name = parameters$representation),
      domain = list(kind = "domain", name = parameters$domain_col)
    )))
  }
  list()
}

.cc_sce_raw_domain <- function(sce, name) {
  column_data <- SummarizedExperiment::colData(sce)
  if (name %in% colnames(column_data)) return(column_data[[name]])
  reduced <- SingleCellExperiment::reducedDims(sce, withDimnames = FALSE)
  if (name %in% names(reduced)) {
    candidate <- reduced[[name]]
    if (length(dim(candidate)) == 2L && ncol(candidate) == 1L) {
      return(candidate[, 1L])
    }
  }
  hits <- which(vapply(reduced, function(x) {
    !is.null(colnames(x)) && name %in% colnames(x)
  }, logical(1L)))
  if (length(hits) == 1L) {
    return(reduced[[hits]][, match(name, colnames(reduced[[hits]]))])
  }
  stop("Cluster labels `", name,
       "` were not found unambiguously in `colData(sce)` or `reducedDims(sce)`.",
       call. = FALSE)
}

.cc_sce_stage_input_sha256 <- function(sce, spec) {
  if (identical(spec$kind, "colData")) {
    return(.cc_sce_vector_sha256(.cc_sce_coldata(sce, spec$name)))
  }
  if (identical(spec$kind, "domain")) {
    return(.cc_sce_vector_sha256(.cc_sce_raw_domain(sce, spec$name)))
  }
  if (identical(spec$kind, "reducedDim")) {
    return(.autok_representation_sha256(.cc_sce_reduced_dim(sce, spec$name)))
  }
  if (identical(spec$kind, "coordinates")) {
    coordinates <- if (identical(spec$source, "spatialCoords")) {
      if (!inherits(sce, "SpatialExperiment")) {
        stop("Native spatial coordinates are no longer available.", call. = FALSE)
      }
      SpatialExperiment::spatialCoords(sce)
    } else {
      .cc_sce_reduced_dim(sce, spec$name, "spatial reduced dimension")
    }
    return(.autok_representation_sha256(coordinates[, seq_len(2L), drop = FALSE]))
  }
  if (identical(spec$kind, "graph")) {
    return(.cc_sce_graph_sha256(sce, spec$name))
  }
  stop("Unknown CellCharterR stage-input kind `", spec$kind, "`.", call. = FALSE)
}

.cc_sce_stage_inputs <- function(sce, specs) {
  lapply(specs, function(spec) c(
    spec, list(sha256 = .cc_sce_stage_input_sha256(sce, spec))
  ))
}

.cc_sce_record_stage <- function(sce, stage, parameters, extra = list(),
                                 input_specs = NULL) {
  info <- .cc_validate_sce_identity(sce)
  if (is.null(input_specs)) input_specs <- .cc_sce_stage_input_specs(stage, parameters)
  metadata <- .cc_sce_metadata(sce)
  metadata$stages[[stage]] <- c(list(
    cell_order_sha256 = info$cell_order_sha256,
    n_cells = ncol(sce),
    parameters = parameters,
    input_fingerprints = .cc_sce_stage_inputs(sce, input_specs)
  ), extra)
  S4Vectors::metadata(sce)$CellCharterR <- metadata
  sce
}

.cc_sce_require_stage <- function(sce, stage) {
  info <- .cc_validate_sce_identity(sce)
  record <- .cc_sce_metadata(sce)$stages[[stage]]
  if (is.null(record)) {
    stop("SCE lacks CellCharterR stage `", stage, "`; run the corresponding SCE adapter first.",
         call. = FALSE)
  }
  if (!identical(record$n_cells, ncol(sce)) ||
      !identical(record$cell_order_sha256, info$cell_order_sha256)) {
    stop("CellCharterR stage `", stage,
         "` is stale after cell reordering or subsetting; recompute it.", call. = FALSE)
  }
  for (label in names(record$input_fingerprints %||% list())) {
    fingerprint <- record$input_fingerprints[[label]]
    current <- try(.cc_sce_stage_input_sha256(sce, fingerprint), silent = TRUE)
    if (inherits(current, "try-error") || !identical(current, fingerprint$sha256)) {
      detail <- if (identical(fingerprint$kind, "reducedDim")) {
        "representation contents"
      } else {
        paste0(fingerprint$kind, " input")
      }
      stop("CellCharterR stage `", stage, "` is stale because its `", label,
           "` ", detail, " changed; recompute it.", call. = FALSE)
    }
  }
  record
}

.cc_sce_graph_matrix <- function(sce, graph_name, require_current = TRUE) {
  if (!graph_name %in% SingleCellExperiment::colPairNames(sce)) {
    stop("SCE colPairs lacks graph `", graph_name, "`.", call. = FALSE)
  }
  stage <- paste0("graph:", graph_name)
  if (require_current && !is.null(.cc_sce_metadata(sce)$stages[[stage]])) {
    .cc_sce_require_stage(sce, stage)
  }
  graph <- SingleCellExperiment::colPair(sce, graph_name, asSparse = TRUE)
  entries <- Matrix::summary(graph)
  if (any(entries$x < 0)) {
    stop("SCE graph `", graph_name, "` contains negative edge weights.",
         call. = FALSE)
  }
  entries <- entries[entries$x != 0, , drop = FALSE]
  graph <- Matrix::sparseMatrix(
    i = entries$i, j = entries$j, x = rep.int(1, nrow(entries)),
    dims = dim(graph)
  )
  if (!identical(dim(graph), c(ncol(sce), ncol(sce)))) {
    stop("SCE graph `", graph_name, "` has dimensions inconsistent with its cells.",
         call. = FALSE)
  }
  graph
}

.cc_sce_assert_sample_local_graph <- function(
    sce, sample, graph_name, target_block_edges = 1000000L) {
  pairs <- SingleCellExperiment::colPair(sce, graph_name)
  starts <- if (length(pairs)) {
    seq.int(1L, length(pairs), by = target_block_edges)
  } else {
    integer()
  }
  crossing_count <- 0L
  first_crossing <- NULL
  for (start in starts) {
    index <- seq.int(start, min(length(pairs), start + target_block_edges - 1L))
    block <- pairs[index]
    metadata <- S4Vectors::mcols(block)
    values <- if (ncol(metadata)) metadata[[1L]] else rep.int(1, length(block))
    keep <- values != 0
    from <- S4Vectors::from(block)[keep]
    to <- S4Vectors::to(block)[keep]
    crossing <- which(sample[from] != sample[to])
    crossing_count <- crossing_count + length(crossing)
    if (is.null(first_crossing) && length(crossing)) {
      first_crossing <- c(from[[crossing[[1L]]]], to[[crossing[[1L]]]])
    }
  }
  if (crossing_count) {
    stop(
      "SCE graph `", graph_name, "` contains ", crossing_count,
      " cross-sample edge(s), including ", first_crossing[[1L]], " (`",
      sample[first_crossing[[1L]]], "`) -> ", first_crossing[[2L]], " (`",
      sample[first_crossing[[2L]]],
      "`); rebuild or correct the graph before aggregation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Build sample-local CellCharter graphs on a SpatialExperiment
#'
#' Delaunay graphs are built independently for every physical sample, then one
#' global long-edge cutoff is calculated across those graph blocks. The final
#' symmetric graph is stored in `colPairs(sce, graph_name)`, so endpoints remain
#' synchronized when the SCE is subset. The recorded stage fingerprints the
#' coordinates, sample assignments, and binary graph topology; changing any of
#' them makes the stage stale and requires graph reconstruction.
#'
#' @param sce A `SingleCellExperiment` or `SpatialExperiment`. Native
#'   `spatialCoords()` are used for a `SpatialExperiment` when `spatial` is the
#'   default `"spatial"`.
#' @param spatial Name of the legacy coordinate matrix in `reducedDims(sce)`.
#'   This is retained for plain `SingleCellExperiment` inputs and non-default
#'   coordinate representations.
#' @param sample_col Name of the physical-sample field in `colData(sce)`.
#' @param graph_name Name used in `colPairs(sce)`.
#' @param percentile Global positive-edge distance percentile to retain.
#' @param progress Whether to display sample-level graph-construction progress.
#' @return The updated `SingleCellExperiment`.
#' @examples
#' example_sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = matrix(
#'     0, nrow = 1, ncol = 6,
#'     dimnames = list("marker1", paste0("cell", seq_len(6)))
#'   )),
#'   colData = S4Vectors::DataFrame(sample = rep("slide1", 6))
#' )
#' SingleCellExperiment::reducedDim(example_sce, "spatial") <- cbind(
#'   x = c(0, 1, 0, 1, 0.5, 1.5), y = c(0, 0, 1, 1, 0.5, 0.5)
#' )
#' example_sce <- cellcharter_graph(example_sce)
#' SingleCellExperiment::colPairNames(example_sce)
#' @export
cellcharter_graph <- function(sce, spatial = "spatial", sample_col = "sample",
                              graph_name = "cellcharter_graph", percentile = 99,
                              progress = interactive()) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  coordinates <- .cc_sce_spatial_coordinates(sce, spatial)
  coordinate_source <- attr(coordinates, "cellcharterr_coordinate_source")
  attr(coordinates, "cellcharterr_coordinate_source") <- NULL
  if (ncol(coordinates) < 2L) stop("Spatial coordinates must have at least two columns.", call. = FALSE)
  sample <- .cc_sce_coldata(sce, sample_col, "sample column")
  sample <- as.character(sample)
  if (any(!nzchar(trimws(sample)))) stop("Sample identifiers must not be blank.", call. = FALSE)
  rows <- split(seq_len(ncol(sce)), factor(sample, levels = unique(sample)))
  coordinates_by_sample <- lapply(
    rows, function(i) coordinates[i, seq_len(2L), drop = FALSE]
  )
  sample_ids <- names(rows)
  bar <- if (progress) utils::txtProgressBar(0L, 2L * length(sample_ids), style = 3L) else NULL
  if (!is.null(bar)) on.exit(close(bar), add = TRUE)
  step <- 0L
  adjacency <- lapply(sample_ids, function(id) {
    value <- tryCatch(
      delaunay_graph(coordinates_by_sample[[id]]),
      error = function(e) stop("Cannot build sample-local graph for `", id, "`: ",
                               conditionMessage(e), call. = FALSE)
    )
    step <<- step + 1L
    if (!is.null(bar)) utils::setTxtProgressBar(bar, step)
    value
  })
  names(adjacency) <- sample_ids
  cutoff <- cellcharter_global_distance_cutoff(adjacency, coordinates_by_sample, percentile)
  adjacency <- lapply(sample_ids, function(id) {
    value <- remove_long_links(
      adjacency[[id]], coordinates_by_sample[[id]],
      percentile = percentile, cutoff = cutoff
    )
    step <<- step + 1L
    if (!is.null(bar)) utils::setTxtProgressBar(bar, step)
    value
  })
  names(adjacency) <- sample_ids

  from <- to <- integer()
  distance <- numeric()
  for (id in names(rows)) {
    entry <- Matrix::summary(adjacency[[id]])
    if (nrow(entry)) {
      global_from <- rows[[id]][entry$i]
      global_to <- rows[[id]][entry$j]
      from <- c(from, global_from)
      to <- c(to, global_to)
      delta <- coordinates[global_from, seq_len(2L), drop = FALSE] -
        coordinates[global_to, seq_len(2L), drop = FALSE]
      distance <- c(distance, sqrt(rowSums(delta^2)))
    }
  }
  pairs <- S4Vectors::SelfHits(from, to, nnode = ncol(sce))
  S4Vectors::mcols(pairs)$distance <- distance
  SingleCellExperiment::colPair(sce, graph_name) <- pairs
  .cc_sce_record_stage(
    sce, paste0("graph:", graph_name),
    list(spatial = spatial, coordinate_source = coordinate_source,
         sample_col = sample_col, graph_name = graph_name,
         percentile = percentile, global_distance_cutoff = cutoff,
         construction = "sample-local Delaunay; globally trimmed"),
    list(sample_sha256 = digest::digest(sample, algo = "sha256"))
  )
}

#' Aggregate CellCharter neighborhoods on a SpatialExperiment
#'
#' Every nonzero graph edge must join observations from the same physical
#' sample. Cross-sample edges are rejected rather than silently discarded. The
#' recorded stage fingerprints the representation, graph topology, sample
#' assignments, and aggregated output using bounded-memory blocks.
#'
#' @param sce A graph-annotated `SingleCellExperiment`.
#' @param representation Input name in `reducedDims(sce)`.
#' @param graph_name Input name in `colPairs(sce)`.
#' @param sample_col Physical-sample field in `colData(sce)`.
#' @param output Output name in `reducedDims(sce)`.
#' @param n_layers Largest exclusive graph hop to aggregate.
#' @param filepath Optional HDF5 output path. The default `NULL` keeps the
#'   aggregated representation in memory and avoids an additional backing-file
#'   dependency. Supply a path only when available RAM is insufficient.
#' @param block_rows HDF5 aggregation row budget.
#' @param progress Whether to display sample- or block-level aggregation progress.
#' @return The updated `SingleCellExperiment`.
#' @examples
#' example_sce <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_sce <- cellcharter_aggregate(
#'   example_sce, output = "example_aggregation", n_layers = 1
#' )
#' dim(SingleCellExperiment::reducedDim(example_sce, "example_aggregation"))
#' @export
cellcharter_aggregate <- function(
    sce, representation = "X_trVAE", graph_name = "cellcharter_graph",
    sample_col = "sample", output = "X_cellcharter", n_layers = 3L,
    filepath = NULL, block_rows = 8192L, progress = interactive()) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  x <- .cc_sce_reduced_dim(sce, representation)
  graph <- .cc_sce_graph_matrix(sce, graph_name)
  sample <- as.character(.cc_sce_coldata(sce, sample_col, "sample column"))
  if (any(!nzchar(trimws(sample)))) {
    stop("Sample identifiers must not be blank.", call. = FALSE)
  }
  .cc_sce_assert_sample_local_graph(sce, sample, graph_name)
  rows <- split(seq_len(ncol(sce)), factor(sample, levels = unique(sample)))
  adjacency <- lapply(rows, function(i) graph[i, i, drop = FALSE])

  if (is.null(filepath)) {
    bar <- if (progress) utils::txtProgressBar(0L, length(rows), style = 3L) else NULL
    if (!is.null(bar)) on.exit(close(bar), add = TRUE)
    pieces <- lapply(seq_along(rows), function(i) {
      value <- .cc_aggregate(
        x[rows[[i]], , drop = FALSE], adjacency[[i]], n_layers
      )
      if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
      value
    })
    names(pieces) <- names(rows)
    value <- matrix(NA_real_, nrow = ncol(sce), ncol = ncol(pieces[[1L]]))
    for (id in names(rows)) value[rows[[id]], ] <- as.matrix(pieces[[id]])
  } else {
    value <- .cc_aggregate_hdf5(
      x, adjacency, rows, filepath = filepath, name = output,
      n_layers = n_layers, block_rows = block_rows, progress = progress
    )
  }
  rownames(value) <- colnames(sce)
  colnames(value) <- paste0("feature_", seq_len(ncol(value)))
  SingleCellExperiment::reducedDim(sce, output) <- value
  .cc_sce_record_stage(
    sce, paste0("aggregation:", output),
    list(representation = representation, graph_name = graph_name,
         sample_col = sample_col, output = output, n_layers = n_layers,
         storage = if (is.null(filepath)) "memory" else "HDF5",
         filepath = if (is.null(filepath)) NULL else normalizePath(filepath))
  )
}

#' Run experimental CellCharter AutoK and attach selected domains to an SCE
#'
#' @section Experimental status:
#' The R Torch GMM and AutoK fitting path is experimental. It is validated on
#' small and tutorial-scale inputs, but can be substantially slower than the
#' Python CellCharter implementation on multi-million-observation datasets.
#' For such datasets, fit domains with the pinned Python implementation and
#' store the resulting labels in `colData(sce)` (or a one-column reducedDim)
#' before using CellCharterR's downstream analyses.
#'
#' @param sce An aggregated `SingleCellExperiment`.
#' @param representation Input name in `reducedDims(sce)`.
#' @param candidate_k Two selectable K bounds or an expanded consecutive list,
#'   as accepted by [cellcharter_autok_fit()].
#' @param checkpoint_dir Durable AutoK checkpoint directory.
#' @param domain_col Output field in `colData(sce)`.
#' @param ... Additional arguments passed to [cellcharter_autok_fit()].
#' @return The updated `SingleCellExperiment`.
#' @examples
#' example_sce <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' demo_fit <- function(x, n_clusters, seed, ...) {
#'   list(k = n_clusters, seed = seed)
#' }
#' demo_predict <- function(model, newdata, ...) {
#'   (seq_len(nrow(newdata)) - 1L) %% model$k
#' }
#' checkpoint <- tempfile("CellCharterR-example-")
#' example_sce <- cellcharter_cluster(
#'   example_sce, representation = "X_cellcharter_tiny",
#'   candidate_k = c(2L, 3L), checkpoint_dir = checkpoint,
#'   domain_col = "example_domain", max_runs = 2L, convergence_tol = 0,
#'   fit_fun = demo_fit, predict_fun = demo_predict,
#'   nll_fun = function(model) as.numeric(model$k)
#' )
#' table(example_sce$example_domain)
#' unlink(checkpoint, recursive = TRUE)
#' @export
cellcharter_cluster <- function(
    sce, representation = "X_cellcharter", candidate_k,
    checkpoint_dir, domain_col = "cellcharter_domain", ...) {
  x <- .cc_sce_reduced_dim(sce, representation)
  aggregation_stage <- .cc_sce_require_stage(
    sce, paste0("aggregation:", representation)
  )
  identity <- .cc_validate_sce_identity(sce, reduced_dims = representation)
  stage_identity <- digest::digest(
    .cc_contract_canonical(list(
      schema = "CellCharterR_sce_autok_input_v1",
      cell_order_sha256 = identity$cell_order_sha256,
      representation = representation,
      aggregation_stage = aggregation_stage
    )),
    algo = "sha256", serialize = TRUE
  )
  fit <- cellcharter_autok_fit(
    x = x, n_clusters = candidate_k, checkpoint_dir = checkpoint_dir,
    data_id = stage_identity, ...
  )
  selected_k <- cellcharter_autok_best_k(fit)
  selected <- fit$best_index[fit$best_index$k == selected_k, , drop = FALSE]
  if (!is.null(fit$embedded_best_artifacts)) {
    artifact <- fit$embedded_best_artifacts[[as.character(selected_k)]]
  } else {
    artifact <- readRDS(selected$path[[1L]])
    if (identical(artifact$format, "cellcharter_autok_round_v1")) {
      artifact <- artifact$artifacts[[as.character(selected_k)]]
    }
  }
  labels <- artifact$labels
  if (length(labels) != ncol(sce) || anyNA(labels)) {
    stop("Selected AutoK labels do not align one-to-one with SCE cells.", call. = FALSE)
  }
  sce[[domain_col]] <- factor(labels, levels = sort(unique(labels)))
  metadata <- .cc_sce_metadata(sce)
  metadata$autok[[domain_col]] <- fit
  S4Vectors::metadata(sce)$CellCharterR <- metadata
  .cc_sce_record_stage(
    sce, paste0("clustering:", domain_col),
    list(representation = representation, candidate_k = candidate_k,
         checkpoint_dir = if (identical(fit$checkpoint_mode, "none")) {
           NULL
         } else {
           normalizePath(checkpoint_dir)
         },
         domain_col = domain_col, selected_k = selected_k,
         representation_sha256 = fit$representation_sha256,
         autok_input_stage_sha256 = stage_identity,
         selected_artifact = if (identical(fit$checkpoint_mode, "none")) {
           NULL
         } else {
           normalizePath(selected$path[[1L]])
         })
  )
}

.cc_sce_store_result <- function(sce, result_name, result, parameters,
                                 input_specs = NULL) {
  metadata <- .cc_sce_metadata(sce)
  identity <- .cc_validate_sce_identity(sce)
  metadata$results[[result_name]] <- list(
    value = result, cell_order_sha256 = identity$cell_order_sha256,
    n_cells = ncol(sce), parameters = parameters,
    parameters_sha256 = .cc_result_parameters_sha256(parameters),
    input_fingerprints = .cc_result_inputs(sce, input_specs)
  )
  S4Vectors::metadata(sce)$CellCharterR <- metadata
  sce
}

#' Attach CellCharter cell-type enrichment to an SCE
#'
#' The stored result fingerprints its group and domain annotations and analysis
#' parameters. Experiment-aware plotting rejects the result if those inputs
#' subsequently change.
#' @param sce A clustered `SingleCellExperiment`.
#' @param group_col Cell-type or group field in `colData(sce)`.
#' @param domain_col Cluster-label variable in `colData(sce)`, the name of a
#'   one-column reducedDim, or an unambiguous column name inside a reducedDim.
#' @param result_name Name under `metadata(sce)$CellCharterR$results`.
#' @param log2_fold_change Return log2 fold change.
#' @param progress Whether to display progress while counting observations.
#' @param chunk_size Maximum number of observations processed per counting
#'   chunk. This bounds temporary memory for very large experiments.
#' @return The updated `SingleCellExperiment`.
#' @examples
#' example_sce <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_sce <- cellcharter_enrichment(
#'   example_sce, group_col = "cell_type_reduced", domain_col = "tiny_domain"
#' )
#' names(S4Vectors::metadata(example_sce)$CellCharterR$results$enrichment$value)
#' @export
cellcharter_enrichment <- function(
    sce, group_col, domain_col = "cellcharter_domain",
    result_name = "enrichment", log2_fold_change = TRUE,
    progress = interactive(), chunk_size = 1000000L) {
  group <- .cc_sce_coldata(sce, group_col, "group column")
  domain <- .cc_sce_domain(sce, domain_col)
  result <- .cc_enrichment(
    group, domain, log2_fold_change,
    progress = progress, chunk_size = chunk_size
  )
  .cc_sce_store_result(
    sce, result_name, result,
    list(group_col = group_col, domain_col = domain_col,
         log2_fold_change = log2_fold_change, progress = progress,
         chunk_size = as.integer(chunk_size)),
    list(
      group = list(kind = "colData", name = group_col),
      domain = list(kind = "domain", name = domain_col)
    )
  )
}

#' Attach CellCharter neighborhood enrichment to an SCE
#'
#' The stored result fingerprints the weighted graph, domain annotations, and
#' analysis parameters so stale results can be detected without materializing
#' the complete graph.
#' @inheritParams cellcharter_enrichment
#' @param graph_name Graph name in `colPairs(sce)`.
#' @param only_inter Exclude within-domain edges.
#' @return The updated `SingleCellExperiment`.
#' @examples
#' example_sce <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_sce <- cellcharter_nhood_enrichment(
#'   example_sce, domain_col = "tiny_domain"
#' )
#' S4Vectors::metadata(example_sce)$CellCharterR$results$nhood_enrichment$value$enrichment
#' @export
cellcharter_nhood_enrichment <- function(
    sce, graph_name = "cellcharter_graph", domain_col = "cellcharter_domain",
    result_name = "nhood_enrichment", log2_fold_change = FALSE,
    only_inter = TRUE) {
  graph <- .cc_sce_graph_matrix(sce, graph_name)
  domain <- .cc_sce_domain(sce, domain_col)
  result <- .cc_nhood_enrichment(graph, domain, log2_fold_change, only_inter)
  .cc_sce_store_result(
    sce, result_name, result,
    list(graph_name = graph_name, domain_col = domain_col,
         log2_fold_change = log2_fold_change, only_inter = only_inter),
    list(
      graph = list(kind = "weighted_graph", name = graph_name),
      domain = list(kind = "domain", name = domain_col)
    )
  )
}

#' Attach differential CellCharter neighborhood enrichment to an SCE
#'
#' The stored result fingerprints the weighted graph, domains, conditions,
#' sample/library assignments, and analysis parameters. Downstream plotting
#' rejects a result after any of those scientific inputs changes.
#' @inheritParams cellcharter_nhood_enrichment
#' @param condition_col Condition field in `colData(sce)`.
#' @param sample_col Library/sample field in `colData(sce)`.
#' @param condition_groups Two conditions to compare.
#' @param n_perms Number of library-level permutations.
#' @param seed Permutation seed.
#' @param permutation_mode CellCharter parity or conventional shuffle mode.
#' @param progress Whether to display progress through the permutations.
#' @return The updated `SingleCellExperiment`.
#' @examples
#' example_sce <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_sce <- cellcharter_diff_nhood_enrichment(
#'   example_sce, domain_col = "tiny_domain",
#'   condition_groups = c("BALBc", "MRL"), n_perms = 5L, seed = 1L
#' )
#' S4Vectors::metadata(example_sce)$CellCharterR$results$diff_nhood_enrichment$value$pvalue
#' @export
cellcharter_diff_nhood_enrichment <- function(
    sce, graph_name = "cellcharter_graph", domain_col = "cellcharter_domain",
    condition_col = "condition", sample_col = "sample",
    result_name = "diff_nhood_enrichment", condition_groups = NULL,
    n_perms = 100L, seed = 0L, log2_fold_change = FALSE, only_inter = TRUE,
    permutation_mode = c("cellcharter_v0.2.0", "shuffle"),
    progress = interactive()) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  if (progress) {
    message("CellCharterR differential neighborhoods: preparing imported graph and library statistics.")
  }
  graph <- .cc_sce_graph_matrix(sce, graph_name)
  domain <- .cc_sce_domain(sce, domain_col)
  condition <- .cc_sce_coldata(sce, condition_col, "condition column")
  sample <- .cc_sce_coldata(sce, sample_col, "sample column")
  if (progress) {
    message("CellCharterR differential neighborhoods: running ",
            format(n_perms, big.mark = ","), " library-level permutations.")
  }
  result <- .cc_diff_nhood_enrichment(
    graph, domain, condition, sample, condition_groups = condition_groups,
    n_perms = n_perms, seed = seed, log2_fold_change = log2_fold_change,
    only_inter = only_inter, permutation_mode = match.arg(permutation_mode),
    progress = progress
  )
  .cc_sce_store_result(
    sce, result_name, result,
    list(graph_name = graph_name, domain_col = domain_col,
         condition_col = condition_col, sample_col = sample_col,
         condition_groups = condition_groups, n_perms = n_perms, seed = seed,
         log2_fold_change = log2_fold_change, only_inter = only_inter,
         permutation_mode = match.arg(permutation_mode), progress = progress),
    list(
      graph = list(kind = "weighted_graph", name = graph_name),
      domain = list(kind = "domain", name = domain_col),
      condition = list(kind = "colData", name = condition_col),
      sample = list(kind = "colData", name = sample_col)
    )
  )
}

#' Run the complete CellCharter SpatialExperiment workflow
#'
#' This convenience function composes the exported SCE adapters. Expensive
#' full datasets should use the staged functions directly so each stage can be
#' checkpointed; this orchestrator is intended for compact examples and tests.
#'
#' @param sce A `SingleCellExperiment`.
#' @param group_col Cell-type/group field.
#' @param candidate_k AutoK candidate range.
#' @param checkpoint_dir Durable AutoK checkpoint directory.
#' @param spatial,representation,sample_col,condition_col,graph_name,aggregated
#'   Named SCE fields used by the workflow.
#' @param domain_col Output domain-label field in `colData(sce)`.
#' @param n_layers Neighborhood aggregation depth.
#' @param condition_groups Differential contrast.
#' @param n_perms Number of differential permutations.
#' @param seed Base seed used by AutoK and permutation analysis.
#' @param autok_args Additional named AutoK arguments, primarily for controlled
#'   testing hooks.
#' @return The updated `SingleCellExperiment`.
#' @examples
#' codex_mouse_spleen_tiny <- cellcharter_mouse_spleen_tiny()
#' demo_fit <- function(x, n_clusters, seed, ...) {
#'   list(k = n_clusters, seed = seed)
#' }
#' demo_predict <- function(model, newdata, ...) {
#'   (seq_len(nrow(newdata)) - 1L) %% model$k
#' }
#' checkpoint <- tempfile("CellCharterR-workflow-")
#' example_sce <- cellcharter_workflow(
#'   codex_mouse_spleen_tiny, group_col = "cell_type_reduced",
#'   candidate_k = c(2L, 3L), checkpoint_dir = checkpoint,
#'   n_layers = 1L, n_perms = 2L, seed = 1L,
#'   condition_groups = c("BALBc", "MRL"),
#'   autok_args = list(
#'     max_runs = 2L, convergence_tol = 0,
#'     fit_fun = demo_fit, predict_fun = demo_predict,
#'     nll_fun = function(model) as.numeric(model$k)
#'   )
#' )
#' example_sce
#' unlink(checkpoint, recursive = TRUE)
#' @export
cellcharter_workflow <- function(
    sce, group_col, candidate_k, checkpoint_dir,
    spatial = "spatial", representation = "X_trVAE", sample_col = "sample",
    condition_col = "condition", graph_name = "cellcharter_graph",
    aggregated = "X_cellcharter", domain_col = "cellcharter_domain",
    n_layers = 3L, condition_groups = NULL, n_perms = 100L, seed = 1L,
    autok_args = list()) {
  sce <- cellcharter_graph(sce, spatial, sample_col, graph_name)
  sce <- cellcharter_aggregate(
    sce, representation, graph_name, sample_col, aggregated, n_layers
  )
  cluster_args <- c(
    list(sce = sce, representation = aggregated, candidate_k = candidate_k,
         checkpoint_dir = checkpoint_dir, domain_col = domain_col, seed = seed),
    autok_args
  )
  sce <- do.call(cellcharter_cluster, cluster_args)
  sce <- cellcharter_enrichment(sce, group_col, domain_col)
  sce <- cellcharter_nhood_enrichment(sce, graph_name, domain_col)
  cellcharter_diff_nhood_enrichment(
    sce, graph_name, domain_col, condition_col, sample_col,
    condition_groups = condition_groups, n_perms = n_perms, seed = seed
  )
}
