# High-level SpatialExperiment API. Numerical kernels use internal `.cc_*`
# names, while these functions keep identity-sensitive workflow state attached
# to one SpatialExperiment.

.cc_spe_coldata <- function(spe, name, context = "column") {
  .cc_validate_spe_identity(spe)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`", context, "` must name one colData column.", call. = FALSE)
  }
  if (!name %in% colnames(SummarizedExperiment::colData(spe))) {
    stop("SPE colData lacks requested ", context, " `", name, "`.", call. = FALSE)
  }
  value <- SummarizedExperiment::colData(spe)[[name]]
  if (length(value) != ncol(spe) || anyNA(value)) {
    stop("SPE ", context, " `", name, "` must contain one non-missing value per spot.",
         call. = FALSE)
  }
  value
}

.cc_spe_domain <- function(spe, name) {
  .cc_validate_spe_identity(spe)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`domain_col` must name one cluster-label variable.", call. = FALSE)
  }

  column_data <- SummarizedExperiment::colData(spe)
  if (name %in% colnames(column_data)) {
    value <- column_data[[name]]
  } else {
    reduced <- SingleCellExperiment::reducedDims(spe, withDimnames = FALSE)
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
             "` were not found in `colData(spe)` or `reducedDims(spe)`.",
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
  if (length(value) != ncol(spe) || anyNA(value)) {
    stop("Cluster labels `", name,
         "` must contain one non-missing value per observation.", call. = FALSE)
  }
  # Domain labels are valid downstream inputs independently of how they were
  # produced. A caller may intentionally replace a previously fitted label
  # column (for example, to select another already-fitted AutoK solution).
  # Each downstream result fingerprints the labels it actually consumes, so a
  # superseded clustering-stage record must not block that analysis.
  value
}

.cc_spe_reduced_dim <- function(spe, name, context = "representation") {
  .cc_validate_spe_identity(spe, reduced_dims = name)
  value <- SingleCellExperiment::reducedDim(spe, name)
  # `is.numeric()` is FALSE for an S4 DelayedMatrix even when its seed stores
  # doubles or integers. Inspect the delayed storage type without realizing the
  # representation, which is essential for multi-million-spot inputs.
  numeric_value <- is.numeric(value) ||
    (inherits(value, "DelayedArray") &&
       DelayedArray::type(value) %in% c("integer", "double"))
  if (!numeric_value || length(dim(value)) != 2L || nrow(value) != ncol(spe) ||
      ncol(value) == 0L) {
    stop("SPE ", context, " `", name, "` must be a non-empty numeric matrix-like object.",
         call. = FALSE)
  }
  if (is.matrix(value) && any(!is.finite(value))) {
    stop("SPE ", context, " `", name, "` contains non-finite values.", call. = FALSE)
  }
  value
}

.cc_spe_spatial_coordinates <- function(spe, name = "spatial") {
  .cc_validate_spe_identity(spe, reduced_dims = NULL)
  if (inherits(spe, "SpatialExperiment") && identical(name, "spatial")) {
    value <- SpatialExperiment::spatialCoords(spe)
    source <- "spatialCoords"
  } else {
    value <- .cc_spe_reduced_dim(spe, name, "spatial reduced dimension")
    source <- paste0("reducedDim:", name)
  }
  if (!is.numeric(value) || length(dim(value)) != 2L ||
      nrow(value) != ncol(spe) || ncol(value) < 2L) {
    stop("Spatial coordinates must be a numeric spot-by-coordinate matrix with at least two columns.",
         call. = FALSE)
  }
  if (any(!is.finite(as.matrix(value[, seq_len(2L), drop = FALSE])))) {
    stop("Spatial coordinates contain non-finite values.", call. = FALSE)
  }
  attr(value, "cellcharterr_coordinate_source") <- source
  value
}

.cc_spe_metadata <- function(spe) {
  value <- S4Vectors::metadata(spe)$CellCharterR
  if (is.null(value)) value <- list()
  if (is.null(value$stages)) value$stages <- list()
  if (is.null(value$results)) value$results <- list()
  value
}

.cc_spe_vector_sha256 <- function(x, target_block_length = 1000000L) {
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

.cc_spe_graph_sha256 <- function(spe, graph_name,
                                 target_block_edges = 1000000L) {
  if (!graph_name %in% SingleCellExperiment::colPairNames(spe)) {
    stop("SPE colPairs lacks graph `", graph_name, "`.", call. = FALSE)
  }
  pairs <- SingleCellExperiment::colPair(spe, graph_name)
  if (!identical(S4Vectors::nnode(pairs), ncol(spe))) {
    stop("SPE graph `", graph_name, "` has dimensions inconsistent with its spots.",
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
      stop("SPE graph `", graph_name,
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

.cc_spe_stage_input_specs <- function(stage, parameters) {
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

.cc_spe_raw_domain <- function(spe, name) {
  column_data <- SummarizedExperiment::colData(spe)
  if (name %in% colnames(column_data)) return(column_data[[name]])
  reduced <- SingleCellExperiment::reducedDims(spe, withDimnames = FALSE)
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
       "` were not found unambiguously in `colData(spe)` or `reducedDims(spe)`.",
       call. = FALSE)
}

.cc_spe_stage_input_sha256 <- function(spe, spec) {
  if (identical(spec$kind, "colData")) {
    return(.cc_spe_vector_sha256(.cc_spe_coldata(spe, spec$name)))
  }
  if (identical(spec$kind, "domain")) {
    return(.cc_spe_vector_sha256(.cc_spe_raw_domain(spe, spec$name)))
  }
  if (identical(spec$kind, "reducedDim")) {
    return(.autok_representation_sha256(.cc_spe_reduced_dim(spe, spec$name)))
  }
  if (identical(spec$kind, "coordinates")) {
    coordinates <- if (identical(spec$source, "spatialCoords")) {
      if (!inherits(spe, "SpatialExperiment")) {
        stop("Native spatial coordinates are no longer available.", call. = FALSE)
      }
      SpatialExperiment::spatialCoords(spe)
    } else {
      .cc_spe_reduced_dim(spe, spec$name, "spatial reduced dimension")
    }
    return(.autok_representation_sha256(coordinates[, seq_len(2L), drop = FALSE]))
  }
  if (identical(spec$kind, "graph")) {
    return(.cc_spe_graph_sha256(spe, spec$name))
  }
  stop("Unknown CellCharterR stage-input kind `", spec$kind, "`.", call. = FALSE)
}

.cc_spe_stage_inputs <- function(spe, specs) {
  lapply(specs, function(spec) c(
    spec, list(sha256 = .cc_spe_stage_input_sha256(spe, spec))
  ))
}

.cc_spe_record_stage <- function(spe, stage, parameters, extra = list(),
                                 input_specs = NULL) {
  info <- .cc_validate_spe_identity(spe)
  if (is.null(input_specs)) input_specs <- .cc_spe_stage_input_specs(stage, parameters)
  metadata <- .cc_spe_metadata(spe)
  metadata$stages[[stage]] <- c(list(
    cell_order_sha256 = info$cell_order_sha256,
    n_cells = ncol(spe),
    parameters = parameters,
    input_fingerprints = .cc_spe_stage_inputs(spe, input_specs)
  ), extra)
  S4Vectors::metadata(spe)$CellCharterR <- metadata
  spe
}

.cc_spe_require_stage <- function(spe, stage) {
  info <- .cc_validate_spe_identity(spe)
  record <- .cc_spe_metadata(spe)$stages[[stage]]
  if (is.null(record)) {
    stop("SPE lacks CellCharterR stage `", stage, "`; run the corresponding SPE adapter first.",
         call. = FALSE)
  }
  if (!identical(record$n_cells, ncol(spe)) ||
      !identical(record$cell_order_sha256, info$cell_order_sha256)) {
    stop("CellCharterR stage `", stage,
         "` is stale after spot reordering or subsetting; recompute it.", call. = FALSE)
  }
  for (label in names(record$input_fingerprints %||% list())) {
    fingerprint <- record$input_fingerprints[[label]]
    current <- try(.cc_spe_stage_input_sha256(spe, fingerprint), silent = TRUE)
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

.cc_spe_graph_matrix <- function(spe, graph_name, require_current = FALSE) {
  if (!graph_name %in% SingleCellExperiment::colPairNames(spe)) {
    stop("SPE colPairs lacks graph `", graph_name, "`.", call. = FALSE)
  }
  stage <- paste0("graph:", graph_name)
  # A graph stored in colPairs() is a first-class scientific input. The
  # SingleCellExperiment subsetting machinery remaps SelfHits endpoints to the
  # current columns, so an older construction-stage fingerprint must not make
  # an otherwise valid imported, reordered, or subset graph unusable. Callers
  # that specifically audit construction provenance can still opt in with
  # require_current = TRUE.
  if (require_current && !is.null(.cc_spe_metadata(spe)$stages[[stage]])) {
    .cc_spe_require_stage(spe, stage)
  }
  graph <- SingleCellExperiment::colPair(spe, graph_name, asSparse = TRUE)
  entries <- Matrix::summary(graph)
  if (any(entries$x < 0)) {
    stop("SPE graph `", graph_name, "` contains negative edge weights.",
         call. = FALSE)
  }
  entries <- entries[entries$x != 0, , drop = FALSE]
  graph <- Matrix::sparseMatrix(
    i = entries$i, j = entries$j, x = rep.int(1, nrow(entries)),
    dims = dim(graph)
  )
  if (!identical(dim(graph), c(ncol(spe), ncol(spe)))) {
    stop("SPE graph `", graph_name, "` has dimensions inconsistent with its spots.",
         call. = FALSE)
  }
  graph
}

.cc_spe_assert_sample_local_graph <- function(
    spe, sample, graph_name, target_block_edges = 1000000L) {
  pairs <- SingleCellExperiment::colPair(spe, graph_name)
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
      "SPE graph `", graph_name, "` contains ", crossing_count,
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
#' symmetric graph is stored in `colPairs(spe, graph_name)`, so endpoints remain
#' synchronized when the SPE is subset. The recorded stage fingerprints the
#' coordinates, sample assignments, and binary graph topology for optional
#' provenance auditing; it does not prevent use of the current graph artifact.
#'
#' @param spe A `SpatialExperiment` with one column per spatial spot.
#' @param spatial Coordinate source. The default `"spatial"` reads the first
#'   two columns of `spatialCoords(spe)`. Any other value names a spot-by-
#'   coordinate matrix in `reducedDims(spe)`.
#' @param sample_col Name of the non-missing physical-sample identifier in
#'   `colData(spe)`. Graphs are constructed independently within each value;
#'   no edge may connect different samples.
#' @param graph_name Name under which the symmetric graph is written to
#'   `colPairs(spe)`. Its `SelfHits` metadata column `distance` stores the
#'   Euclidean length of every retained directed edge. Provenance is written to
#'   `metadata(spe)$CellCharterR$stages[[paste0("graph:", graph_name)]]`.
#' @param percentile Numeric percentile from 0 to 100 used to calculate one
#'   global cutoff over all positive sample-local Delaunay edge lengths. Edges
#'   with distance less than or equal to the cutoff are retained.
#' @param progress A single logical value indicating whether to display the two
#'   sample-level phases: Delaunay construction and long-edge removal.
#' @return The updated `SpatialExperiment`, with the graph in
#'   `colPairs(spe, graph_name)`.
#' @examples
#' example_spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = matrix(
#'     0, nrow = 1, ncol = 6,
#'     dimnames = list("marker1", paste0("spot", seq_len(6)))
#'   )),
#'   colData = S4Vectors::DataFrame(sample = rep("slide1", 6)),
#'   spatialCoords = matrix(
#'     c(0, 1, 0, 1, 0.5, 1.5, 0, 0, 1, 1, 0.5, 0.5),
#'     ncol = 2,
#'     dimnames = list(paste0("spot", seq_len(6)), c("x", "y"))
#'   )
#' )
#' example_spe <- cellcharter_graph(example_spe)
#' SingleCellExperiment::colPairNames(example_spe)
#' @export
cellcharter_graph <- function(spe, spatial = "spatial", sample_col = "sample",
                              graph_name = "cellcharter_graph", percentile = 99,
                              progress = interactive()) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  coordinates <- .cc_spe_spatial_coordinates(spe, spatial)
  coordinate_source <- attr(coordinates, "cellcharterr_coordinate_source")
  attr(coordinates, "cellcharterr_coordinate_source") <- NULL
  if (ncol(coordinates) < 2L) stop("Spatial coordinates must have at least two columns.", call. = FALSE)
  sample <- .cc_spe_coldata(spe, sample_col, "sample column")
  sample <- as.character(sample)
  if (any(!nzchar(trimws(sample)))) stop("Sample identifiers must not be blank.", call. = FALSE)
  rows <- split(seq_len(ncol(spe)), factor(sample, levels = unique(sample)))
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
  pairs <- S4Vectors::SelfHits(from, to, nnode = ncol(spe))
  S4Vectors::mcols(pairs)$distance <- distance
  SingleCellExperiment::colPair(spe, graph_name) <- pairs
  .cc_spe_record_stage(
    spe, paste0("graph:", graph_name),
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
#' @param spe A graph-annotated `SpatialExperiment` with a numeric input
#'   representation and sample-local graph.
#' @param representation Name of the numeric spot-by-feature input in
#'   `reducedDims(spe)`.
#' @param graph_name Name of the graph in `colPairs(spe)`. Only its nonzero
#'   topology is used; stored edge magnitudes are ignored. Cross-sample edges
#'   are rejected.
#' @param sample_col Name of the non-missing physical-sample identifier in
#'   `colData(spe)` used to isolate aggregation between samples.
#' @param output Name assigned to the aggregated spot-by-feature matrix in
#'   `reducedDims(spe)`. Provenance is stored in
#'   `metadata(spe)$CellCharterR$stages[[paste0("aggregation:", output)]]`.
#' @param n_layers Either one non-negative whole number, interpreted as all
#'   layers from zero through that value, or an explicit vector of non-negative
#'   whole-number layers that includes zero. The output concatenates the native
#'   representation at layer zero and means over spots at each requested exact
#'   graph hop; an empty hop contributes zeros.
#' @param filepath Optional HDF5 output path. The default `NULL` keeps the
#'   aggregated representation in memory and avoids an additional backing-file
#'   dependency. When supplied, `reducedDim(spe, output)` is HDF5-backed by the
#'   dataset named `output` in this file, so the file must remain available.
#' @param block_rows Positive number of spots read and written per block when
#'   `filepath` is supplied; it does not affect in-memory aggregation.
#' @param progress A single logical value indicating whether to display sample-
#'   or HDF5-block-level progress.
#' @return The updated `SpatialExperiment`, with one input-feature block per
#'   requested layer stored in `reducedDim(spe, output)`.
#' @examples
#' example_spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_spe <- cellcharter_aggregate(
#'   example_spe, output = "example_aggregation", n_layers = 1
#' )
#' dim(SingleCellExperiment::reducedDim(example_spe, "example_aggregation"))
#' @export
cellcharter_aggregate <- function(
    spe, representation = "X_trVAE", graph_name = "cellcharter_graph",
    sample_col = "sample", output = "X_cellcharter", n_layers = 3L,
    filepath = NULL, block_rows = 8192L, progress = interactive()) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  x <- .cc_spe_reduced_dim(spe, representation)
  graph <- .cc_spe_graph_matrix(spe, graph_name)
  sample <- as.character(.cc_spe_coldata(spe, sample_col, "sample column"))
  if (any(!nzchar(trimws(sample)))) {
    stop("Sample identifiers must not be blank.", call. = FALSE)
  }
  .cc_spe_assert_sample_local_graph(spe, sample, graph_name)
  rows <- split(seq_len(ncol(spe)), factor(sample, levels = unique(sample)))
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
    value <- matrix(NA_real_, nrow = ncol(spe), ncol = ncol(pieces[[1L]]))
    for (id in names(rows)) value[rows[[id]], ] <- as.matrix(pieces[[id]])
  } else {
    value <- .cc_aggregate_hdf5(
      x, adjacency, rows, filepath = filepath, name = output,
      n_layers = n_layers, block_rows = block_rows, progress = progress
    )
  }
  rownames(value) <- colnames(spe)
  colnames(value) <- paste0("feature_", seq_len(ncol(value)))
  SingleCellExperiment::reducedDim(spe, output) <- value
  .cc_spe_record_stage(
    spe, paste0("aggregation:", output),
    list(representation = representation, graph_name = graph_name,
         sample_col = sample_col, output = output, n_layers = n_layers,
         storage = if (is.null(filepath)) "memory" else "HDF5",
         filepath = if (is.null(filepath)) NULL else normalizePath(filepath))
  )
}

#' Run CellCharter AutoK and attach selected domains to a SpatialExperiment
#'
#' @section Backends:
#' The compiled CPU backend (`backend = "cpp"` in `...`) is the supported,
#' scalable implementation. The R Torch backend (`backend = "torch"`) remains
#' experimental. Both use full-covariance Gaussian mixtures; consult
#' [cellcharter_autok_fit()] for their execution and checkpoint controls.
#'
#' @param spe An aggregated `SpatialExperiment`.
#' @param representation Name of the numeric spot-by-feature matrix in
#'   `reducedDims(spe)`. The matrix is a first-class input and may have been
#'   created by CellCharterR, imported from AnnData, or attached independently;
#'   historical stage metadata is not required.
#' @param candidate_k Two selectable K bounds or an expanded consecutive list,
#'   as accepted by [cellcharter_autok_fit()]. AutoK additionally fits the two
#'   neighboring sentinel K values needed for stability calculation.
#' @param checkpoint_dir Non-empty path used by AutoK. Its exact write behavior
#'   is controlled by `checkpoint_mode` passed through `...`; mode `"none"`
#'   performs no writes even though this path is still required by AutoK.
#' @param domain_col Name of the factor written to `colData(spe)`. Labels are
#'   the selected model's zero-based component identifiers. The complete AutoK
#'   object is stored in `metadata(spe)$CellCharterR$autok[[domain_col]]`, and
#'   stage provenance in `metadata(spe)$CellCharterR$stages[[paste0(
#'   "clustering:", domain_col)]]`.
#' @param ... Named arguments passed to [cellcharter_autok_fit()], including
#'   `backend`, `max_runs`, `convergence_tol`, `seed`, `fit_args`, `workers`,
#'   `threads_per_fit`, `checkpoint_mode`, `resume`, and `progress`. Do not pass
#'   `x`, `n_clusters`, `checkpoint_dir`, or `data_id`; these are supplied from
#'   the current SPE representation.
#' @return The updated `SpatialExperiment` with selected labels in
#'   `colData(spe)[[domain_col]]` and the AutoK fit in metadata.
#' @examples
#' example_spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' demo_fit <- function(x, n_clusters, seed, ...) {
#'   list(k = n_clusters, seed = seed)
#' }
#' demo_predict <- function(model, newdata, ...) {
#'   (seq_len(nrow(newdata)) - 1L) %% model$k
#' }
#' checkpoint <- tempfile("CellCharterR-example-")
#' example_spe <- cellcharter_cluster(
#'   example_spe, representation = "X_cellcharter_tiny",
#'   candidate_k = c(2L, 3L), checkpoint_dir = checkpoint,
#'   domain_col = "example_domain", max_runs = 2L, convergence_tol = 0,
#'   fit_fun = demo_fit, predict_fun = demo_predict,
#'   nll_fun = function(model) as.numeric(model$k)
#' )
#' table(example_spe$example_domain)
#' unlink(checkpoint, recursive = TRUE)
#' @export
cellcharter_cluster <- function(
    spe, representation = "X_cellcharter", candidate_k,
    checkpoint_dir, domain_col = "cellcharter_domain", ...) {
  x <- .cc_spe_reduced_dim(spe, representation)
  identity <- .cc_validate_spe_identity(spe, reduced_dims = representation)
  # AutoK fingerprints the complete current representation and binds resume
  # compatibility to that content and the numerical execution contract. Never
  # mix mutable producer-stage metadata or storage paths into data identity.
  fit <- cellcharter_autok_fit(
    x = x, n_clusters = candidate_k, checkpoint_dir = checkpoint_dir,
    data_id = NULL,
    ...
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
  if (length(labels) != ncol(spe) || anyNA(labels)) {
    stop("Selected AutoK labels do not align one-to-one with SPE spots.", call. = FALSE)
  }
  spe[[domain_col]] <- factor(labels, levels = sort(unique(labels)))
  metadata <- .cc_spe_metadata(spe)
  metadata$autok[[domain_col]] <- fit
  S4Vectors::metadata(spe)$CellCharterR <- metadata
  .cc_spe_record_stage(
    spe, paste0("clustering:", domain_col),
    list(representation = representation, candidate_k = candidate_k,
         checkpoint_dir = if (identical(fit$checkpoint_mode, "none")) {
           NULL
         } else {
           normalizePath(checkpoint_dir)
         },
         domain_col = domain_col, selected_k = selected_k,
         representation_sha256 = fit$representation_sha256,
         autok_input_data_sha256 = fit$representation_sha256,
         observation_order_sha256 = identity$cell_order_sha256,
         selected_artifact = if (identical(fit$checkpoint_mode, "none")) {
           NULL
         } else {
           normalizePath(selected$path[[1L]])
         })
  )
}

.cc_spe_store_result <- function(spe, result_name, result, parameters,
                                 input_specs = NULL) {
  metadata <- .cc_spe_metadata(spe)
  identity <- .cc_validate_spe_identity(spe)
  metadata$results[[result_name]] <- list(
    value = result, cell_order_sha256 = identity$cell_order_sha256,
    n_cells = ncol(spe), parameters = parameters,
    parameters_sha256 = .cc_result_parameters_sha256(parameters),
    input_fingerprints = .cc_result_inputs(spe, input_specs)
  )
  S4Vectors::metadata(spe)$CellCharterR <- metadata
  spe
}

#' Attach CellCharter group enrichment to a SpatialExperiment
#'
#' The stored result fingerprints its group and domain annotations and analysis
#' parameters. Experiment-aware plotting rejects the result if those inputs
#' subsequently change.
#'
#' @param spe A `SpatialExperiment` containing one group and one domain label
#'   per spot.
#' @param group_col Name of the non-missing cell-type, program, or other group
#'   annotation in `colData(spe)`.
#' @param domain_col Name of the non-missing domain labels in `colData(spe)`,
#'   the name of a one-column matrix in `reducedDims(spe)`, or an unambiguous
#'   column name within a reduced dimension.
#' @param result_name Name used to store the freshness record at
#'   `metadata(spe)$CellCharterR$results[[result_name]]`. The scientific result
#'   is its `$value`, a list with `observed`, `expected`, `enrichment`, and
#'   `params`; reusing a name replaces the previous record.
#' @param log2_fold_change A single logical value. `TRUE` returns
#'   `log2(observed / expected)`; `FALSE` returns `observed / expected`.
#'   Observed values are within-domain group proportions and expected values
#'   are cohort-wide group proportions.
#' @param progress A single logical value indicating whether to display
#'   chunk-level counting progress.
#' @param chunk_size Maximum number of observations processed per counting
#'   chunk. This bounds temporary count-index memory without changing results.
#' @return The updated `SpatialExperiment`; no assay or reduced dimension is
#'   modified.
#' @examples
#' example_spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_spe <- cellcharter_enrichment(
#'   example_spe, group_col = "cell_type_reduced", domain_col = "tiny_domain"
#' )
#' names(S4Vectors::metadata(example_spe)$CellCharterR$results$enrichment$value)
#' @export
cellcharter_enrichment <- function(
    spe, group_col, domain_col = "cellcharter_domain",
    result_name = "enrichment", log2_fold_change = TRUE,
    progress = interactive(), chunk_size = 1000000L) {
  group <- .cc_spe_coldata(spe, group_col, "group column")
  domain <- .cc_spe_domain(spe, domain_col)
  result <- .cc_enrichment(
    group, domain, log2_fold_change,
    progress = progress, chunk_size = chunk_size
  )
  .cc_spe_store_result(
    spe, result_name, result,
    list(group_col = group_col, domain_col = domain_col,
         log2_fold_change = log2_fold_change, progress = progress,
         chunk_size = as.integer(chunk_size)),
    list(
      group = list(kind = "colData", name = group_col),
      domain = list(kind = "domain", name = domain_col)
    )
  )
}

#' Attach CellCharter neighborhood enrichment to a SpatialExperiment
#'
#' The stored result fingerprints the weighted graph, domain annotations, and
#' analysis parameters so stale results can be detected without materializing
#' the complete graph.
#' @inheritParams cellcharter_enrichment
#' @param graph_name Name of the symmetric graph in `colPairs(spe)`. Nonzero
#'   edge values define adjacency; direction is reported as source domain by
#'   target domain.
#' @param log2_fold_change A single logical value. `TRUE` returns
#'   `log2(observed / expected)`; `FALSE` reproduces the CellCharter default
#'   `observed - expected` difference.
#' @param only_inter A single logical value. When `TRUE`, same-domain edges are
#'   removed before calculating proportions and the diagonal is returned as
#'   `NA`.
#' @return The updated `SpatialExperiment`. The record at
#'   `metadata(spe)$CellCharterR$results[[result_name]]$value` contains directed
#'   `observed`, analytical `expected`, `enrichment`, and `params` matrices;
#'   no graph or annotation is modified.
#' @examples
#' example_spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_spe <- cellcharter_nhood_enrichment(
#'   example_spe, domain_col = "tiny_domain"
#' )
#' S4Vectors::metadata(example_spe)$CellCharterR$results$nhood_enrichment$value$enrichment
#' @export
cellcharter_nhood_enrichment <- function(
    spe, graph_name = "cellcharter_graph", domain_col = "cellcharter_domain",
    result_name = "nhood_enrichment", log2_fold_change = FALSE,
    only_inter = TRUE) {
  graph <- .cc_spe_graph_matrix(spe, graph_name)
  domain <- .cc_spe_domain(spe, domain_col)
  result <- .cc_nhood_enrichment(graph, domain, log2_fold_change, only_inter)
  .cc_spe_store_result(
    spe, result_name, result,
    list(graph_name = graph_name, domain_col = domain_col,
         log2_fold_change = log2_fold_change, only_inter = only_inter),
    list(
      graph = list(kind = "weighted_graph", name = graph_name),
      domain = list(kind = "domain", name = domain_col)
    )
  )
}

#' Attach differential neighborhood enrichment to a SpatialExperiment
#'
#' The stored result fingerprints the weighted graph, domains, conditions,
#' sample/library assignments, and analysis parameters. Downstream plotting
#' rejects a result after any of those scientific inputs changes.
#' @inheritParams cellcharter_nhood_enrichment
#' @param condition_col Name of the non-missing experimental-condition
#'   annotation in `colData(spe)`.
#' @param sample_col Name of the physical sample or library in `colData(spe)`.
#'   Each library must belong to one condition and graph edges may not cross
#'   libraries; permutations operate on these library labels, not spots.
#' @param condition_groups Two distinct observed conditions, or `NULL` to use
#'   all present condition levels when there are exactly two. Their factor-level
#'   order, not the order supplied here, defines the stored numerator-minus-
#'   denominator effect.
#' @param n_perms Positive integer number of library-level permutations.
#' @param seed Finite whole-number seed used for permutation sampling.
#' @param permutation_mode Either `"cellcharter_v0.2.0"` for the pinned
#'   with-replacement CellCharter procedure or `"shuffle"` for a conventional
#'   size-preserving without-replacement assignment.
#' @param progress A single logical value indicating whether to report setup and
#'   display permutation progress.
#' @return The updated `SpatialExperiment`. The result at
#'   `metadata(spe)$CellCharterR$results[[result_name]]$value` contains each
#'   condition's `observed`, `expected`, and `enrichment`, plus `effect`,
#'   `permutation_counts`, directional empirical `pvalue`, and `params`.
#' @examples
#' example_spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' example_spe <- cellcharter_diff_nhood_enrichment(
#'   example_spe, domain_col = "tiny_domain",
#'   condition_groups = c("BALBc", "MRL"), n_perms = 5L, seed = 1L
#' )
#' S4Vectors::metadata(example_spe)$CellCharterR$results$diff_nhood_enrichment$value$pvalue
#' @export
cellcharter_diff_nhood_enrichment <- function(
    spe, graph_name = "cellcharter_graph", domain_col = "cellcharter_domain",
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
  graph <- .cc_spe_graph_matrix(spe, graph_name)
  domain <- .cc_spe_domain(spe, domain_col)
  condition <- .cc_spe_coldata(spe, condition_col, "condition column")
  sample <- .cc_spe_coldata(spe, sample_col, "sample column")
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
  .cc_spe_store_result(
    spe, result_name, result,
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
#' This convenience function composes the exported SpatialExperiment stages. Expensive
#' full datasets should use the staged functions directly so each stage can be
#' checkpointed; this orchestrator is intended for compact examples and tests.
#'
#' @param spe A `SpatialExperiment` containing spatial coordinates, an input
#'   representation, sample identifiers, conditions, and group annotations.
#' @param group_col Name of the group annotation in `colData(spe)` used for
#'   observed-versus-expected enrichment.
#' @param candidate_k Candidate K bounds or consecutive sequence forwarded to
#'   [cellcharter_cluster()].
#' @param checkpoint_dir Non-empty AutoK working directory forwarded to
#'   [cellcharter_cluster()].
#' @param spatial Coordinate source forwarded to [cellcharter_graph()].
#' @param representation Name of the initial representation in
#'   `reducedDims(spe)` forwarded to [cellcharter_aggregate()].
#' @param sample_col Name of the physical sample/library in `colData(spe)` used
#'   for graph isolation, aggregation isolation, and differential permutations.
#' @param condition_col Name of the experimental condition in `colData(spe)`.
#' @param graph_name Name written to `colPairs(spe)` by the graph stage and read
#'   by aggregation and both neighborhood analyses.
#' @param aggregated Name written to `reducedDims(spe)` by aggregation and used
#'   as the clustering representation.
#' @param domain_col Name written to `colData(spe)` for selected zero-based
#'   cluster labels and used by all downstream analyses.
#' @param n_layers Whole-number layer specification passed to
#'   [cellcharter_aggregate()]; either a maximum hop or an explicit vector that
#'   includes layer zero.
#' @param condition_groups Two conditions included in the differential
#'   neighborhood contrast, or `NULL` when exactly two condition levels are
#'   present. Factor-level order determines effect direction.
#' @param n_perms Positive number of library-level differential permutations.
#' @param seed Base seed supplied to AutoK and the differential permutation
#'   analysis.
#' @param autok_args Named list of additional [cellcharter_autok_fit()]
#'   arguments. Names already supplied by this orchestrator (`spe`,
#'   `representation`, `candidate_k`, `checkpoint_dir`, `domain_col`, and
#'   `seed`) must not be repeated.
#' @return The updated `SpatialExperiment`. Graph, aggregated representation,
#'   domains, AutoK metadata, and enrichment results are stored in the exact
#'   compartments documented by the individual stage functions.
#' @examples
#' codex_mouse_spleen_tiny <- cellcharter_mouse_spleen_tiny()
#' demo_fit <- function(x, n_clusters, seed, ...) {
#'   list(k = n_clusters, seed = seed)
#' }
#' demo_predict <- function(model, newdata, ...) {
#'   (seq_len(nrow(newdata)) - 1L) %% model$k
#' }
#' checkpoint <- tempfile("CellCharterR-workflow-")
#' example_spe <- cellcharter_workflow(
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
#' example_spe
#' unlink(checkpoint, recursive = TRUE)
#' @export
cellcharter_workflow <- function(
    spe, group_col, candidate_k, checkpoint_dir,
    spatial = "spatial", representation = "X_trVAE", sample_col = "sample",
    condition_col = "condition", graph_name = "cellcharter_graph",
    aggregated = "X_cellcharter", domain_col = "cellcharter_domain",
    n_layers = 3L, condition_groups = NULL, n_perms = 100L, seed = 1L,
    autok_args = list()) {
  spe <- cellcharter_graph(spe, spatial, sample_col, graph_name)
  spe <- cellcharter_aggregate(
    spe, representation, graph_name, sample_col, aggregated, n_layers
  )
  cluster_args <- c(
    list(spe = spe, representation = aggregated, candidate_k = candidate_k,
         checkpoint_dir = checkpoint_dir, domain_col = domain_col, seed = seed),
    autok_args
  )
  spe <- do.call(cellcharter_cluster, cluster_args)
  spe <- cellcharter_enrichment(spe, group_col, domain_col)
  spe <- cellcharter_nhood_enrichment(spe, graph_name, domain_col)
  cellcharter_diff_nhood_enrichment(
    spe, graph_name, domain_col, condition_col, sample_col,
    condition_groups = condition_groups, n_perms = n_perms, seed = seed
  )
}
