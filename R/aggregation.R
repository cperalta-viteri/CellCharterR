#' Remove the longest spatial graph edges
#'
#' @param adjacency symmetric sparse connectivity matrix.
#' @param coordinates two-column coordinate matrix in the same row order.
#' @param percentile edges longer than this global positive-distance percentile
#'   are removed. CellCharter tutorials use 99. When `adjacency` is a combined
#'   block-diagonal graph, one cutoff is therefore calculated across all blocks.
#' @param cutoff Optional precomputed distance cutoff. Use this when graphs are
#'   stored per sample but the reference cutoff was calculated globally across
#'   all samples. When supplied, `percentile` is recorded but not recomputed.
#' @return A sparse adjacency matrix with only retained edges. The cutoff is
#'   stored in the `distance_cutoff` attribute.
#' @export
remove_long_links <- function(adjacency, coordinates, percentile = 99, cutoff = NULL) {
  topology <- .cc_adjacency_topology(adjacency)
  coordinates <- .cc_coordinates(coordinates, nrow(topology))
  if (!is.numeric(percentile) || length(percentile) != 1L ||
      is.na(percentile) || !is.finite(percentile) ||
      percentile < 0 || percentile > 100) {
    stop("`percentile` must be one finite number between 0 and 100.", call. = FALSE)
  }

  if (!is.null(cutoff) &&
      (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) ||
       !is.finite(cutoff) || cutoff < 0)) {
    stop("`cutoff` must be NULL or one finite non-negative number.", call. = FALSE)
  }

  ij <- Matrix::summary(topology)
  if (nrow(ij) == 0L) {
    answer <- .cc_binary_graph(topology)
    attr(answer, "distance_cutoff") <- if (is.null(cutoff)) NA_real_ else cutoff
    return(answer)
  }

  # v0.2.0 takes the percentile over all positive entries in the sparse
  # distance matrix. For an undirected graph this includes both directions.
  distance <- sqrt(rowSums((coordinates[ij$i, , drop = FALSE] - coordinates[ij$j, , drop = FALSE])^2))
  positive_distance <- distance[distance > 0]
  if (length(positive_distance) == 0L) {
    answer <- .cc_binary_graph(topology)
    attr(answer, "distance_cutoff") <- if (is.null(cutoff)) NA_real_ else cutoff
    return(answer)
  }

  if (is.null(cutoff)) {
    cutoff <- stats::quantile(positive_distance, percentile / 100, names = FALSE, type = 7)
  }
  # Zero-distance links are not part of the percentile in the Python reference,
  # but are retained because they are not longer than the positive cutoff.
  keep <- distance <= cutoff
  answer <- Matrix::sparseMatrix(
    i = ij$i[keep], j = ij$j[keep], x = 1,
    dims = dim(topology), dimnames = dimnames(topology), symmetric = FALSE
  )
  attr(answer, "distance_cutoff") <- as.numeric(cutoff)
  answer
}

#' Compute the pinned global long-link cutoff from sample-local graphs
#'
#' CellCharter 0.2.0 builds Delaunay graphs independently by physical sample,
#' assembles a block-diagonal distance matrix, and then takes one percentile of
#' all its positive entries. This helper reproduces the same cutoff without
#' materialising the cohort-wide block-diagonal graph.
#'
#' @param adjacency_by_sample Named list of sample-local sparse graphs.
#' @param coordinates_by_sample Named list of aligned coordinate matrices.
#' @param percentile Distance percentile, normally 99.
#' @return The scalar global distance cutoff.
#' @export
cellcharter_global_distance_cutoff <- function(
    adjacency_by_sample, coordinates_by_sample, percentile = 99) {
  if (!is.list(adjacency_by_sample) || !is.list(coordinates_by_sample) ||
      !length(adjacency_by_sample) ||
      !identical(names(adjacency_by_sample), names(coordinates_by_sample))) {
    stop("Graph and coordinate lists must be non-empty and identically named.", call. = FALSE)
  }
  if (!is.numeric(percentile) || length(percentile) != 1L || is.na(percentile) ||
      !is.finite(percentile) || percentile < 0 || percentile > 100) {
    stop("`percentile` must be one finite number between 0 and 100.", call. = FALSE)
  }
  distances <- lapply(names(adjacency_by_sample), function(sample_id) {
    graph <- .cc_adjacency_topology(adjacency_by_sample[[sample_id]])
    coordinates <- .cc_coordinates(coordinates_by_sample[[sample_id]], nrow(graph))
    ij <- Matrix::summary(graph)
    if (!nrow(ij)) return(numeric())
    value <- sqrt(rowSums((coordinates[ij$i, , drop = FALSE] -
                           coordinates[ij$j, , drop = FALSE])^2))
    value[value > 0]
  })
  distances <- unlist(distances, use.names = FALSE)
  if (!length(distances)) return(NA_real_)
  as.numeric(stats::quantile(distances, percentile / 100, names = FALSE, type = 7))
}

.cc_coordinates <- function(coordinates, expected_rows = NULL) {
  if (is.null(dim(coordinates)) || length(dim(coordinates)) != 2L) {
    stop("`coordinates` must be a two-dimensional numeric object.", call. = FALSE)
  }
  coordinates <- as.matrix(coordinates)
  if (!is.numeric(coordinates) || ncol(coordinates) < 2L) {
    stop("`coordinates` must contain numeric x/y values in its first two columns.", call. = FALSE)
  }
  coordinates <- coordinates[, 1:2, drop = FALSE]
  if (!is.null(expected_rows) && nrow(coordinates) != expected_rows) {
    stop("`coordinates` and `adjacency` must contain the same number of cells.", call. = FALSE)
  }
  if (any(!is.finite(coordinates))) {
    stop("`coordinates` must contain only finite x/y values.", call. = FALSE)
  }
  coordinates
}

.cc_adjacency_topology <- function(adjacency, expected_rows = NULL) {
  if (!inherits(adjacency, "sparseMatrix")) {
    stop("`adjacency` must be a sparse Matrix object.", call. = FALSE)
  }
  if (length(dim(adjacency)) != 2L || nrow(adjacency) != ncol(adjacency)) {
    stop("`adjacency` must be square.", call. = FALSE)
  }
  if (!is.null(expected_rows) && nrow(adjacency) != expected_rows) {
    stop("`adjacency` must have one row and column per observation in `x`.", call. = FALSE)
  }

  # Coerce symmetric/diagonal sparse classes to a general sparse representation
  # before inspecting stored entries. Pattern matrices remain sparse.
  adjacency <- methods::as(adjacency, "generalMatrix")
  entries <- Matrix::summary(adjacency)
  if (nrow(entries) > 0L &&
      (anyNA(entries$x) || any(!is.finite(entries$x)) || any(entries$x < 0))) {
    stop("`adjacency` must contain finite, non-negative values.", call. = FALSE)
  }
  topology <- Matrix::drop0(adjacency != 0)
  Matrix::diag(topology) <- FALSE
  topology <- Matrix::drop0(topology)
  if (!isTRUE(Matrix::isSymmetric(topology))) {
    stop("`adjacency` must have symmetric connectivity.", call. = FALSE)
  }
  topology
}

.cc_binary_graph <- function(topology) {
  ij <- Matrix::summary(topology)
  Matrix::sparseMatrix(
    i = ij$i, j = ij$j, x = rep.int(1, nrow(ij)),
    dims = dim(topology), dimnames = dimnames(topology), symmetric = FALSE
  )
}

.cc_layers <- function(n_layers) {
  if (!is.numeric(n_layers) || length(n_layers) == 0L ||
      anyNA(n_layers) || any(!is.finite(n_layers)) ||
      any(n_layers < 0) || any(n_layers != floor(n_layers))) {
    stop("`n_layers` must contain finite, non-negative whole numbers.", call. = FALSE)
  }
  layers <- if (length(n_layers) == 1L) {
    seq.int(0L, as.integer(n_layers))
  } else {
    sort(unique(as.integer(n_layers)))
  }
  if (layers[[1L]] != 0L) {
    stop("`n_layers` must include 0 so the native representation is retained.", call. = FALSE)
  }
  layers
}

.cc_row_normalise <- function(a) {
  degree <- Matrix::rowSums(a)
  inverse_degree <- ifelse(degree == 0, 0, 1 / degree)
  Matrix::Diagonal(x = inverse_degree) %*% a
}

.cc_hop <- function(current, adjacency, visited) {
  # Binarise path counts before the sparse `>` operation. This is the intended
  # logical-not operation in CellCharter v0.2.0: a node reached by several paths
  # must not overcome a TRUE entry in `visited` and re-enter a later ring.
  candidate <- Matrix::drop0((current %*% adjacency) > 0)
  candidate <- Matrix::drop0(candidate > visited)
  list(
    current = candidate,
    visited = Matrix::drop0((visited + candidate) > 0)
  )
}

.cc_contiguous_runs <- function(rows) {
  # HDF5Array writes rectangular blocks. Samples are usually contiguous after
  # AnnData concatenation, but this also handles interleaved row indices safely.
  rows <- sort(as.integer(rows))
  split(rows, cumsum(c(TRUE, diff(rows) != 1L)))
}

#' CellCharter-style graph-neighbourhood aggregation
#'
#' Aggregates the original representation at 0 hop and mean neighbour features
#' at each requested *exclusive* graph hop. A graph must be supplied separately
#' for each slide/sample so no artificial cross-section edges can appear.
#'
#' @param x observation-by-feature dense or delayed matrix.
#' @param adjacency sample-local sparse adjacency matrix.
#' @param n_layers largest neighbourhood distance, or explicit layer vector.
#' @return A (possibly delayed) observation-by-feature matrix with concatenated
#'   0, 1, ..., `n_layers` features.
#' @export
cellcharter_aggregate <- function(x, adjacency, n_layers = 3L) {
  if (is.null(dim(x)) || length(dim(x)) != 2L || nrow(x) == 0L || ncol(x) == 0L) {
    stop("`x` must be a non-empty observation-by-feature matrix-like object.", call. = FALSE)
  }
  if (is.matrix(x) && (!is.numeric(x) || any(!is.finite(x)))) {
    stop("A dense `x` must contain only finite numeric values.", call. = FALSE)
  }
  layers <- .cc_layers(n_layers)
  adjacency <- .cc_adjacency_topology(adjacency, nrow(x))
  current <- adjacency
  visited <- Matrix::drop0((current + Matrix::Diagonal(nrow(adjacency))) > 0)
  delayed_input <- inherits(x, "DelayedArray")
  output <- list(x)
  for (hop in seq_len(max(layers))) {
    if (hop > 1L) {
      state <- .cc_hop(current, adjacency, visited)
      current <- state$current; visited <- state$visited
    }
    if (hop %in% layers) {
      aggregated <- .cc_row_normalise(current) %*% x
      # Keep every block in the same delayed family. Mixing a DelayedMatrix
      # native block with a Matrix hop block makes cbind dispatch to bindCOLS
      # without a compatible method on real HDF5-backed tutorial input.
      if (delayed_input) aggregated <- DelayedArray::DelayedArray(aggregated)
      output[[length(output) + 1L]] <- aggregated
    }
  }
  if (delayed_input) {
    do.call(DelayedArray::acbind, output)
  } else {
    do.call(cbind, output)
  }
}

#' Write CellCharter aggregate features to HDF5 without cohort materialisation
#'
#' Each input sample is aggregated alone and assigned into an HDF5-backed
#' matrix. The delayed target is never realised in RAM. Keep writing serial;
#' independent sample computation can be prepared in parallel, but concurrent
#' writers to one HDF5 file are unsafe.
#'
#' @param x observation-by-feature delayed matrix.
#' @param adjacency_by_sample named list of sparse sample-local graphs.
#' @param rows_by_sample named list of matching cohort row indices.
#' @param filepath output HDF5 file.
#' @param name HDF5 dataset name.
#' @param n_layers number of neighbourhood layers.
#' @return An HDF5-backed `DelayedMatrix`.
#' @export
cellcharter_aggregate_hdf5 <- function(x, adjacency_by_sample, rows_by_sample,
                                       filepath, name = "X_cellcharter", n_layers = 3L) {
  if (!requireNamespace("HDF5Array", quietly = TRUE)) stop("Install HDF5Array for disk-backed output.")
  layers <- .cc_layers(n_layers)
  sample_names <- names(rows_by_sample)
  if (!is.list(adjacency_by_sample) || !is.list(rows_by_sample) ||
      length(rows_by_sample) == 0L || is.null(sample_names) ||
      any(!nzchar(sample_names)) || anyDuplicated(sample_names) ||
      !identical(names(adjacency_by_sample), sample_names)) {
    stop("Graph and row lists must be non-empty, uniquely named, and have identical names.", call. = FALSE)
  }
  checked_rows <- lapply(rows_by_sample, function(rows) {
    if (!is.numeric(rows) || length(rows) == 0L || anyNA(rows) ||
        any(!is.finite(rows)) || any(rows != floor(rows)) ||
        any(rows < 1L) || any(rows > nrow(x)) || anyDuplicated(rows)) {
      stop("Each `rows_by_sample` entry must contain unique, valid row indices.", call. = FALSE)
    }
    as.integer(rows)
  })
  all_rows <- unlist(checked_rows, use.names = FALSE)
  if (length(all_rows) != nrow(x) || anyDuplicated(all_rows) ||
      !identical(sort(all_rows), seq_len(nrow(x)))) {
    stop("`rows_by_sample` must partition all rows of `x` exactly once.", call. = FALSE)
  }
  rows_by_sample <- checked_rows
  n_features <- ncol(x) * length(layers)
  # A realization sink creates the on-disk array first. Unlike
  # writeHDF5Array(matrix(0, ...)), it never allocates a cohort-sized zero
  # matrix in R memory.
  target <- HDF5Array::HDF5RealizationSink(
    dim = c(nrow(x), n_features), filepath = filepath, name = name,
    type = "double"
  )
  for (sample_id in names(rows_by_sample)) {
    rows <- rows_by_sample[[sample_id]]
    # The selected sample is the only data block allowed into memory.
    aggregate <- cellcharter_aggregate(x[rows, , drop = FALSE], adjacency_by_sample[[sample_id]], layers)
    # Only contiguous runs can be written as an HDF5 rectangular block. This
    # preserves original cohort order while avoiding any global realization.
    for (run in .cc_contiguous_runs(rows)) {
      local_rows <- match(run, rows)
      viewport <- S4Arrays::ArrayViewport(
        refdim = c(nrow(x), n_features),
        # One IRanges entry describes the selected extent in each dimension.
        ranges = IRanges::IRanges(
          start = c(run[[1]], 1L), width = c(length(run), n_features)
        )
      )
      DelayedArray::write_block(target, viewport, as.matrix(aggregate[local_rows, , drop = FALSE]))
    }
  }
  HDF5Array::HDF5Array(filepath, name)
}
#' Build a Delaunay graph for one physical sample
#'
#' This intentionally accepts coordinates from only one slide/sample. Call it
#' separately per sample and never run Delaunay on a concatenated cohort: doing
#' so would create biologically meaningless edges across tissue sections. The
#' Qhull options mirror SciPy's two-dimensional defaults. For coincident points,
#' Qhull retains one representative and the omitted duplicates remain isolated.
#'
#' @param coordinates numeric matrix with x and y in columns one and two.
#' @return Symmetric binary sparse adjacency matrix.
#' @export
delaunay_graph <- function(coordinates) {
  coordinates <- .cc_coordinates(coordinates)
  n <- nrow(coordinates)
  graph_dimnames <- list(rownames(coordinates), rownames(coordinates))
  if (n == 0L) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(0L, 0L), dimnames = graph_dimnames))
  }
  if (n == 1L) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(1L, 1L), dimnames = graph_dimnames))
  }
  if (n == 2L) {
    return(Matrix::sparseMatrix(
      i = c(1L, 2L), j = c(2L, 1L), x = 1,
      dims = c(2L, 2L), dimnames = graph_dimnames
    ))
  }
  centered <- sweep(coordinates, 2L, colMeans(coordinates), FUN = "-")
  if (qr(centered)$rank < 2L) {
    stop("Delaunay triangulation requires at least three non-collinear coordinate pairs.", call. = FALSE)
  }
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("Install `geometry` for Qhull-compatible Delaunay graph construction.")
  }
  # Squidpy 1.3.0 calls scipy.spatial.Delaunay, which uses Qhull with these
  # options for two-dimensional inputs. Qhull retains one representative of
  # coincident coordinates and leaves the other observations isolated.
  triangles <- tryCatch(
    suppressWarnings(geometry::delaunayn(
      coordinates, options = "Qbb Qc Qz Q12 Qt"
    )),
    error = function(e) {
      stop("Delaunay triangulation failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  used <- sort(unique(as.integer(triangles)))
  omitted <- setdiff(seq_len(n), used)
  if (length(omitted)) {
    keys <- paste(coordinates[, 1L], coordinates[, 2L], sep = "\r")
    if (any(!keys[omitted] %in% keys[used])) {
      stop("Qhull omitted a non-duplicate coordinate; graph construction is not reliable for this geometry.", call. = FALSE)
    }
  }
  edges <- rbind(
    triangles[, c(1L, 2L), drop = FALSE],
    triangles[, c(1L, 3L), drop = FALSE],
    triangles[, c(2L, 3L), drop = FALSE]
  )
  graph <- Matrix::sparseMatrix(
    i = c(edges[, 1L], edges[, 2L]),
    j = c(edges[, 2L], edges[, 1L]), x = 1,
    dims = c(n, n), dimnames = graph_dimnames
  )
  Matrix::drop0(graph > 0)
}
