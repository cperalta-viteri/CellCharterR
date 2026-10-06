#' Remove the longest spatial graph edges
#'
#' @param adjacency symmetric sparse connectivity matrix.
#' @param coordinates two-column coordinate matrix in the same row order.
#' @param percentile edges longer than this within-sample distance percentile are
#'   removed. CellCharter tutorials use 99.
#' @return A sparse adjacency matrix with only retained edges.
#' @export
remove_long_links <- function(adjacency, coordinates, percentile = 99) {
  stopifnot(inherits(adjacency, "sparseMatrix"), nrow(adjacency) == nrow(coordinates))
  ij <- Matrix::summary(adjacency)
  # Work on unique undirected edges: otherwise every distance is counted twice.
  ij <- ij[ij$i < ij$j, , drop = FALSE]
  distance <- sqrt(rowSums((coordinates[ij$i, , drop = FALSE] - coordinates[ij$j, , drop = FALSE])^2))
  cutoff <- stats::quantile(distance, percentile / 100, names = FALSE, type = 7)
  keep <- distance <= cutoff
  Matrix::sparseMatrix(
    i = c(ij$i[keep], ij$j[keep]), j = c(ij$j[keep], ij$i[keep]), x = 1,
    dims = dim(adjacency), symmetric = FALSE
  )
}

.cc_row_normalise <- function(a) {
  degree <- Matrix::rowSums(a)
  inverse_degree <- ifelse(degree == 0, 0, 1 / degree)
  Matrix::Diagonal(x = inverse_degree) %*% a
}

.cc_hop <- function(current, adjacency, visited) {
  # `current` contains paths exactly h hops long. Removing previously visited
  # nodes makes each layer exclusive, exactly as CellCharter's sparse algebra.
  candidate <- current %*% adjacency
  candidate <- Matrix::drop0(candidate > visited)
  list(current = candidate, visited = Matrix::drop0(visited + candidate))
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
  stopifnot(nrow(x) == nrow(adjacency), inherits(adjacency, "sparseMatrix"))
  layers <- if (length(n_layers) == 1L) seq.int(0L, as.integer(n_layers)) else sort(unique(as.integer(n_layers)))
  if (min(layers) != 0L) stop("`n_layers` must include 0 so the native representation is retained.")
  adjacency <- Matrix::drop0(Matrix::tril(adjacency, k = -1) + Matrix::triu(adjacency, k = 1))
  current <- adjacency > 0
  visited <- current + Matrix::Diagonal(nrow(adjacency))
  output <- list(x)
  for (hop in seq_len(max(layers))) {
    if (hop > 1L) {
      state <- .cc_hop(current, adjacency, visited)
      current <- state$current; visited <- state$visited
    }
    if (hop %in% layers) output[[length(output) + 1L]] <- .cc_row_normalise(current) %*% x
  }
  do.call(cbind, output)
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
  stopifnot(identical(names(adjacency_by_sample), names(rows_by_sample)))
  n_features <- ncol(x) * (as.integer(n_layers) + 1L)
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
    aggregate <- cellcharter_aggregate(x[rows, , drop = FALSE], adjacency_by_sample[[sample_id]], n_layers)
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
#' so would create biologically meaningless edges across tissue sections.
#'
#' @param coordinates numeric matrix with x and y in columns one and two.
#' @return Symmetric binary sparse adjacency matrix.
#' @export
delaunay_graph <- function(coordinates) {
  if (!requireNamespace("deldir", quietly = TRUE)) stop("Install `deldir` for Delaunay graph construction.")
  coordinates <- as.matrix(coordinates)
  if (ncol(coordinates) < 2L || any(!is.finite(coordinates[, 1:2])))
    stop("`coordinates` must contain finite x/y values in its first two columns.")
  # `delsgs` records endpoint indices, preserving the submitted row order.
  edges <- deldir::deldir(coordinates[, 1], coordinates[, 2])$delsgs
  Matrix::sparseMatrix(
    i = c(edges$ind1, edges$ind2), j = c(edges$ind2, edges$ind1), x = 1,
    dims = c(nrow(coordinates), nrow(coordinates))
  )
}
