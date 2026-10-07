test_that("long-link trimming uses one global cutoff for block-diagonal graphs", {
  first <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3), j = c(2, 1, 3, 2), x = 1,
    dims = c(3, 3)
  )
  second <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3), j = c(2, 1, 3, 2), x = 1,
    dims = c(3, 3)
  )
  graph <- Matrix::bdiag(first, second)

  # The samples deliberately reuse a coordinate origin. Only the supplied
  # block-diagonal connectivity, never coordinate proximity, defines edges.
  coordinates <- rbind(
    c(0, 0), c(1, 0), c(3, 0),
    c(0, 0), c(10, 0), c(110, 0)
  )
  trimmed <- remove_long_links(graph, coordinates, percentile = 75)

  # Directed positive distances are 1,1,2,2,10,10,100,100. The global 75th
  # percentile is 32.5, so only the length-100 edge is removed. A separate
  # per-sample cutoff would incorrectly remove the length-2 edge as well.
  expect_equal(as.matrix(trimmed)[1:3, 1:3], as.matrix(first))
  expect_equal(as.matrix(trimmed)[4:6, 4:6],
               matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, byrow = TRUE))
  expect_equal(sum(trimmed[1:3, 4:6]), 0)
  expect_equal(sum(trimmed[4:6, 1:3]), 0)
})

test_that("long-link trimming preserves isolates and handles empty graphs", {
  graph <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(3, 3))
  coordinates <- cbind(1:3, 0)
  trimmed <- remove_long_links(graph, coordinates)

  expect_s4_class(trimmed, "sparseMatrix")
  expect_equal(dim(trimmed), c(3L, 3L))
  expect_equal(Matrix::nnzero(trimmed), 0L)
})

test_that("sample-local graphs can share a precomputed global cutoff", {
  graph_a <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3), j = c(2, 1, 3, 2), x = 1, dims = c(3, 3)
  )
  graphs <- list(a = graph_a, b = graph_a)
  coordinates <- list(
    a = cbind(c(0, 1, 2), 0),
    b = cbind(c(0, 10, 20), 0)
  )
  cutoff <- cellcharter_global_distance_cutoff(graphs, coordinates, percentile = 50)
  expect_equal(cutoff, 5.5)
  trimmed <- lapply(names(graphs), function(id) {
    remove_long_links(graphs[[id]], coordinates[[id]], cutoff = cutoff)
  })
  expect_equal(sum(trimmed[[1]]), 4)
  expect_equal(sum(trimmed[[2]]), 0)
  expect_equal(attr(trimmed[[1]], "distance_cutoff"), cutoff)
})

test_that("graph inputs are validated before trimming", {
  symmetric <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(2, 1), x = 1, dims = c(2, 2)
  )
  asymmetric <- Matrix::sparseMatrix(i = 1, j = 2, x = 1, dims = c(2, 2))

  expect_error(remove_long_links(asymmetric, matrix(0, 2, 2)), "symmetric")
  expect_error(remove_long_links(symmetric, matrix(0, 3, 2)), "same number")
  expect_error(remove_long_links(symmetric, matrix(c(0, 0, Inf, 0), 2, 2, byrow = TRUE)), "finite")
  expect_error(remove_long_links(symmetric, matrix(0, 2, 2), -1), "between 0 and 100")
  expect_error(remove_long_links(symmetric, matrix(0, 2, 2), c(90, 99)), "one finite number")
})

test_that("Delaunay construction has deterministic small and degenerate policies", {
  empty <- delaunay_graph(matrix(numeric(), nrow = 0, ncol = 2))
  singleton <- delaunay_graph(matrix(c(2, 3), nrow = 1))
  pair <- delaunay_graph(rbind(c(0, 0), c(1, 0)))

  expect_equal(dim(empty), c(0L, 0L))
  expect_equal(Matrix::nnzero(singleton), 0L)
  expect_equal(as.matrix(pair), matrix(c(0, 1, 1, 0), 2, byrow = TRUE))
  duplicate <- delaunay_graph(rbind(c(0, 0), c(0, 0), c(1, 0), c(0, 1)))
  expect_identical(which(Matrix::rowSums(duplicate) == 0), 2L)
  expect_true(all(Matrix::rowSums(duplicate)[3:4] > 0))
  expect_error(
    delaunay_graph(rbind(c(0, 0), c(1, 1), c(2, 2))),
    "non-collinear"
  )
})

test_that("Delaunay construction supports tied axes without cross-sample edges", {
  skip_if_not_installed("geometry")
  coordinates <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  graph <- delaunay_graph(coordinates)

  expect_true(Matrix::isSymmetric(graph))
  expect_false(any(Matrix::diag(graph)))
  expect_true(all(Matrix::rowSums(graph) >= 2))

  combined <- Matrix::bdiag(graph, graph)
  expect_equal(sum(combined[1:4, 5:8]), 0)
  expect_equal(sum(combined[5:8, 1:4]), 0)
})

test_that("R Qhull matches the permanent Squidpy 1.3.0 graph fixture", {
  skip_if_not_installed("geometry")
  # Expected edges and percentile were generated with the pinned
  # scipy.spatial.Delaunay/vertex_neighbor_vertices implementation. The points
  # deliberately avoid co-circular ambiguity so this is a stable parity test.
  coordinates <- rbind(
    c(0, 0), c(2, 0), c(0, 3), c(2, 2),
    c(5, 1), c(4, 5), c(1, 6)
  )
  expected_edges <- rbind(
    c(1, 2), c(1, 3), c(1, 4), c(2, 4), c(2, 5), c(3, 4),
    c(3, 7), c(4, 5), c(4, 6), c(4, 7), c(5, 6), c(6, 7)
  )
  expected <- Matrix::sparseMatrix(
    i = c(expected_edges[, 1], expected_edges[, 2]),
    j = c(expected_edges[, 2], expected_edges[, 1]), x = 1,
    dims = c(7, 7)
  )
  observed <- delaunay_graph(coordinates)

  expect_identical(as.matrix(observed != 0), as.matrix(expected != 0))
  expect_equal(
    cellcharter_global_distance_cutoff(
      list(sample = observed), list(sample = coordinates), percentile = 75
    ),
    3.273096063992282,
    tolerance = 1e-12
  )
})
