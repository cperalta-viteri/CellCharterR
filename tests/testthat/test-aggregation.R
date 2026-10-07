test_that("aggregation uses exclusive shortest-path rings", {
  path <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1, dims = c(5, 5)
  )
  x <- matrix(c(10, 20, 30, 40, 50), ncol = 1)
  observed <- as.matrix(cellcharter_aggregate(x, path, n_layers = 3))

  expected <- cbind(
    native = c(10, 20, 30, 40, 50),
    hop1 = c(20, 20, 30, 30, 0),
    hop2 = c(30, 40, 10, 20, 0),
    hop3 = c(40, 0, 0, 10, 0)
  )
  expect_equal(unname(observed), unname(expected))
})

test_that("mean aggregation matches the pinned CellCharter v0.2.0 fixture", {
  graph <- Matrix::Matrix(matrix(c(
    0, 1, 1, 1, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 1, 1, 0,
    1, 0, 0, 0, 1, 0, 0, 0, 1,
    1, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0
  ), 9, 9, byrow = TRUE), sparse = TRUE)
  x <- cbind(2^(0:8), 2^(8:0))

  hop1 <- cbind(
    c(7.5, 16.5, 193 / 3, 91, 4.5, 2, 66, 34, 8),
    c(60, 132, 262 / 3, 91, 144, 128, 33, 34, 32)
  )
  hop2 <- cbind(
    c(120, 28 / 3, 26 / 3, 3, 262 / 3, 1, 1, 1, 8.5),
    c(3.75, 112 / 3, 176 / 3, 96, 193 / 3, 256, 256, 256, 136)
  )
  observed <- as.matrix(cellcharter_aggregate(x, graph, 2))

  expect_equal(observed, cbind(x, hop1, hop2), tolerance = 1e-6)
})

test_that("cycles do not reintroduce focal or previously visited cells", {
  triangle <- Matrix::Matrix(
    matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3, byrow = TRUE),
    sparse = TRUE
  )
  observed <- as.matrix(cellcharter_aggregate(matrix(1:3, ncol = 1), triangle, 3))

  expect_equal(observed[, 1], 1:3)
  expect_equal(observed[, 2], c(2.5, 2, 1.5))
  expect_equal(observed[, 3:4, drop = FALSE], matrix(0, 3, 2))
})

test_that("zero-neighbor rows are zero and samples remain isolated", {
  paired <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4), j = c(2, 1, 4, 3), x = 1,
    dims = c(5, 5)
  )
  x <- matrix(c(1, 2, 100, 200, 999), ncol = 1)
  observed <- as.matrix(cellcharter_aggregate(x, paired, 2))

  expect_equal(observed[, 2], c(2, 1, 200, 100, 0))
  expect_equal(observed[, 3], rep(0, 5))
})

test_that("explicit layer vectors select only requested exclusive rings", {
  path <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3), j = c(2, 1, 3, 2), x = c(4, 4, 8, 8),
    dims = c(3, 3)
  )
  Matrix::diag(path) <- 10
  x <- cbind(a = 1:3, b = 11:13)
  observed <- as.matrix(cellcharter_aggregate(x, path, n_layers = c(2, 0, 2)))

  expect_equal(dim(observed), c(3L, 4L))
  expect_equal(observed[, 1:2], x)
  expect_equal(observed[, 3:4], rbind(x[3, ], c(0, 0), x[1, ]))
})

test_that("aggregation validates topology, dimensions, layers, and dense values", {
  graph <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(2, 1), x = 1, dims = c(2, 2)
  )
  asymmetric <- Matrix::sparseMatrix(i = 1, j = 2, x = 1, dims = c(2, 2))

  expect_error(cellcharter_aggregate(matrix(1, 3, 1), graph), "one row and column")
  expect_error(cellcharter_aggregate(matrix(1, 2, 1), asymmetric), "symmetric")
  expect_error(cellcharter_aggregate(matrix(1, 2, 1), graph, c(1, 2)), "include 0")
  expect_error(cellcharter_aggregate(matrix(1, 2, 1), graph, -1), "non-negative")
  expect_error(cellcharter_aggregate(matrix(1, 2, 1), graph, 1.5), "whole numbers")
  expect_error(cellcharter_aggregate(matrix(c(1, NA), 2, 1), graph), "finite numeric")
})

test_that("HDF5 aggregation supports explicit layers and interleaved sample rows", {
  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("S4Arrays")
  skip_if_not_installed("IRanges")

  x <- cbind(a = c(1, 100, 2, 200), b = c(11, 110, 12, 120))
  pair <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(2, 1), x = 1, dims = c(2, 2)
  )
  graphs <- list(A = pair, B = pair)
  rows <- list(A = c(1, 3), B = c(2, 4))
  filepath <- tempfile(fileext = ".h5")
  on.exit(unlink(filepath), add = TRUE)

  observed <- cellcharter_aggregate_hdf5(
    x, graphs, rows, filepath, n_layers = c(0, 1)
  )
  expected <- rbind(
    c(1, 11, 2, 12), c(100, 110, 200, 120),
    c(2, 12, 1, 11), c(200, 120, 100, 110)
  )

  expect_equal(dim(observed), c(4L, 4L))
  expect_equal(as.matrix(observed), expected)
})

test_that("direct aggregation accepts an HDF5-backed native representation", {
  skip_if_not_installed("HDF5Array")
  filepath <- tempfile(fileext = ".h5")
  on.exit(unlink(filepath), add = TRUE)
  x <- HDF5Array::writeHDF5Array(matrix(1:8, nrow = 4), filepath = filepath)
  pair <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 1, 3, 4, 4, 3),
    j = c(2, 1, 1, 2, 4, 3, 3, 4), x = 1, dims = c(4, 4)
  )

  observed <- cellcharter_aggregate(x, pair, n_layers = 1L)

  expect_s4_class(observed, "DelayedMatrix")
  expect_equal(as.matrix(observed)[, 1:2], as.matrix(x))
  expect_equal(as.matrix(observed)[, 3:4], as.matrix(x)[c(2, 1, 4, 3), ])
})
