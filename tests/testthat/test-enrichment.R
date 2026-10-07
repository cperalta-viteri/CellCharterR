test_that("cell-type enrichment matches CellCharter v0.2.0 proportions", {
  cell_type <- factor(c("A", "A", "B", "B", "A", "B"), levels = c("A", "B"))
  domain <- factor(c("D1", "D1", "D2", "D2", "D2", "D1"), levels = c("D1", "D2"))

  result <- cellcharter_enrichment(cell_type, domain)
  expected_observed <- matrix(
    c(2 / 3, 1 / 3, 1 / 3, 2 / 3), nrow = 2L,
    dimnames = list(c("A", "B"), c("D1", "D2"))
  )
  expected_null <- matrix(
    0.5, nrow = 2L, ncol = 2L,
    dimnames = dimnames(expected_observed)
  )

  expect_equal(result$observed, expected_observed)
  expect_equal(result$expected, expected_null)
  expect_equal(result$enrichment, log2(expected_observed / expected_null))
  expect_equal(result$params$reference, "CellCharter 0.2.0")
})

test_that("asymmetric inter-domain expectation uses target cell fractions", {
  # Undirected edges: A1-B, A2-B, A2-C, and B-C.
  adjacency <- Matrix::sparseMatrix(
    i = c(1, 3, 2, 3, 2, 4, 3, 4),
    j = c(3, 1, 3, 2, 4, 2, 4, 3),
    x = 1, dims = c(4, 4)
  )
  cluster <- factor(c("A", "A", "B", "C"), levels = c("A", "B", "C"))

  result <- cellcharter_nhood_enrichment(adjacency, cluster)
  expected_observed <- matrix(
    c(NA, 2 / 3, 1 / 3,
      2 / 3, NA, 1 / 3,
      1 / 2, 1 / 2, NA),
    nrow = 3L, byrow = TRUE,
    dimnames = list(levels(cluster), levels(cluster))
  )
  expected_null <- matrix(
    rep(c(1 / 2, 1 / 4, 1 / 4), each = 3L),
    nrow = 3L, byrow = FALSE,
    dimnames = dimnames(expected_observed)
  )
  diag(expected_null) <- NA_real_

  expect_equal(result$observed, expected_observed)
  expect_equal(result$expected, expected_null)
  expect_equal(result$enrichment, expected_observed - expected_null)
  expect_false(isTRUE(all.equal(result$enrichment, t(result$enrichment))))
})

test_that("differential permutations are deterministic and library-level", {
  # Four disconnected three-cell libraries. In condition BALBc, A is the graph
  # center; in MRL, B is the center. No edge crosses a physical library.
  one_library <- function(center) {
    if (center == "A") {
      Matrix::sparseMatrix(
        i = c(1, 2, 1, 3), j = c(2, 1, 3, 1), x = 1, dims = c(3, 3)
      )
    } else {
      Matrix::sparseMatrix(
        i = c(2, 1, 2, 3), j = c(1, 2, 3, 2), x = 1, dims = c(3, 3)
      )
    }
  }
  adjacency <- Matrix::bdiag(
    one_library("A"), one_library("A"),
    one_library("B"), one_library("B")
  )
  cluster <- factor(rep(c("A", "B", "C"), 4L), levels = c("A", "B", "C"))
  library <- factor(rep(c("BALBc-1", "BALBc-2", "MRL-1", "MRL-2"), each = 3L))
  condition <- factor(rep(c("BALBc", "BALBc", "MRL", "MRL"), each = 3L),
                      levels = c("BALBc", "MRL"))

  first <- cellcharter_diff_nhood_enrichment(
    adjacency, cluster, condition, library,
    condition_groups = c("BALBc", "MRL"), n_perms = 25L, seed = 81L
  )
  second <- cellcharter_diff_nhood_enrichment(
    adjacency, cluster, condition, library,
    condition_groups = c("BALBc", "MRL"), n_perms = 25L, seed = 81L
  )

  expect_equal(first, second)
  expect_equal(first$effect, first$enrichment$BALBc - first$enrichment$MRL)
  expect_true(all(first$permutation_counts >= 0L & first$permutation_counts <= 25L))
  expect_true(all(first$pvalue >= 0 & first$pvalue <= 1))
  expect_equal(unname(diag(first$pvalue)), rep(0, 3L))
  expect_equal(first$params$library_counts, c(BALBc = 2, MRL = 2))
  expect_identical(first$params$library_unit, "sample")
  expect_equal(first$params$permutation_mode, "cellcharter_v0.2.0")
})

test_that("differential enrichment rejects split libraries and cross-library edges", {
  adjacency <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(2, 1), x = 1, dims = c(2, 2)
  )
  cluster <- factor(c("A", "B"))
  condition <- factor(c("BALBc", "MRL"))

  expect_error(
    cellcharter_diff_nhood_enrichment(
      adjacency, cluster, condition, factor(c("same", "same")), n_perms = 2L
    ),
    "exactly one condition"
  )
  expect_error(
    cellcharter_diff_nhood_enrichment(
      adjacency, cluster, condition, factor(c("lib1", "lib2")), n_perms = 2L
    ),
    "edges between physical libraries"
  )
})
