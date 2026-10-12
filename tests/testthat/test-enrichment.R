test_that("cell-type enrichment matches CellCharter v0.2.0 proportions", {
  cell_type <- factor(c("A", "A", "B", "B", "A", "B"), levels = c("A", "B"))
  domain <- factor(c("D1", "D1", "D2", "D2", "D2", "D1"), levels = c("D1", "D2"))

  result <- .cc_enrichment(cell_type, domain)
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

  result <- .cc_nhood_enrichment(adjacency, cluster)
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

  first <- .cc_diff_nhood_enrichment(
    adjacency, cluster, condition, library,
    condition_groups = c("BALBc", "MRL"), n_perms = 25L, seed = 81L
  )
  second <- .cc_diff_nhood_enrichment(
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

test_that("v0.2.0 differential null exactly replays unequal 3-vs-6 libraries", {
  one_library <- function(center) {
    center_index <- match(center, c("A", "B", "C"))
    other_index <- setdiff(seq_len(3L), center_index)
    Matrix::sparseMatrix(
      i = c(center_index, other_index[[1L]], center_index, other_index[[2L]]),
      j = c(other_index[[1L]], center_index, other_index[[2L]], center_index),
      x = 1, dims = c(3L, 3L)
    )
  }
  centers <- c("A", "A", "B", "A", "B", "B", "C", "C", "C")
  adjacency <- do.call(Matrix::bdiag, lapply(centers, one_library))
  cluster <- factor(rep(c("A", "B", "C"), length(centers)),
                    levels = c("A", "B", "C"))
  library_ids <- c(
    "BALBc-1", "BALBc-2", "BALBc-3",
    "MRL-4", "MRL-5", "MRL-6", "MRL-7", "MRL-8", "MRL-9"
  )
  library <- factor(rep(library_ids, each = 3L), levels = library_ids)
  condition <- factor(
    rep(c(rep("BALBc", 3L), rep("MRL", 6L)), each = 3L),
    levels = c("BALBc", "MRL")
  )
  n_perms <- 37L
  seed <- 9182L

  result <- .cc_diff_nhood_enrichment(
    adjacency, cluster, condition, library,
    condition_groups = c("BALBc", "MRL"),
    n_perms = n_perms, seed = seed,
    permutation_mode = "cellcharter_v0.2.0"
  )

  enrichment_for <- function(ids) {
    rows <- which(as.character(library) %in% ids)
    .cc_nhood_enrichment(
      adjacency[rows, rows, drop = FALSE], droplevels(cluster[rows]),
      log2_fold_change = FALSE, only_inter = TRUE
    )$enrichment
  }
  manual_counts <- matrix(
    0L, nrow = nrow(result$effect), ncol = ncol(result$effect),
    dimnames = dimnames(result$effect)
  )
  set.seed(seed)
  for (permutation in seq_len(n_perms)) {
    # v0.2.0 draws all nine entries with replacement even though only the
    # BALBc-sized prefix controls membership. Unused suffix draws therefore
    # still advance the RNG before the next permutation.
    draws <- sample(library_ids, length(library_ids), replace = TRUE)
    selected <- unique(draws[seq_len(3L)])
    unselected <- setdiff(library_ids, selected)
    null_effect <- enrichment_for(unselected) - enrichment_for(selected)
    positive <- !is.na(result$effect) & result$effect > 0
    negative <- !is.na(result$effect) & result$effect < 0
    manual_counts[positive] <- manual_counts[positive] + as.integer(
      is.na(null_effect[positive]) | null_effect[positive] >= result$effect[positive]
    )
    manual_counts[negative] <- manual_counts[negative] + as.integer(
      is.na(null_effect[negative]) | null_effect[negative] <= result$effect[negative]
    )
  }
  manual_pvalue <- manual_counts / n_perms
  manual_pvalue[is.na(result$effect) | result$effect == 0] <- 0

  expect_identical(result$permutation_counts, manual_counts)
  expect_equal(result$pvalue, manual_pvalue, tolerance = 0)
  expect_identical(
    result$permutation_counts,
    matrix(
      c(0L, 4L, 7L, 2L, 0L, 2L, 5L, 2L, 0L),
      nrow = 3L, byrow = TRUE, dimnames = dimnames(result$effect)
    )
  )
  expect_identical(result$params$condition_level_order, c("BALBc", "MRL"))
  expect_identical(result$params$contrast_numerator, "BALBc")
  expect_identical(result$params$contrast_denominator, "MRL")
  expect_identical(result$params$library_counts, c(BALBc = 3L, MRL = 6L))
  expect_identical(result$params$library_order, library_ids)
  expect_identical(result$params$permutation_selected_draws, 3L)
  expect_identical(result$params$empirical_pvalue_denominator, n_perms)
  expect_identical(result$params$empirical_tail, "directional one-sided inclusive")
  expect_identical(result$params$seed, seed)

  # Supplying plot order in reverse must not silently change the inferential
  # null from a three-draw BALBc prefix to a six-draw MRL prefix.
  reverse_request <- .cc_diff_nhood_enrichment(
    adjacency, cluster, condition, library,
    condition_groups = c("MRL", "BALBc"),
    n_perms = n_perms, seed = seed,
    permutation_mode = "cellcharter_v0.2.0"
  )
  expect_equal(reverse_request$effect, result$effect)
  expect_identical(reverse_request$permutation_counts, result$permutation_counts)
  expect_equal(reverse_request$pvalue, result$pvalue, tolerance = 0)
  expect_identical(reverse_request$params$condition_level_order, c("BALBc", "MRL"))
  expect_identical(reverse_request$params$permutation_selected_draws, 3L)

  reverse_display <- CellCharterR:::.cc_diff_nhood_display(
    result, c("MRL", "BALBc")
  )
  expect_equal(reverse_display$effect, -result$effect)
  expect_identical(reverse_display$permutation_counts, result$permutation_counts)
  expect_identical(reverse_display$pvalue, result$pvalue)
  expect_identical(reverse_display$observed, result$observed)
  expect_identical(reverse_display$expected, result$expected)
  expect_identical(reverse_display$enrichment, result$enrichment)
  expect_identical(reverse_display$params$condition_level_order, c("BALBc", "MRL"))
  expect_identical(reverse_display$params$contrast_numerator, "BALBc")
  expect_identical(reverse_display$params$contrast_denominator, "MRL")
  expect_identical(reverse_display$params$display_groups, c("MRL", "BALBc"))
  expect_identical(reverse_display$params$display_effect_multiplier, -1)
  expect_true(reverse_display$params$inference_reused)
})

test_that("differential enrichment rejects split libraries and cross-library edges", {
  adjacency <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(2, 1), x = 1, dims = c(2, 2)
  )
  cluster <- factor(c("A", "B"))
  condition <- factor(c("BALBc", "MRL"))

  expect_error(
    .cc_diff_nhood_enrichment(
      adjacency, cluster, condition, factor(c("same", "same")), n_perms = 2L
    ),
    "exactly one condition"
  )
  expect_error(
    .cc_diff_nhood_enrichment(
      adjacency, cluster, condition, factor(c("lib1", "lib2")), n_perms = 2L
    ),
    "edges between physical libraries"
  )
})
