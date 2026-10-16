.diff_network_fixture <- function() {
  effect <- matrix(
    c(NA, 4, -2, 3, NA, -1, 0, Inf, NA),
    nrow = 3L, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  pvalue <- matrix(
    c(0, 0.01, 0.04, 0.02, 0, 0.20, 0, 0.9, 0),
    nrow = 3L, byrow = TRUE, dimnames = dimnames(effect)
  )
  list(
    effect = effect, pvalue = pvalue,
    permutation_counts = matrix(seq_len(9L), 3L),
    params = list(
      condition_level_order = c("healthy", "disease"),
      contrast_numerator = "healthy", contrast_denominator = "disease",
      empirical_tail = "directional one-sided inclusive"
    )
  )
}

test_that("differential network data preserves direction and panel signs", {
  edges <- .cc_diff_nhood_network_data(.diff_network_fixture())
  keys <- paste(edges$source_domain, edges$neighbor_domain, sep = "->")

  expect_setequal(keys, c("A->B", "A->C", "B->A", "B->C"))
  expect_equal(edges$signed_effect[match("A->B", keys)], 4)
  expect_equal(edges$signed_effect[match("B->A", keys)], 3)
  expect_equal(edges$signed_effect[match("A->C", keys)], -2)
  expect_equal(edges$signed_effect[match("B->C", keys)], -1)
  expect_true(all(edges$display_effect == abs(edges$signed_effect)))
  expect_true(all(as.character(edges$favoured_condition[edges$signed_effect > 0]) ==
                    "healthy"))
  expect_true(all(as.character(edges$favoured_condition[edges$signed_effect < 0]) ==
                    "disease"))
  expect_false(any(keys %in% c("A->A", "B->B", "C->C", "C->A", "C->B")))
})

test_that("p adjustment precedes display filtering and top_n is per panel", {
  all_edges <- .cc_diff_nhood_network_data(
    .diff_network_fixture(), p_adjust_method = "BH"
  )
  filtered <- .cc_diff_nhood_network_data(
    .diff_network_fixture(), p_adjust_method = "BH",
    effect_threshold = 1.5, top_n = 1L
  )
  all_keys <- paste(all_edges$source_domain, all_edges$neighbor_domain, sep = "->")
  filtered_keys <- paste(filtered$source_domain, filtered$neighbor_domain, sep = "->")
  expect_setequal(filtered_keys, c("A->B", "A->C"))
  expect_equal(
    filtered$adjusted_pvalue,
    all_edges$adjusted_pvalue[match(filtered_keys, all_keys)]
  )
  expect_identical(
    stats::setNames(as.integer(table(filtered$favoured_condition)),
                    levels(filtered$favoured_condition)),
    c(healthy = 1L, disease = 1L)
  )
})

test_that("differential network arguments and matrices are validated", {
  result <- .diff_network_fixture()
  bad <- result
  dimnames(bad$pvalue)[[2L]] <- c("B", "A", "C")
  expect_error(.cc_diff_nhood_network_data(bad), "identical.*dimnames")
  bad <- result
  bad$pvalue[1, 2] <- 2
  expect_error(.cc_diff_nhood_network_data(bad), "between zero and one")
  bad <- result
  bad$params$condition_level_order <- "healthy"
  expect_error(.cc_diff_nhood_network_data(bad), "two valid")
  expect_error(.cc_diff_nhood_network_data(result, top_n = 0), "positive whole")
  expect_error(.cc_diff_nhood_network_data(result, effect_threshold = -1),
               "non-negative")
})

test_that("plot_graph_diff uses stable shared coordinates and black arrows", {
  skip_if_not_installed("ggplot2")
  plot <- cellcharter_plot_graph_diff(
    .diff_network_fixture(), top_n = 1L, edge_colours = "black"
  )
  expect_s3_class(plot, "ggplot")
  edges <- attr(plot, "cellcharter_network_data")
  nodes <- attr(plot, "cellcharter_node_coordinates")
  expect_equal(unique(as.character(nodes$favoured_condition)),
               c("healthy", "disease"))
  wide <- split(nodes[, c("domain", "x", "y")], nodes$favoured_condition)
  rownames(wide[[1L]]) <- rownames(wide[[2L]]) <- NULL
  expect_equal(wide[[1L]], wide[[2L]])
  expect_true(all(is.finite(edges$x)) && all(is.finite(edges$xend)))
  expect_identical(unname(plot$scales$get_scales("colour")$palette(2)),
                   c("black", "black"))
  expect_gte(plot$theme$text$size, 15)
})

test_that("plot_graph_diff extracts and validates a stored experiment result", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
  identity <- .cc_validate_sce_identity(spe)
  metadata <- S4Vectors::metadata(spe)$CellCharterR
  metadata$results$diff_nhood_enrichment <- list(
    value = .diff_network_fixture(), n_cells = ncol(spe),
    cell_order_sha256 = identity$cell_order_sha256, parameters = list()
  )
  S4Vectors::metadata(spe)$CellCharterR <- metadata
  direct <- .cc_diff_nhood_network_data(.diff_network_fixture())
  expect_warning(
    plot <- cellcharter_plot_graph_diff(spe),
    "unverified legacy freshness"
  )
  stored <- attr(plot, "cellcharter_network_data")
  for (column in names(direct)) {
    expect_equal(stored[[column]], direct[[column]], info = column)
  }

  stale <- spe[, rev(seq_len(ncol(spe)))]
  expect_error(cellcharter_plot_graph_diff(stale), "stale")
})

test_that("reciprocal curvature is deterministic", {
  result <- .diff_network_fixture()
  first <- .cc_diff_nhood_network_data(result)
  second <- .cc_diff_nhood_network_data(result)
  first_keys <- paste(first$source_domain, first$neighbor_domain, sep = "->")
  expect_equal(first$curvature, second$curvature)
  expect_equal(first$curvature[first_keys %in% c("A->B", "B->A")],
               c(0.18, 0.18))
})
