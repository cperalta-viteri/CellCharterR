.result_freshness_spe <- function() {
  sample <- rep(paste0("s", 1:4), each = 4L)
  condition <- factor(
    rep(c("healthy", "healthy", "disease", "disease"), each = 4L),
    levels = c("healthy", "disease")
  )
  ids <- paste0("spot", seq_along(sample))
  spatial <- cbind(x = seq_along(ids), y = 0)
  rownames(spatial) <- ids
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(
      0, nrow = 1L, ncol = length(ids), dimnames = list("g1", ids)
    )),
    colData = S4Vectors::DataFrame(
      sample = sample,
      condition = condition,
      group = rep(c("T", "B"), length.out = length(ids)),
      domain = factor(rep(c("D1", "D2"), each = 2L, length.out = length(ids))),
      row.names = ids
    ),
    spatialCoords = spatial
  )
  from <- unlist(lapply(seq(1L, length(ids), by = 4L), function(i) {
    c(i, i + 1L, i + 1L, i + 2L, i + 2L, i + 3L)
  }), use.names = FALSE)
  to <- unlist(lapply(seq(1L, length(ids), by = 4L), function(i) {
    c(i + 1L, i, i + 2L, i + 1L, i + 3L, i + 2L)
  }), use.names = FALSE)
  SingleCellExperiment::colPair(spe, "graph") <- S4Vectors::SelfHits(
    from, to, nnode = ncol(spe), value = rep(1, length(from))
  )
  spe
}

.result_freshness_zonation_spe <- function() {
  block <- rep(paste0("p", 1:4), each = 6L)
  zone <- rep(rep(c("Pre", "Early", "Late"), each = 2L), 4L)
  scores <- matrix(
    -1, nrow = length(block), ncol = 4L,
    dimnames = list(NULL, c("Pre-EMT", "Early-EMT", "Late-EMT", "Immune"))
  )
  scores[cbind(seq_len(nrow(scores)), match(zone, c("Pre", "Early", "Late")))] <- 3
  scores[, "Immune"] <- ifelse(zone == "Late", 2, -0.5)
  ids <- paste0("z", seq_along(block))
  rownames(scores) <- ids
  spatial <- cbind(x = seq_along(ids), y = 0)
  rownames(spatial) <- ids
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(
      0, 1L, length(ids), dimnames = list("g1", ids)
    )),
    colData = S4Vectors::DataFrame(block = block, row.names = ids),
    spatialCoords = spatial
  )
  SingleCellExperiment::reducedDim(spe, "AUCell") <- scores
  spe
}

test_that("enrichment results are bound to annotations and parameters", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_enrichment(
    .result_freshness_spe(), group_col = "group", domain_col = "domain",
    progress = FALSE
  )
  record <- S4Vectors::metadata(spe)$CellCharterR$results$enrichment
  expect_match(record$parameters_sha256, "^[0-9a-f]{64}$")
  expect_named(record$input_fingerprints, c("group", "domain"))
  expect_s3_class(cellcharter_plot_enrichment(spe), "ggplot")

  changed_group <- spe
  changed_group$group[[1L]] <- "changed"
  expect_error(cellcharter_plot_enrichment(changed_group), "stale.*group")

  changed_domain <- spe
  changed_domain$domain[[1L]] <- "D2"
  expect_error(cellcharter_plot_enrichment(changed_domain), "stale.*domain")

  changed_parameters <- spe
  metadata <- S4Vectors::metadata(changed_parameters)$CellCharterR
  metadata$results$enrichment$parameters$log2_fold_change <- FALSE
  S4Vectors::metadata(changed_parameters)$CellCharterR <- metadata
  expect_error(cellcharter_plot_enrichment(changed_parameters),
               "stale.*parameters")
})

test_that("complete legacy result records remain readable after upgrade", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_enrichment(
    .result_freshness_spe(), group_col = "group", domain_col = "domain",
    progress = FALSE
  )
  metadata <- S4Vectors::metadata(spe)$CellCharterR
  legacy <- metadata$results$enrichment
  legacy$n_cells <- NULL
  legacy$cell_order_sha256 <- NULL
  legacy$input_fingerprints <- NULL
  legacy$parameters_sha256 <- NULL
  metadata$results$enrichment <- legacy
  S4Vectors::metadata(spe)$CellCharterR <- metadata

  expect_warning(
    plot <- cellcharter_plot_enrichment(spe),
    "unverified legacy freshness"
  )
  expect_s3_class(plot, "ggplot")

  malformed <- spe
  metadata <- S4Vectors::metadata(malformed)$CellCharterR
  metadata$results$enrichment$n_cells <- ncol(malformed)
  S4Vectors::metadata(malformed)$CellCharterR <- metadata
  expect_error(
    cellcharter_plot_enrichment(malformed),
    "incomplete identity metadata"
  )
})

test_that("neighborhood results detect topology, weight, and domain changes", {
  spe <- cellcharter_nhood_enrichment(
    .result_freshness_spe(), graph_name = "graph", domain_col = "domain"
  )
  expect_named(
    S4Vectors::metadata(spe)$CellCharterR$results$nhood_enrichment$
      input_fingerprints,
    c("graph", "domain")
  )

  changed_weight <- spe
  pairs <- SingleCellExperiment::colPair(changed_weight, "graph")
  S4Vectors::mcols(pairs)[[1L]][[1L]] <- 2
  SingleCellExperiment::colPair(changed_weight, "graph") <- pairs
  expect_error(
    .cc_spe_require_result(changed_weight, "nhood_enrichment"),
    "stale.*graph"
  )

  changed_topology <- spe
  pairs <- SingleCellExperiment::colPair(changed_topology, "graph")
  SingleCellExperiment::colPair(changed_topology, "graph") <- pairs[-1L]
  expect_error(
    .cc_spe_require_result(changed_topology, "nhood_enrichment"),
    "stale.*graph"
  )

  changed_domain <- spe
  changed_domain$domain[[1L]] <- "D2"
  expect_error(
    .cc_spe_require_result(changed_domain, "nhood_enrichment"),
    "stale.*domain"
  )
})

test_that("differential plots reject changed graph and design variables", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_diff_nhood_enrichment(
    .result_freshness_spe(), graph_name = "graph", domain_col = "domain",
    condition_col = "condition", sample_col = "sample",
    condition_groups = c("healthy", "disease"), n_perms = 3L,
    permutation_mode = "shuffle", progress = FALSE
  )
  expect_named(
    S4Vectors::metadata(spe)$CellCharterR$results$diff_nhood_enrichment$
      input_fingerprints,
    c("graph", "domain", "condition", "sample")
  )

  changed_condition <- spe
  changed_condition$condition[[1L]] <- "disease"
  expect_error(cellcharter_plot_graph_diff(changed_condition),
               "stale.*condition")

  changed_sample <- spe
  changed_sample$sample[[1L]] <- "another_sample"
  expect_error(cellcharter_plot_graph_diff(changed_sample), "stale.*sample")

  changed_graph <- spe
  pairs <- SingleCellExperiment::colPair(changed_graph, "graph")
  S4Vectors::mcols(pairs)[[1L]][[1L]] <- 3
  SingleCellExperiment::colPair(changed_graph, "graph") <- pairs
  expect_error(cellcharter_plot_graph_diff(changed_graph), "stale.*graph")
})

test_that("zonation plots reject changed scores, blocks, and parameters", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_zonation(
    .result_freshness_zonation_spe(), block = "block",
    min_spots_per_zone = 1L, sample_size = 1L, n_iterations = 2L
  )
  expect_s3_class(
    cellcharter_plot_zonation(spe, compare = "none", show_individual = FALSE),
    "ggplot"
  )

  changed_scores <- spe
  scores <- SingleCellExperiment::reducedDim(changed_scores, "AUCell")
  scores[1L, 1L] <- scores[1L, 1L] + 0.5
  SingleCellExperiment::reducedDim(changed_scores, "AUCell") <- scores
  expect_error(
    cellcharter_plot_zonation(
      changed_scores, compare = "none", show_individual = FALSE
    ),
    "stale.*scores"
  )

  changed_block <- spe
  changed_block$block[[1L]] <- "another_patient"
  expect_error(
    cellcharter_plot_zonation(
      changed_block, compare = "none", show_individual = FALSE
    ),
    "stale.*block"
  )

  changed_parameters <- spe
  metadata <- S4Vectors::metadata(changed_parameters)$CellCharterR
  metadata$results$zonation$parameters$sample_size <- 2L
  S4Vectors::metadata(changed_parameters)$CellCharterR <- metadata
  expect_error(
    cellcharter_plot_zonation(
      changed_parameters, compare = "none", show_individual = FALSE
    ),
    "stale.*parameters"
  )
})
