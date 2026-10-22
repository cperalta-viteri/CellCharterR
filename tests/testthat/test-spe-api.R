.tiny_spe_fixture <- function() {
  sample <- rep(c("A1", "A2", "B1", "B2"), each = 4L)
  condition <- rep(c("healthy", "healthy", "lupus", "lupus"), each = 4L)
  local_xy <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1.2, 1.1))
  spatial <- do.call(rbind, lapply(seq_len(4L), function(i) {
    sweep(local_xy, 2L, c(i * 10, i * 3), "+")
  }))
  ids <- paste0("spot", seq_along(sample))
  representation <- cbind(
    x = rep(c(-2, -1, 1, 2), 4L),
    y = rep(c(-1, 1, -1, 1), 4L)
  )
  rownames(spatial) <- rownames(representation) <- ids
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(protein = matrix(seq_len(32L), nrow = 2L,
                                   dimnames = list(c("m1", "m2"), ids))),
    colData = S4Vectors::DataFrame(
      sample = sample, condition = factor(condition, levels = c("healthy", "lupus")),
      cell_type = rep(c("T", "B"), 8L), row.names = ids
    ),
    reducedDims = S4Vectors::SimpleList(X_trVAE = representation),
    spatialCoords = spatial,
    sample_id = sample
  )
  spe
}

.tiny_fit <- function(x, n_clusters, seed, ...) list(k = n_clusters, seed = seed)
.tiny_predict <- function(model, newdata, ...) {
  (seq_len(nrow(newdata)) - 1L) %% model$k
}
.tiny_nll <- function(model) as.numeric(model$k)

test_that("colData domain labels need no CellCharterR clustering stage", {
  spe <- .tiny_spe_fixture()
  spe$python_domain <- factor(rep(c("domain_b", "domain_a"), 8L))
  enriched <- cellcharter_enrichment(
    spe, group_col = "cell_type", domain_col = "python_domain"
  )
  expect_true("enrichment" %in% names(S4Vectors::metadata(enriched)$CellCharterR$results))
})

test_that("enrichment accepts repeated observation labels within a sample", {
  spe <- .tiny_spe_fixture()
  spe$python_domain <- factor(rep(c("domain_b", "domain_a"), 8L))
  colnames(spe) <- rep(c("slide_A", "slide_B"), each = 8L)
  coordinates <- SpatialExperiment::spatialCoords(spe)
  rownames(coordinates) <- colnames(spe)
  SpatialExperiment::spatialCoords(spe) <- coordinates

  enriched <- cellcharter_enrichment(
    spe, group_col = "cell_type", domain_col = "python_domain",
    progress = FALSE
  )

  expect_true("enrichment" %in%
                names(S4Vectors::metadata(enriched)$CellCharterR$results))
})

test_that("imported colPair graphs need no CellCharterR graph stage", {
  spe <- .tiny_spe_fixture()
  spe$python_domain <- factor(rep(c("domain_a", "domain_b"), 8L))
  from <- unlist(lapply(seq(1L, ncol(spe), by = 4L), function(i) {
    c(i, i + 1L, i + 1L, i + 2L, i + 2L, i + 3L)
  }), use.names = FALSE)
  to <- unlist(lapply(seq(1L, ncol(spe), by = 4L), function(i) {
    c(i + 1L, i, i + 2L, i + 1L, i + 3L, i + 2L)
  }), use.names = FALSE)
  SingleCellExperiment::colPair(spe, "spatial_connectivities") <-
    S4Vectors::SelfHits(from, to, nnode = ncol(spe), value = 1)
  observed <- cellcharter_diff_nhood_enrichment(
    spe, graph_name = "spatial_connectivities",
    domain_col = "python_domain", condition_col = "condition",
    sample_col = "sample", condition_groups = c("healthy", "lupus"),
    n_perms = 3L, permutation_mode = "shuffle", progress = FALSE
  )
  expect_true("diff_nhood_enrichment" %in%
                names(S4Vectors::metadata(observed)$CellCharterR$results))
})

test_that("domain labels can be read from reducedDims", {
  spe <- .tiny_spe_fixture()
  labels <- matrix(rep(0:1, 8L), ncol = 1L,
                   dimnames = list(colnames(spe), "python_domain"))
  SingleCellExperiment::reducedDim(spe, "python_labels") <- labels

  by_column <- cellcharter_enrichment(
    spe, group_col = "cell_type", domain_col = "python_domain"
  )
  by_dimension <- cellcharter_enrichment(
    spe, group_col = "cell_type", domain_col = "python_labels",
    result_name = "from_dimension"
  )
  expect_true("enrichment" %in% names(S4Vectors::metadata(by_column)$CellCharterR$results))
  expect_true("from_dimension" %in% names(S4Vectors::metadata(by_dimension)$CellCharterR$results))
})

test_that("public object API completes the tiny workflow without sample leakage", {
  skip_if_not_installed("geometry")
  checkpoint <- tempfile("spe-autok-")
  spe <- cellcharter_workflow(
    .tiny_spe_fixture(), group_col = "cell_type", candidate_k = c(2L, 3L),
    checkpoint_dir = checkpoint, n_layers = 1L, n_perms = 4L, seed = 17L,
    autok_args = list(
      max_runs = 2L, convergence_tol = 0,
      fit_fun = .tiny_fit, predict_fun = .tiny_predict, nll_fun = .tiny_nll
    )
  )

  expect_s4_class(spe, "SpatialExperiment")
  expect_true("cellcharter_graph" %in% SingleCellExperiment::colPairNames(spe))
  expect_true("X_cellcharter" %in% SingleCellExperiment::reducedDimNames(spe))
  expect_identical(
    rownames(SingleCellExperiment::reducedDim(spe, "X_cellcharter")), colnames(spe)
  )
  expect_true("cellcharter_domain" %in% colnames(SingleCellExperiment::colData(spe)))

  graph <- SingleCellExperiment::colPair(spe, "cellcharter_graph")
  sample <- as.character(SingleCellExperiment::colData(spe)$sample)
  expect_true(all(sample[S4Vectors::from(graph)] == sample[S4Vectors::to(graph)]))
  results <- S4Vectors::metadata(spe)$CellCharterR$results
  expect_true(all(c("enrichment", "nhood_enrichment", "diff_nhood_enrichment") %in%
                    names(results)))
  expect_identical(results$diff_nhood_enrichment$value$params$library_unit, "sample")
})

test_that("SPE clustering records and enforces representation content identity", {
  skip_if_not_installed("geometry")
  checkpoint <- tempfile("spe-content-identity-")
  on.exit(unlink(checkpoint, recursive = TRUE), add = TRUE)
  spe <- cellcharter_graph(.tiny_spe_fixture(), progress = FALSE)
  spe <- cellcharter_aggregate(spe, n_layers = 1L, progress = FALSE)
  clustered <- cellcharter_cluster(
    spe, candidate_k = c(2L, 3L), checkpoint_dir = checkpoint,
    max_runs = 2L, convergence_tol = 0,
    fit_fun = .tiny_fit, predict_fun = .tiny_predict,
    nll_fun = .tiny_nll, progress = FALSE
  )
  stage <- S4Vectors::metadata(clustered)$CellCharterR$stages[[
    "clustering:cellcharter_domain"
  ]]
  expect_match(stage$parameters$representation_sha256, "^[0-9a-f]{64}$")
  expect_match(stage$parameters$autok_input_stage_sha256, "^[0-9a-f]{64}$")

  changed <- spe
  value <- SingleCellExperiment::reducedDim(changed, "X_cellcharter")
  value[[1L]] <- value[[1L]] + 0.25
  SingleCellExperiment::reducedDim(changed, "X_cellcharter") <- value
  expect_error(
    cellcharter_cluster(
      changed, candidate_k = c(2L, 3L), checkpoint_dir = checkpoint,
      max_runs = 2L, convergence_tol = 0,
      fit_fun = .tiny_fit, predict_fun = .tiny_predict,
      nll_fun = .tiny_nll, progress = FALSE
    ),
    "representation contents"
  )
})

test_that("native SpatialExperiment coordinates survive the complete workflow", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("geometry")
  checkpoint <- tempfile("spe-autok-")
  spe <- .tiny_spe_fixture()
  expected_coordinates <- SpatialExperiment::spatialCoords(spe)
  result <- cellcharter_workflow(
    spe, group_col = "cell_type", candidate_k = c(2L, 3L),
    checkpoint_dir = checkpoint, n_layers = 1L, n_perms = 4L, seed = 17L,
    autok_args = list(
      max_runs = 2L, convergence_tol = 0,
      fit_fun = .tiny_fit, predict_fun = .tiny_predict, nll_fun = .tiny_nll
    )
  )

  expect_s4_class(result, "SpatialExperiment")
  expect_equal(SpatialExperiment::spatialCoords(result), expected_coordinates, tolerance = 0)
  expect_false("spatial" %in% SingleCellExperiment::reducedDimNames(result))
  expect_identical(
    S4Vectors::metadata(result)$CellCharterR$stages[["graph:cellcharter_graph"]]$
      parameters$coordinate_source,
    "spatialCoords"
  )
  graph <- SingleCellExperiment::colPair(result, "cellcharter_graph")
  sample <- as.character(SummarizedExperiment::colData(result)$sample)
  expect_true(all(sample[S4Vectors::from(graph)] == sample[S4Vectors::to(graph)]))

  index <- c(1:3, 5:7, 9:11, 13:15)
  subset <- result[, index]
  expect_s4_class(subset, "SpatialExperiment")
  expect_equal(
    SpatialExperiment::spatialCoords(subset), expected_coordinates[index, , drop = FALSE],
    tolerance = 0
  )
  subset_graph <- SingleCellExperiment::colPair(subset, "cellcharter_graph")
  expect_true(all(S4Vectors::from(subset_graph) <= ncol(subset)))
  expect_true(all(S4Vectors::to(subset_graph) <= ncol(subset)))
})

test_that("disk-backed numeric reduced dimensions pass validation without realization", {
  spe <- .tiny_spe_fixture()
  delayed <- DelayedArray::DelayedArray(
    SingleCellExperiment::reducedDim(spe, "X_trVAE")
  )
  SingleCellExperiment::reducedDim(spe, "X_delayed") <- delayed

  observed <- .cc_spe_reduced_dim(spe, "X_delayed")

  expect_s4_class(observed, "DelayedMatrix")
  expect_identical(DelayedArray::type(observed), "double")
  expect_equal(dim(observed), c(ncol(spe), 2L))
})

test_that("enrichment plotting accepts stored results and direct result lists", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
  result <- S4Vectors::metadata(spe)$CellCharterR$results$enrichment$value

  # The regenerated packaged workflow carries current freshness metadata, so
  # plotting its stored result must not be treated as a legacy handoff.
  expect_no_warning(
    from_spe <- cellcharter_plot_enrichment(spe)
  )
  from_result <- cellcharter_plot_enrichment(result, title = "Tiny enrichment")

  expect_s3_class(from_spe, "ggplot")
  expect_s3_class(from_result, "ggplot")
  expected_positive <- sum(
    is.finite(result$enrichment) & result$enrichment > 0
  )
  expect_equal(nrow(from_result$data), expected_positive)
  expect_true(all(from_result$data$enrichment > 0))
  expect_equal(from_result$labels$x, "Cell type")
  expect_equal(from_result$labels$y, "CellCharter domain")
  expect_identical(from_result$labels$title, "Tiny enrichment")
})

test_that("enrichment plotting suppresses depletion and caps positive display values", {
  skip_if_not_installed("ggplot2")
  values <- matrix(
    c(-Inf, -3, -0.3, 0, 0.3, 2), nrow = 3L,
    dimnames = list(c("A", "B", "C"), c("0", "1"))
  )

  plotted <- cellcharter_plot_enrichment(
    values, size_threshold = 0.58, colour_threshold = 0.58,
    cluster_domains = FALSE
  )

  expect_equal(sort(plotted$data$enrichment), c(0.3, 2))
  expect_equal(sort(plotted$data$point_size), c(0.3, 0.58))
  expect_equal(sort(plotted$data$colour_value), c(0.3, 0.58))
  expect_identical(levels(plotted$data$group), rownames(values))
  expect_identical(levels(plotted$data$domain), rev(colnames(values)))
})

test_that("enrichment plotting never labels raw ratios as log2 fold changes", {
  skip_if_not_installed("ggplot2")
  ratio_result <- list(
    enrichment = matrix(
      c(0.5, 2), nrow = 1L,
      dimnames = list("program" = NULL, "domain" = NULL)
    ),
    params = list(log2_fold_change = FALSE)
  )
  rownames(ratio_result$enrichment) <- "program"
  colnames(ratio_result$enrichment) <- c("A", "B")
  expect_error(
    cellcharter_plot_enrichment(ratio_result),
    "requires log2 fold-change enrichment", fixed = TRUE
  )
})

test_that("enrichment plotting clusters program profiles with hclust", {
  skip_if_not_installed("ggplot2")
  values <- matrix(
    c(2, -1, 0.2, -2, 1.5, 0.1, 0.1, 0.2, 2),
    nrow = 3L, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("D1", "D2", "D3"))
  )
  states <- ifelse(values < 0, -1, ifelse(values < 0.58, 0, 1))
  expected <- rownames(values)[stats::hclust(
    stats::dist(states), method = "complete"
  )$order]
  plotted <- cellcharter_plot_enrichment(
    values, cluster_domains = FALSE, cluster_programs = TRUE
  )
  expect_identical(levels(plotted$data$group), expected)
  expect_error(
    cellcharter_plot_enrichment(values, cluster_programs = NA),
    "cluster_programs"
  )
})

test_that("object API rejects positional identity and stale subset metadata", {
  skip_if_not_installed("geometry")
  spe <- .tiny_spe_fixture()
  expect_error(
    .cc_align_observation_ids(rev(colnames(spe)), colnames(spe)),
    "order mismatch", fixed = TRUE
  )

  spe <- cellcharter_graph(spe)
  subset <- spe[, c(1:3, 5:7, 9:11, 13:15)]
  pairs <- SingleCellExperiment::colPair(subset, "cellcharter_graph")
  expect_true(all(S4Vectors::from(pairs) <= ncol(subset)))
  expect_true(all(S4Vectors::to(pairs) <= ncol(subset)))
  expect_error(
    cellcharter_aggregate(subset),
    "stale after spot reordering or subsetting", fixed = TRUE
  )
})

test_that("graph stages reject changed coordinates, samples, and topology", {
  skip_if_not_installed("geometry")
  original <- cellcharter_graph(.tiny_spe_fixture(), progress = FALSE)

  changed_coordinates <- original
  coordinates <- SpatialExperiment::spatialCoords(changed_coordinates)
  coordinates[1L, 1L] <- coordinates[1L, 1L] + 0.25
  SpatialExperiment::spatialCoords(changed_coordinates) <- coordinates
  expect_error(
    cellcharter_aggregate(changed_coordinates, n_layers = 1L, progress = FALSE),
    "`coordinates`.*changed"
  )

  changed_samples <- original
  changed_samples$sample[[1L]] <- "A2"
  expect_error(
    cellcharter_aggregate(changed_samples, n_layers = 1L, progress = FALSE),
    "`sample`.*changed"
  )

  changed_graph <- original
  pairs <- SingleCellExperiment::colPair(changed_graph, "cellcharter_graph")
  SingleCellExperiment::colPair(changed_graph, "cellcharter_graph") <- pairs[-1L]
  expect_error(
    cellcharter_aggregate(changed_graph, n_layers = 1L, progress = FALSE),
    "`graph`.*changed"
  )
})

test_that("aggregation and clustering stages reject changed scientific inputs", {
  skip_if_not_installed("geometry")
  spe <- cellcharter_graph(.tiny_spe_fixture(), progress = FALSE)
  spe <- cellcharter_aggregate(spe, n_layers = 1L, progress = FALSE)

  changed_source <- spe
  source <- SingleCellExperiment::reducedDim(changed_source, "X_trVAE")
  source[1L, 1L] <- source[1L, 1L] + 0.25
  SingleCellExperiment::reducedDim(changed_source, "X_trVAE") <- source
  expect_error(
    .cc_spe_require_stage(changed_source, "aggregation:X_cellcharter"),
    "`representation` representation contents changed"
  )

  changed_output <- spe
  output <- SingleCellExperiment::reducedDim(changed_output, "X_cellcharter")
  output[1L, 1L] <- output[1L, 1L] + 0.25
  SingleCellExperiment::reducedDim(changed_output, "X_cellcharter") <- output
  expect_error(
    .cc_spe_require_stage(changed_output, "aggregation:X_cellcharter"),
    "`output` representation contents changed"
  )

  checkpoint <- tempfile("spe-stage-domain-")
  on.exit(unlink(checkpoint, recursive = TRUE), add = TRUE)
  clustered <- cellcharter_cluster(
    spe, candidate_k = c(2L, 3L), checkpoint_dir = checkpoint,
    max_runs = 2L, convergence_tol = 0,
    fit_fun = .tiny_fit, predict_fun = .tiny_predict,
    nll_fun = .tiny_nll, progress = FALSE
  )
  clustered$cellcharter_domain[[1L]] <- clustered$cellcharter_domain[[2L]]
  expect_error(
    cellcharter_enrichment(
      clustered, group_col = "cell_type", domain_col = "cellcharter_domain"
    ),
    "`domain` domain input changed"
  )
})

test_that("aggregation rejects cross-sample imported graph edges", {
  spe <- .tiny_spe_fixture()
  SingleCellExperiment::colPair(spe, "imported_graph") <-
    S4Vectors::SelfHits(1L, 5L, nnode = ncol(spe), value = 1)
  expect_error(
    cellcharter_aggregate(
      spe, graph_name = "imported_graph", n_layers = 1L, progress = FALSE
    ),
    "contains 1 cross-sample edge"
  )
})

test_that("graph errors identify a degenerate physical sample", {
  skip_if_not_installed("geometry")
  spe <- .tiny_spe_fixture()
  spatial <- SpatialExperiment::spatialCoords(spe)
  spatial[1:4, ] <- cbind(1:4, 1:4)
  SpatialExperiment::spatialCoords(spe) <- spatial
  expect_error(cellcharter_graph(spe), "sample-local graph for `A1`")
})

test_that("graph and aggregation APIs report progress", {
  skip_if_not_installed("geometry")
  expect_output(
    spe <- cellcharter_graph(.tiny_spe_fixture(), progress = TRUE),
    "100%",
    fixed = TRUE
  )
  expect_output(
    spe <- cellcharter_aggregate(spe, n_layers = 1L, progress = TRUE),
    "100%",
    fixed = TRUE
  )
  expect_true("X_cellcharter" %in% SingleCellExperiment::reducedDimNames(spe))
  expect_error(cellcharter_graph(.tiny_spe_fixture(), progress = NA), "progress")
  expect_error(cellcharter_aggregate(spe, progress = NA), "progress")
})

test_that("SPE-prefixed aliases forward to the canonical API", {
  skip_if_not_installed("geometry")
  expect_no_warning(result <- cellcharter_spe_graph(.tiny_spe_fixture()))
  expect_s4_class(result, "SpatialExperiment")
  expect_true("cellcharter_graph" %in% SingleCellExperiment::colPairNames(result))
})

test_that("SpatialExperiment APIs consistently expose `spe`", {
  functions <- list(
    cellcharter_graph,
    cellcharter_aggregate,
    cellcharter_cluster,
    cellcharter_enrichment,
    cellcharter_nhood_enrichment,
    cellcharter_diff_nhood_enrichment,
    cellcharter_workflow,
    cellcharter_cluster_fixed,
    cellcharter_spe_graph,
    cellcharter_spe_aggregate,
    cellcharter_spe_cluster,
    cellcharter_spe_enrichment,
    cellcharter_spe_nhood_enrichment,
    cellcharter_spe_diff_nhood_enrichment,
    cellcharter_spe_workflow
  )
  expect_true(all(vapply(
    functions,
    function(fun) identical(names(formals(fun))[[1L]], "spe"),
    logical(1L)
  )))

  skip_if_not_installed("geometry")
  expect_s4_class(
    cellcharter_graph(spe = .tiny_spe_fixture(), progress = FALSE),
    "SpatialExperiment"
  )
})
