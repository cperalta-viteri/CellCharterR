.tiny_sce_fixture <- function() {
  sample <- rep(c("A1", "A2", "B1", "B2"), each = 4L)
  condition <- rep(c("healthy", "healthy", "lupus", "lupus"), each = 4L)
  local_xy <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1.2, 1.1))
  spatial <- do.call(rbind, lapply(seq_len(4L), function(i) {
    sweep(local_xy, 2L, c(i * 10, i * 3), "+")
  }))
  ids <- paste0("cell", seq_along(sample))
  representation <- cbind(
    x = rep(c(-2, -1, 1, 2), 4L),
    y = rep(c(-1, 1, -1, 1), 4L)
  )
  rownames(spatial) <- rownames(representation) <- ids
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(protein = matrix(seq_len(32L), nrow = 2L,
                                   dimnames = list(c("m1", "m2"), ids))),
    colData = S4Vectors::DataFrame(
      sample = sample, condition = factor(condition, levels = c("healthy", "lupus")),
      cell_type = rep(c("T", "B"), 8L), row.names = ids
    ),
    reducedDims = S4Vectors::SimpleList(spatial = spatial, X_trVAE = representation)
  )
  sce
}

.tiny_spe_fixture <- function() {
  sce <- .tiny_sce_fixture()
  spatial <- SingleCellExperiment::reducedDim(sce, "spatial")
  reduced <- SingleCellExperiment::reducedDims(sce)
  reduced$spatial <- NULL
  SpatialExperiment::SpatialExperiment(
    assays = SummarizedExperiment::assays(sce),
    rowData = SummarizedExperiment::rowData(sce),
    colData = SummarizedExperiment::colData(sce),
    reducedDims = reduced,
    spatialCoords = spatial,
    sample_id = as.character(SummarizedExperiment::colData(sce)$sample)
  )
}

.tiny_fit <- function(x, n_clusters, seed, ...) list(k = n_clusters, seed = seed)
.tiny_predict <- function(model, newdata, ...) {
  (seq_len(nrow(newdata)) - 1L) %% model$k
}
.tiny_nll <- function(model) as.numeric(model$k)

test_that("colData domain labels need no CellCharterR clustering stage", {
  sce <- .tiny_sce_fixture()
  sce$python_domain <- factor(rep(c("domain_b", "domain_a"), 8L))
  enriched <- cellcharter_enrichment(
    sce, group_col = "cell_type", domain_col = "python_domain"
  )
  expect_true("enrichment" %in% names(S4Vectors::metadata(enriched)$CellCharterR$results))
})

test_that("domain labels can be read from reducedDims", {
  sce <- .tiny_sce_fixture()
  labels <- matrix(rep(0:1, 8L), ncol = 1L,
                   dimnames = list(colnames(sce), "python_domain"))
  SingleCellExperiment::reducedDim(sce, "python_labels") <- labels

  by_column <- cellcharter_enrichment(
    sce, group_col = "cell_type", domain_col = "python_domain"
  )
  by_dimension <- cellcharter_enrichment(
    sce, group_col = "cell_type", domain_col = "python_labels",
    result_name = "from_dimension"
  )
  expect_true("enrichment" %in% names(S4Vectors::metadata(by_column)$CellCharterR$results))
  expect_true("from_dimension" %in% names(S4Vectors::metadata(by_dimension)$CellCharterR$results))
})

test_that("public object API completes the tiny workflow without sample leakage", {
  skip_if_not_installed("geometry")
  checkpoint <- tempfile("sce-autok-")
  sce <- cellcharter_workflow(
    .tiny_sce_fixture(), group_col = "cell_type", candidate_k = c(2L, 3L),
    checkpoint_dir = checkpoint, n_layers = 1L, n_perms = 4L, seed = 17L,
    autok_args = list(
      max_runs = 2L, convergence_tol = 0,
      fit_fun = .tiny_fit, predict_fun = .tiny_predict, nll_fun = .tiny_nll
    )
  )

  expect_s4_class(sce, "SingleCellExperiment")
  expect_true("cellcharter_graph" %in% SingleCellExperiment::colPairNames(sce))
  expect_true("X_cellcharter" %in% SingleCellExperiment::reducedDimNames(sce))
  expect_identical(
    rownames(SingleCellExperiment::reducedDim(sce, "X_cellcharter")), colnames(sce)
  )
  expect_true("cellcharter_domain" %in% colnames(SingleCellExperiment::colData(sce)))

  graph <- SingleCellExperiment::colPair(sce, "cellcharter_graph")
  sample <- as.character(SingleCellExperiment::colData(sce)$sample)
  expect_true(all(sample[S4Vectors::from(graph)] == sample[S4Vectors::to(graph)]))
  results <- S4Vectors::metadata(sce)$CellCharterR$results
  expect_true(all(c("enrichment", "nhood_enrichment", "diff_nhood_enrichment") %in%
                    names(results)))
  expect_identical(results$diff_nhood_enrichment$value$params$library_unit, "sample")
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
  sce <- .tiny_sce_fixture()
  delayed <- DelayedArray::DelayedArray(
    SingleCellExperiment::reducedDim(sce, "X_trVAE")
  )
  SingleCellExperiment::reducedDim(sce, "X_delayed") <- delayed

  observed <- .cc_sce_reduced_dim(sce, "X_delayed")

  expect_s4_class(observed, "DelayedMatrix")
  expect_identical(DelayedArray::type(observed), "double")
  expect_equal(dim(observed), c(ncol(sce), 2L))
})

test_that("enrichment plotting accepts stored results and direct result lists", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
  result <- S4Vectors::metadata(spe)$CellCharterR$results$enrichment$value

  from_spe <- cellcharter_plot_enrichment(spe)
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

test_that("object API rejects positional identity and stale subset metadata", {
  skip_if_not_installed("geometry")
  sce <- .tiny_sce_fixture()
  expect_error(
    .cc_align_observation_ids(rev(colnames(sce)), colnames(sce)),
    "order mismatch", fixed = TRUE
  )

  sce <- cellcharter_graph(sce)
  subset <- sce[, c(1:3, 5:7, 9:11, 13:15)]
  pairs <- SingleCellExperiment::colPair(subset, "cellcharter_graph")
  expect_true(all(S4Vectors::from(pairs) <= ncol(subset)))
  expect_true(all(S4Vectors::to(pairs) <= ncol(subset)))
  expect_error(
    cellcharter_aggregate(subset),
    "stale after cell reordering or subsetting", fixed = TRUE
  )
})

test_that("graph errors identify a degenerate physical sample", {
  skip_if_not_installed("geometry")
  sce <- .tiny_sce_fixture()
  spatial <- SingleCellExperiment::reducedDim(sce, "spatial")
  spatial[1:4, ] <- cbind(1:4, 1:4)
  SingleCellExperiment::reducedDim(sce, "spatial") <- spatial
  expect_error(cellcharter_graph(sce), "sample-local graph for `A1`")
})

test_that("legacy SCE-prefixed functions warn and forward", {
  skip_if_not_installed("geometry")
  expect_warning(
    result <- cellcharter_sce_graph(.tiny_spe_fixture()),
    "deprecated"
  )
  expect_s4_class(result, "SpatialExperiment")
  expect_true("cellcharter_graph" %in% SingleCellExperiment::colPairNames(result))
})
