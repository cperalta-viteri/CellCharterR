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

test_that("exported SCE adapters complete the tiny workflow without sample leakage", {
  skip_if_not_installed("geometry")
  checkpoint <- tempfile("sce-autok-")
  sce <- cellcharter_sce_workflow(
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
  result <- cellcharter_sce_workflow(
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

test_that("SCE adapters reject positional identity and stale subset metadata", {
  skip_if_not_installed("geometry")
  sce <- .tiny_sce_fixture()
  expect_error(
    .cc_align_observation_ids(rev(colnames(sce)), colnames(sce)),
    "order mismatch", fixed = TRUE
  )

  sce <- cellcharter_sce_graph(sce)
  subset <- sce[, c(1:3, 5:7, 9:11, 13:15)]
  pairs <- SingleCellExperiment::colPair(subset, "cellcharter_graph")
  expect_true(all(S4Vectors::from(pairs) <= ncol(subset)))
  expect_true(all(S4Vectors::to(pairs) <= ncol(subset)))
  expect_error(
    cellcharter_sce_aggregate(subset),
    "stale after cell reordering or subsetting", fixed = TRUE
  )
})

test_that("SCE graph errors identify a degenerate physical sample", {
  skip_if_not_installed("geometry")
  sce <- .tiny_sce_fixture()
  spatial <- SingleCellExperiment::reducedDim(sce, "spatial")
  spatial[1:4, ] <- cbind(1:4, 1:4)
  SingleCellExperiment::reducedDim(sce, "spatial") <- spatial
  expect_error(cellcharter_sce_graph(sce), "sample-local graph for `A1`")
})
