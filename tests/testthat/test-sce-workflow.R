test_that("H5AD metadata decoding preserves category order and missing values", {
  skip_if_not_installed("rhdf5")
  path <- tempfile(fileext = ".h5ad")
  rhdf5::h5createFile(path)
  rhdf5::h5createGroup(path, "obs")
  rhdf5::h5write(c("c1", "c2", "c3"), path, "obs/_index")
  rhdf5::h5write(c(1.5, 2.5, 3.5), path, "obs/area")
  rhdf5::h5createGroup(path, "obs/sample")
  rhdf5::h5write(c("A", "B"), path, "obs/sample/categories")
  rhdf5::h5write(c(0L, 1L, -1L), path, "obs/sample/codes")

  value <- .cc_h5ad_dataframe(path, "obs")
  expect_s4_class(value, "DataFrame")
  expect_equal(rownames(value), c("c1", "c2", "c3"))
  expect_equal(levels(value$sample), c("A", "B"))
  expect_equal(as.character(value$sample), c("A", "B", NA))
  expect_equal(value$area, c(1.5, 2.5, 3.5))
})

test_that("file checksum is deterministic", {
  skip_if_not_installed("digest")
  path <- tempfile()
  writeBin(charToRaw("CellCharterR"), path)
  expect_identical(cellcharter_file_sha256(path), cellcharter_file_sha256(path))
  expect_match(cellcharter_file_sha256(path), "^[0-9a-f]{64}$")
})

test_that("observation identity validation rejects ambiguous identifiers", {
  expect_error(
    .cc_validate_observation_ids(NULL, "Fixture IDs"),
    "missing; positional alignment is not accepted", fixed = TRUE
  )
  expect_error(
    .cc_validate_observation_ids(c("c1", NA_character_), "Fixture IDs"),
    "missing values", fixed = TRUE
  )
  expect_error(
    .cc_validate_observation_ids(c("c1", "  "), "Fixture IDs"),
    "blank values", fixed = TRUE
  )
  expect_error(
    .cc_validate_observation_ids(c("c1", "c1"), "Fixture IDs"),
    "duplicated identifier `c1`", fixed = TRUE
  )
})

test_that("observation identity alignment distinguishes order and set mismatches", {
  expected <- c("c1", "c2", "c3")
  aligned <- .cc_align_observation_ids(expected, expected)
  expect_identical(aligned$index, seq_along(expected))
  expect_false(aligned$reordered)

  expect_error(
    .cc_align_observation_ids(c("c2", "c1", "c3"), expected),
    "order mismatch", fixed = TRUE
  )
  expect_error(
    .cc_align_observation_ids(c("c1", "c2", "wrong"), expected),
    "identifier set mismatch", fixed = TRUE
  )

  reordered <- .cc_align_observation_ids(
    c("c2", "c1", "c3"), expected, allow_reorder = TRUE
  )
  expect_identical(reordered$index, c(2L, 1L, 3L))
  expect_true(reordered$reordered)
})

test_that("sample-scoped CODEX identities allow only unambiguous local IDs", {
  skip_if_not_installed("rhdf5")
  raw <- tempfile(fileext = ".h5ad")
  clustered <- tempfile(fileext = ".h5ad")
  for (path in c(raw, clustered)) {
    rhdf5::h5createFile(path)
    rhdf5::h5createGroup(path, "obs")
  }
  rhdf5::h5write(c("0", "1", "0-1", "1-1"), raw, "obs/_index")
  rhdf5::h5write(c("0", "1", "0", "1"), clustered, "obs/_index")
  raw_col_data <- S4Vectors::DataFrame(
    sample = c("A", "A", "B", "B"),
    row.names = c("0", "1", "0-1", "1-1")
  )
  clustered_col_data <- S4Vectors::DataFrame(
    sample = c("A", "A", "B", "B"),
    row.names = c("0", "1", "0", "1")
  )

  identity <- .cc_codex_h5ad_identity(
    raw, clustered, raw_col_data, clustered_col_data
  )
  expect_identical(
    identity$alignment, "exact sample-scoped local observation identifiers"
  )
  expect_false(anyDuplicated(identity$canonical_clustered_ids) > 0L)

  clustered_col_data$sample[[4L]] <- "A"
  expect_error(
    .cc_codex_h5ad_identity(raw, clustered, raw_col_data, clustered_col_data),
    "duplicated identifier", fixed = TRUE
  )
})

test_that("HDF5 representation attachment validates and records cell order", {
  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("SingleCellExperiment")
  path <- tempfile(fileext = ".h5")
  source <- matrix(c(20, 10, 30, 2, 1, 3), ncol = 2L)
  HDF5Array::writeHDF5Array(source, filepath = path, name = "embedding")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 1L, ncol = 3L))
  )
  colnames(sce) <- c("c1", "c2", "c3")

  expect_error(
    cellcharter_attach_hdf5_reduced_dim(sce, path, "embedding"),
    "positional alignment is not accepted", fixed = TRUE
  )
  expect_error(
    cellcharter_attach_hdf5_reduced_dim(
      sce, path, "embedding", cell_ids = c("c2", "c1", "c3")
    ),
    "order mismatch", fixed = TRUE
  )

  attached <- cellcharter_attach_hdf5_reduced_dim(
    sce, path, "embedding", cell_ids = c("c2", "c1", "c3"),
    allow_reorder = TRUE
  )
  expect_identical(
    rownames(SingleCellExperiment::reducedDim(attached, "embedding")),
    colnames(attached)
  )
  expect_equal(
    unname(as.matrix(
      SingleCellExperiment::reducedDim(attached, "embedding")
    )[, 1L]),
    c(10, 20, 30)
  )
  expect_true(
    S4Vectors::metadata(attached)$CellCharterR$identity_alignment$embedding$reordered
  )
})

test_that("packaged tiny H5AD pair reconstructs the portable tiny SPE", {
  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("SpatialExperiment")
  raw <- system.file(
    "extdata", "codex_mouse_spleen_tiny_raw.h5ad",
    package = "CellCharterR", mustWork = TRUE
  )
  clustered <- system.file(
    "extdata", "codex_mouse_spleen_tiny_clustered.h5ad",
    package = "CellCharterR", mustWork = TRUE
  )
  expected <- cellcharter_mouse_spleen_tiny()

  observed <- read_codex_spe(raw, clustered)
  expect_s4_class(observed, "SpatialExperiment")
  expect_s4_class(expected, "SpatialExperiment")
  observed_inputs <- S4Vectors::metadata(observed)$CellCharterR$inputs
  expect_identical(colnames(observed), colnames(expected))
  expect_identical(rownames(observed), rownames(expected))
  expect_identical(observed_inputs$representation_dataset, "obsm/X_trVAE")
  expect_identical(observed_inputs$spatial_dataset, "obsm/spatial")
  expect_identical(
    observed_inputs$raw_sha256,
    cellcharter_file_sha256(raw)
  )
  expect_identical(
    observed_inputs$clustered_sha256,
    cellcharter_file_sha256(clustered)
  )
  expect_identical(
    as.character(observed$reference_cluster_cellcharter),
    as.character(expected$reference_cluster_cellcharter)
  )
  expect_equal(
    as.matrix(SummarizedExperiment::assay(observed, "protein")),
    as.matrix(SummarizedExperiment::assay(expected, "protein")),
    tolerance = 0
  )
  expect_equal(
    as.matrix(SingleCellExperiment::reducedDim(observed, "X_trVAE")),
    as.matrix(SingleCellExperiment::reducedDim(expected, "X_trVAE")),
    tolerance = 0
  )
  expect_false("spatial" %in% SingleCellExperiment::reducedDimNames(observed))
  expect_equal(
    SpatialExperiment::spatialCoords(observed),
    SpatialExperiment::spatialCoords(expected), tolerance = 0
  )
  index <- seq(1L, ncol(observed), by = 7L)
  sampled <- observed[, index]
  expect_s4_class(sampled, "SpatialExperiment")
  expect_identical(colnames(sampled), colnames(observed)[index])
  expect_equal(
    SpatialExperiment::spatialCoords(sampled),
    SpatialExperiment::spatialCoords(observed)[index, , drop = FALSE], tolerance = 0
  )
})

test_that("version-neutral teaching workflow restores graph and results", {
  spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
  expect_s4_class(spe, "SpatialExperiment")
  expect_true(methods::validObject(spe, test = TRUE))
  expect_identical(dim(spe), c(29L, 1000L))
  expect_true("i-niche" %in% colnames(SummarizedExperiment::colData(spe)))
  expect_false("i.niche" %in% colnames(SummarizedExperiment::colData(spe)))
  expect_true("cellcharter_graph" %in% SingleCellExperiment::colPairNames(spe))
  expect_true(all(c("X_trVAE", "X_cellcharter", "X_cellcharter_tiny") %in%
                    SingleCellExperiment::reducedDimNames(spe)))
  expect_true("tiny_domain" %in% colnames(SummarizedExperiment::colData(spe)))
  pair <- SingleCellExperiment::colPair(spe, "cellcharter_graph")
  sample_id <- as.character(SummarizedExperiment::colData(spe)$sample)
  expect_true(all(sample_id[S4Vectors::from(pair)] == sample_id[S4Vectors::to(pair)]))
  info <- S4Vectors::metadata(spe)$CellCharterR
  expect_true(all(c("enrichment", "nhood_enrichment", "diff_nhood_enrichment") %in%
                    names(info$results)))
  expect_identical(info$autok$tiny_domain$candidate_k, 2:10)
  expect_identical(info$autok$tiny_domain$max_runs, 3L)
  expect_error(cellcharter_mouse_spleen_tiny(with_results = NA), "TRUE or FALSE")
})
