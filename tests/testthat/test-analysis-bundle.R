.analysis_bundle_fixture <- function(source_dir) {
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

  counts <- matrix(seq_len(30), nrow = 5L,
                   dimnames = list(paste0("gene", 1:5), paste0("cell", 1:6)))
  assay <- HDF5Array::writeHDF5Array(
    counts,
    filepath = file.path(source_dir, "source-assay.h5"),
    name = "counts",
    with.dimnames = TRUE
  )
  pca_values <- matrix(seq_len(18) / 10, nrow = 6L,
                       dimnames = list(colnames(counts), paste0("PC", 1:3)))
  pca <- HDF5Array::writeHDF5Array(
    pca_values,
    filepath = file.path(source_dir, "source-pca.h5"),
    name = "PCA",
    with.dimnames = TRUE
  )
  # Keep a delayed operation in the seed tree; portability must not depend on
  # serializing an absolute HDF5ArraySeed path from the source object.
  pca <- DelayedArray::DelayedArray(pca) + 0

  spatial <- cbind(x = c(0, 1, 2, 0, 1, 2), y = c(0, 0, 0, 1, 1, 1))
  rownames(spatial) <- colnames(counts)
  col_data <- S4Vectors::DataFrame(
    sample_id = rep(c("A", "B"), each = 3L),
    cell_type = c("T", "T", "B", "B", "M", "M"),
    cellcharter = factor(c(1, 1, 2, 2, 3, 3)),
    row.names = colnames(counts)
  )
  row_data <- S4Vectors::DataFrame(
    highly_variable = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    row.names = rownames(counts)
  )
  spe <- SpatialExperiment::SpatialExperiment(
    assays = S4Vectors::SimpleList(counts = assay),
    rowData = row_data,
    colData = col_data,
    reducedDims = S4Vectors::SimpleList(PCA = pca),
    spatialCoords = spatial
  )
  graph <- S4Vectors::SelfHits(
    from = c(1L, 2L, 4L, 5L), to = c(2L, 3L, 5L, 6L),
    nnode = ncol(spe), distance = c(1, 1, 1, 1)
  )
  SingleCellExperiment::colPair(spe, "cellcharter_graph") <- graph

  cc <- list(
    stages = list(
      "graph:cellcharter_graph" = list(seed = 19L),
      "clustering:cellcharter" = list(seed = 42L)
    ),
    autok = list(cellcharter = structure(list(
      best_k = 3L,
      stability = c(`2` = 0.75, `3` = 0.91, `4` = 0.72),
      model = list(
        weights = c(0.3, 0.4, 0.3),
        means = matrix(seq_len(9) / 9, nrow = 3L),
        covariances = array(diag(3), dim = c(3L, 3L, 3L))
      )
    ), class = "cellcharter_autok")),
    identity_alignment = list(PCA = list(
      source = normalizePath(file.path(source_dir, "source-pca.h5")),
      dataset = "PCA"
    )),
    results = list(enrichment = list(value = list(
      observed = matrix(1:6, nrow = 2L,
                        dimnames = list(c("T", "B"), paste0("domain", 1:3))),
      expected = matrix(rep(1, 6), nrow = 2L,
                        dimnames = list(c("T", "B"), paste0("domain", 1:3))),
      enrichment = matrix(c(1, -1, 2, -2, 0.5, -0.5), nrow = 2L,
                          dimnames = list(c("T", "B"), paste0("domain", 1:3))),
      params = list(log2_fold_change = TRUE)
    )))
  )
  S4Vectors::metadata(spe)$CellCharterR <- cc

  list(spe = spe, counts = counts, pca = pca_values, spatial = spatial,
       graph = graph, cellcharter = cc)
}

.expect_analysis_bundle_contents <- function(observed, expected) {
  expect_s4_class(observed, "SpatialExperiment")
  expect_true(methods::validObject(observed, test = TRUE))
  expect_identical(dim(observed), c(5L, 6L))
  expect_identical(rownames(observed), rownames(expected$counts))
  expect_identical(colnames(observed), colnames(expected$counts))
  expect_equal(as.matrix(SummarizedExperiment::assay(observed, "counts")),
               expected$counts)
  expect_equal(as.matrix(SingleCellExperiment::reducedDim(observed, "PCA")),
               expected$pca)
  expect_equal(SpatialExperiment::spatialCoords(observed), expected$spatial)
  expect_identical(
    as.character(SummarizedExperiment::colData(observed)$cellcharter),
    as.character(c(1, 1, 2, 2, 3, 3))
  )
  observed_graph <- SingleCellExperiment::colPair(observed, "cellcharter_graph")
  expect_identical(S4Vectors::from(observed_graph), S4Vectors::from(expected$graph))
  expect_identical(S4Vectors::to(observed_graph), S4Vectors::to(expected$graph))
  expect_equal(S4Vectors::mcols(observed_graph)$distance,
               S4Vectors::mcols(expected$graph)$distance)
  cc <- S4Vectors::metadata(observed)$CellCharterR
  expect_equal(cc$autok$cellcharter$model$means,
               expected$cellcharter$autok$cellcharter$model$means)
  expect_equal(cc$results$enrichment$value$enrichment,
               expected$cellcharter$results$enrichment$value$enrichment)
  expect_identical(cc$identity_alignment$PCA$source, "source-pca.h5")
}

.copy_analysis_bundle <- function(from, to) {
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(from, recursive = TRUE, all.files = TRUE,
                        full.names = TRUE, no.. = TRUE, include.dirs = TRUE)
  relative <- substring(entries, nchar(from) + 2L)
  directories <- dir.exists(entries)
  if (any(directories)) {
    vapply(file.path(to, relative[directories]), dir.create, logical(1),
           recursive = TRUE, showWarnings = FALSE)
  }
  files <- !directories
  all(file.copy(entries[files], file.path(to, relative[files]),
                overwrite = TRUE, copy.mode = TRUE))
}

test_that("analysis bundle is independent, relocatable, and complete", {
  source_dir <- tempfile("analysis-source-")
  bundle <- tempfile("analysis-bundle-")
  moved_parent <- tempfile("analysis-relocated-")
  moved <- file.path(moved_parent, "renamed-analysis")
  on.exit(unlink(c(source_dir, bundle, moved_parent), recursive = TRUE,
                 force = TRUE), add = TRUE)
  fixture <- .analysis_bundle_fixture(source_dir)

  cellcharter_save_analysis(fixture$spe, bundle)
  expect_true(dir.exists(bundle))
  dir.create(moved_parent)
  expect_true(file.rename(bundle, moved))
  unlink(source_dir, recursive = TRUE, force = TRUE)

  restored <- cellcharter_load_analysis(moved)
  .expect_analysis_bundle_contents(restored, fixture)
  expect_true(DelayedArray::is_sparse(SummarizedExperiment::assay(restored)) ||
                inherits(SummarizedExperiment::assay(restored), "DelayedArray"))

  # A downstream result remains usable without graph/model recomputation.
  plot <- cellcharter_plot_enrichment(restored)
  expect_s3_class(plot, "ggplot")
  expect_true(nrow(plot$data) > 0L)
})

test_that("analysis bundle loads in a fresh R process after relocation", {
  installed_path <- tryCatch(path.package("CellCharterR"), error = function(e) "")
  if (!nzchar(installed_path) ||
      !file.exists(file.path(installed_path, "Meta", "package.rds"))) {
    skip("Fresh-process acceptance requires an installed CellCharterR package")
  }
  source_dir <- tempfile("analysis-fresh-source-")
  original <- tempfile("analysis-fresh-original-")
  relocated <- tempfile("analysis-fresh-relocated-")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(source_dir, original, relocated, script), recursive = TRUE,
                 force = TRUE), add = TRUE)
  fixture <- .analysis_bundle_fixture(source_dir)
  cellcharter_save_analysis(fixture$spe, original)
  expect_true(file.rename(original, relocated))
  unlink(source_dir, recursive = TRUE, force = TRUE)

  writeLines(c(
    "suppressPackageStartupMessages(library(CellCharterR))",
    "x <- cellcharter_load_analysis(commandArgs(TRUE)[1])",
    "stopifnot(methods::validObject(x, test = TRUE))",
    "stopifnot(identical(dim(x), c(5L, 6L)))",
    "stopifnot(identical(as.numeric(SummarizedExperiment::assay(x, 'counts')), as.numeric(matrix(1:30, 5))))",
    "stopifnot(identical(dim(SingleCellExperiment::reducedDim(x, 'PCA')), c(6L, 3L)))",
    "stopifnot(nrow(SingleCellExperiment::colPair(x, 'cellcharter_graph')) == 4L)",
    "stopifnot(S4Vectors::metadata(x)$CellCharterR$autok$cellcharter$best_k == 3L)"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(script), shQuote(relocated)),
                    stdout = TRUE, stderr = TRUE)
  expect_identical(attr(status, "status"), NULL, info = paste(status, collapse = "\n"))
})

test_that("analysis bundle relocates real AutoK checkpoint references", {
  source_dir <- tempfile("analysis-autok-source-")
  checkpoint_dir <- file.path(source_dir, "autok")
  fit_dir <- file.path(checkpoint_dir, "fits", "k_003")
  bundle <- tempfile("analysis-autok-bundle-")
  moved <- paste0(bundle, "-moved")
  on.exit(unlink(c(source_dir, bundle, moved), recursive = TRUE, force = TRUE),
          add = TRUE)
  dir.create(fit_dir, recursive = TRUE)
  fit_path <- file.path(fit_dir, "run_001.rds")
  expected_model <- list(n_clusters = 3L, means = matrix(1:9, nrow = 3L))
  saveRDS(list(model = expected_model), fit_path)

  fixture <- .analysis_bundle_fixture(source_dir)
  fit <- structure(list(
    candidate_k = 3L, stability_mean = 1, best_k = 3L,
    checkpoint_dir = normalizePath(checkpoint_dir),
    fit_index = data.frame(k = 3L, run = 1L, path = normalizePath(fit_path)),
    best_index = data.frame(k = 3L, run = 1L, path = normalizePath(fit_path)),
    selected_artifact = normalizePath(fit_path),
    state_path = normalizePath(fit_path)
  ), class = "cellcharter_autok")
  S4Vectors::metadata(fixture$spe)$CellCharterR$autok$cellcharter <- fit
  alt <- SingleCellExperiment::SingleCellExperiment(
    assays = list(score = matrix(
      seq_len(ncol(fixture$spe)), nrow = 1L,
      dimnames = list("score", colnames(fixture$spe))
    ))
  )
  S4Vectors::metadata(alt)$CellCharterR <- list(cache_path = normalizePath(fit_path))
  SingleCellExperiment::altExp(fixture$spe, "path_reuse") <- alt

  cellcharter_save_analysis(fixture$spe, bundle, progress = FALSE)
  expect_length(
    list.files(bundle, pattern = "^run_001[.]rds$", recursive = TRUE),
    1L
  )
  expect_true(file.rename(bundle, moved))
  unlink(source_dir, recursive = TRUE, force = TRUE)
  restored <- cellcharter_load_analysis(moved)
  restored_fit <- S4Vectors::metadata(restored)$CellCharterR$autok$cellcharter

  expect_true(dir.exists(restored_fit$checkpoint_dir))
  expect_true(file.exists(restored_fit$best_index$path[[1L]]))
  expect_identical(restored_fit$selected_artifact,
                   restored_fit$best_index$path[[1L]])
  expect_identical(restored_fit$state_path,
                   restored_fit$best_index$path[[1L]])
  expect_identical(
    S4Vectors::metadata(SingleCellExperiment::altExp(
      restored, "path_reuse"
    ))$CellCharterR$cache_path,
    restored_fit$best_index$path[[1L]]
  )
  expect_equal(cellcharter_autok_model(restored_fit), expected_model)
})

test_that("analysis bundle preserves rowRanges and recursive alternative experiments", {
  source_dir <- tempfile("analysis-nested-source-")
  bundle <- tempfile("analysis-nested-bundle-")
  moved <- paste0(bundle, "-moved")
  on.exit(unlink(c(source_dir, bundle, moved), recursive = TRUE, force = TRUE),
          add = TRUE)
  fixture <- .analysis_bundle_fixture(source_dir)
  ranges <- GenomicRanges::GRanges(
    seqnames = rep(c("chr1", "chr2"), length.out = nrow(fixture$spe)),
    ranges = IRanges::IRanges(start = seq_len(nrow(fixture$spe)) * 10L, width = 5L),
    strand = rep(c("+", "-"), length.out = nrow(fixture$spe))
  )
  S4Vectors::mcols(ranges)$highly_variable <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
  SummarizedExperiment::rowRanges(fixture$spe) <- ranges
  alt <- SingleCellExperiment::SingleCellExperiment(
    assays = list(score = matrix(
      seq_len(12), nrow = 2L,
      dimnames = list(c("program1", "program2"), colnames(fixture$spe))
    ))
  )
  nested <- SingleCellExperiment::SingleCellExperiment(
    assays = list(nested = matrix(
      seq_len(6), nrow = 1L,
      dimnames = list("nested1", colnames(fixture$spe))
    ))
  )
  SingleCellExperiment::altExp(alt, "nested") <- nested
  SingleCellExperiment::altExp(fixture$spe, "AUCell") <- alt
  SingleCellExperiment::mainExpName(fixture$spe) <- "protein"

  cellcharter_save_analysis(fixture$spe, bundle, progress = FALSE)
  expect_true(file.rename(bundle, moved))
  unlink(source_dir, recursive = TRUE, force = TRUE)
  restored <- cellcharter_load_analysis(moved)

  expect_identical(SingleCellExperiment::altExpNames(restored), "AUCell")
  expect_identical(
    SingleCellExperiment::altExpNames(SingleCellExperiment::altExp(restored, "AUCell")),
    "nested"
  )
  expect_equal(
    as.matrix(SummarizedExperiment::assay(
      SingleCellExperiment::altExp(restored, "AUCell"), "score"
    )),
    as.matrix(SummarizedExperiment::assay(alt, "score"))
  )
  expect_identical(SingleCellExperiment::mainExpName(restored), "protein")
  expect_identical(
    as.character(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(restored))),
    as.character(GenomicRanges::seqnames(ranges))
  )
  expect_equal(BiocGenerics::start(SummarizedExperiment::rowRanges(restored)),
               BiocGenerics::start(ranges))
})

test_that("analysis bundle validates missing and modified payloads", {
  source_dir <- tempfile("analysis-corrupt-source-")
  bundle <- tempfile("analysis-corrupt-")
  missing <- paste0(bundle, "-missing")
  modified <- paste0(bundle, "-modified")
  on.exit(unlink(c(source_dir, bundle, missing, modified), recursive = TRUE,
                 force = TRUE), add = TRUE)
  fixture <- .analysis_bundle_fixture(source_dir)
  cellcharter_save_analysis(fixture$spe, bundle)

  expect_true(.copy_analysis_bundle(bundle, missing))
  missing_payload <- list.files(missing, pattern = "\\.(h5|hdf5)$",
                                recursive = TRUE, full.names = TRUE)[1]
  expect_false(is.na(missing_payload))
  unlink(missing_payload)
  expect_error(cellcharter_load_analysis(missing),
               "missing|required|payload|manifest|checksum", ignore.case = TRUE)

  expect_true(.copy_analysis_bundle(bundle, modified))
  modified_payload <- list.files(modified, pattern = "\\.(h5|hdf5)$",
                                 recursive = TRUE, full.names = TRUE)[1]
  expect_false(is.na(modified_payload))
  writeBin(as.raw(c(1L, 2L, 3L, 4L)), modified_payload, useBytes = TRUE)
  expect_error(cellcharter_load_analysis(modified),
               "checksum|modified|corrupt|size", ignore.case = TRUE)
})

test_that("analysis bundle refuses unsafe overwrite and opaque external paths", {
  source_dir <- tempfile("analysis-path-source-")
  bundle <- tempfile("analysis-path-bundle-")
  external <- tempfile("analysis-external-")
  on.exit(unlink(c(source_dir, bundle, external), recursive = TRUE,
                 force = TRUE), add = TRUE)
  writeLines("external", external)
  fixture <- .analysis_bundle_fixture(source_dir)
  cellcharter_save_analysis(fixture$spe, bundle)
  expect_error(cellcharter_save_analysis(fixture$spe, bundle),
               "exist|overwrite", ignore.case = TRUE)

  S4Vectors::metadata(fixture$spe)$third_party <- list(cache_file = external)
  expect_error(
    cellcharter_save_analysis(fixture$spe, tempfile("unsafe-bundle-")),
    "external|absolute|portable|metadata|path", ignore.case = TRUE
  )

  for (unsafe in c("../../source.h5", "~/source.h5",
                   "\\\\server\\share\\source.h5")) {
    unsafe_source <- tempfile("analysis-unsafe-source-")
    on.exit(unlink(unsafe_source, recursive = TRUE, force = TRUE), add = TRUE)
    candidate <- .analysis_bundle_fixture(unsafe_source)$spe
    S4Vectors::metadata(candidate)$third_party <- list(cache_file = unsafe)
    destination <- tempfile("unsafe-bundle-")
    on.exit(unlink(destination, recursive = TRUE, force = TRUE), add = TRUE)
    expect_error(
      cellcharter_save_analysis(candidate, destination),
      "external|portable|metadata|path", ignore.case = TRUE,
      info = paste("unsafe path was accepted:", unsafe)
    )
  }
})
