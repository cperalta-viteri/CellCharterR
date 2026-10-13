test_that("H5AD import creates an independent movable SPE bundle", {
  source_file <- tempfile(fileext = ".h5ad")
  bundle <- tempfile(pattern = "cellcharterr-spe-")
  moved <- paste0(bundle, "-moved")
  on.exit(unlink(c(source_file, bundle, moved), recursive = TRUE), add = TRUE)

  rhdf5::h5createFile(source_file)
  for (group in c("obs", "var", "obsm", "obsm/score_aucell", "obsp", "uns")) {
    rhdf5::h5createGroup(source_file, group)
  }
  rhdf5::h5write(matrix(seq_len(15), nrow = 3), source_file, "X")
  rhdf5::h5write(paste0("cell", 1:5), source_file, "obs/_index")
  rhdf5::h5createGroup(source_file, "obs/slide_id_concat")
  rhdf5::h5write(c("A", "B"), source_file, "obs/slide_id_concat/categories")
  rhdf5::h5write(c(0L, 0L, 1L, 1L, 1L), source_file, "obs/slide_id_concat/codes")
  rhdf5::h5write(seq_len(5) / 10, source_file, "obs/score")
  rhdf5::h5write(paste0("gene", 1:3), source_file, "var/_index")
  rhdf5::h5write(c(TRUE, FALSE, TRUE), source_file, "var/highly_variable")
  rhdf5::h5write(rbind(1:5, 11:15), source_file, "obsm/spatial")
  rhdf5::h5write(matrix(seq_len(10), nrow = 2), source_file, "obsm/X_pca")
  rhdf5::h5write(matrix(seq_len(15), nrow = 3), source_file, "obsm/X_cellcharter")
  rhdf5::h5write(seq_len(5) / 5, source_file, "obsm/score_aucell/program")
  rhdf5::h5write("preserved", source_file, "uns/provenance")

  spe <- cellcharter_import_h5ad(source_file, bundle)
  expect_s4_class(spe, "SpatialExperiment")
  expect_identical(dim(spe), c(3L, 5L))
  expect_true(file.exists(file.path(bundle, "arrays.h5")))
  expect_true(file.exists(file.path(bundle, "spe.rds")))
  expect_identical(as.character(rhdf5::h5read(
    file.path(bundle, "arrays.h5"), "uns/provenance"
  )), "preserved")
  expect_error(
    .transition_build_spe(
      file.path(bundle, "arrays.h5"), "X", "expression", "score", NULL,
      "obsm/spatial", character(), NULL, "slide_id_concat"
    ),
    "sample_column must be imported"
  )

  unlink(source_file)
  expect_true(file.rename(bundle, moved))
  reloaded <- cellcharter_load_spe_bundle(moved)
  expect_identical(as.numeric(SummarizedExperiment::assay(reloaded)),
                   as.numeric(matrix(seq_len(15), nrow = 3)))
  expect_identical(dim(SingleCellExperiment::reducedDim(reloaded, "PCA")), c(5L, 2L))
  expect_identical(as.numeric(SummarizedExperiment::colData(reloaded)$score),
                   seq_len(5) / 10)
  expect_true(methods::validObject(reloaded, test = TRUE))
})
