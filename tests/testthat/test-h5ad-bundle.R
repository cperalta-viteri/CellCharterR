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
  rhdf5::h5write(c(0L, 0L, 0L, 1L, 1L), source_file, "obs/slide_id_concat/codes")
  requested_obs <- c(
    "x_array", "y_array", "x_pixel", "y_pixel", "x_pixel_fullres",
    "y_pixel_fullres", "CXCL12_iCAF", "CXCL8_iCAF", "Conventional_DC",
    "Cytotoxic_lymphocytes", "Dysfunctional_TC", "ECM-CAF",
    "IFN_Macrophage", "IL1B_Macrophage", "Neutrophils",
    "SPP1_Macrophage", "TLS_humoral", "Treg", "epithelial_cell",
    "cellcharter", "risk_group"
  )
  for (i in seq_along(requested_obs)) {
    rhdf5::h5write(seq_len(5) + i / 100, source_file,
                   paste0("obs/", requested_obs[[i]]))
  }
  rhdf5::h5write(paste0("gene", 1:3), source_file, "var/_index")
  requested_var <- c("variance", "mean", "variance_trend",
                     "residual_variance", "highly_variable")
  for (i in seq_along(requested_var)) {
    rhdf5::h5write(seq_len(3) + i / 100, source_file,
                   paste0("var/", requested_var[[i]]))
  }
  rhdf5::h5write(rbind(1:5, 11:15), source_file, "obsm/spatial")
  rhdf5::h5write(matrix(seq_len(10), nrow = 2), source_file, "obsm/X_pca")
  rhdf5::h5write(matrix(seq_len(15), nrow = 3), source_file, "obsm/X_cellcharter")
  rhdf5::h5write(seq_len(5) / 5, source_file, "obsm/score_aucell/program")
  rhdf5::h5createGroup(source_file, "obsp/spatial_connectivities")
  rhdf5::h5write(c(0.5, 0.5, 1, 1, 0.7, 0.7), source_file,
                 "obsp/spatial_connectivities/data")
  rhdf5::h5write(c(1L, 0L, 2L, 1L, 4L, 3L), source_file,
                 "obsp/spatial_connectivities/indices")
  rhdf5::h5write(c(0L, 1L, 3L, 4L, 5L, 6L), source_file,
                 "obsp/spatial_connectivities/indptr")
  rhdf5::h5write(c(5L, 5L), source_file,
                 "obsp/spatial_connectivities/shape")
  rhdf5::h5createGroup(source_file, "obsp/spatial_distances")
  rhdf5::h5write(c(10, 10, 20, 20, 30, 30), source_file,
                 "obsp/spatial_distances/data")
  rhdf5::h5write(c(1L, 0L, 2L, 1L, 4L, 3L), source_file,
                 "obsp/spatial_distances/indices")
  rhdf5::h5write(c(0L, 1L, 3L, 4L, 5L, 6L), source_file,
                 "obsp/spatial_distances/indptr")
  rhdf5::h5write(c(5L, 5L), source_file, "obsp/spatial_distances/shape")
  expect_error(
    .transition_read_obsp(source_file, "/obsp/spatial_connectivities", 5L),
    "lacks an encoding-type"
  )
  file_handle <- rhdf5::H5Fopen(source_file)
  obs_handle <- rhdf5::H5Gopen(file_handle, "obs")
  rhdf5::h5writeAttribute(c("slide_id_concat", requested_obs), obs_handle,
                          "column-order")
  rhdf5::H5Gclose(obs_handle)
  var_handle <- rhdf5::H5Gopen(file_handle, "var")
  rhdf5::h5writeAttribute(requested_var, var_handle, "column-order")
  rhdf5::H5Gclose(var_handle)
  for (graph_name in c("spatial_connectivities", "spatial_distances")) {
    group_handle <- rhdf5::H5Gopen(file_handle, paste0("obsp/", graph_name))
    rhdf5::h5writeAttribute("csr_matrix", group_handle, "encoding-type")
    rhdf5::h5writeAttribute("0.1.0", group_handle, "encoding-version")
    rhdf5::H5Gclose(group_handle)
  }
  rhdf5::H5Fclose(file_handle)
  rhdf5::h5createGroup(source_file, "uns/slides")
  rhdf5::h5createGroup(source_file, "uns/spatial_neighbors")
  rhdf5::h5write("slide-a", source_file, "uns/slides/name")
  rhdf5::h5write("delaunay", source_file, "uns/spatial_neighbors/method")

  spe <- cellcharter_import_h5ad(source_file, bundle)
  expect_s4_class(spe, "SpatialExperiment")
  expect_identical(dim(spe), c(3L, 5L))
  expect_true(file.exists(file.path(bundle, "arrays.h5")))
  expect_true(file.exists(file.path(bundle, "spe.rds")))
  expect_setequal(SingleCellExperiment::colPairNames(spe),
                  c("spatial_connectivities", "spatial_distances"))
  expect_identical(colnames(SummarizedExperiment::colData(spe)),
                   c("slide_id_concat", requested_obs, "sample_id"))
  expect_identical(colnames(SummarizedExperiment::rowData(spe)), requested_var)
  expect_equal(unname(SpatialExperiment::spatialCoords(spe)),
               unname(cbind(x = 1:5, y = 11:15)))
  expect_equal(unname(as.matrix(SingleCellExperiment::reducedDim(spe, "PCA"))),
               unname(t(matrix(seq_len(10), nrow = 2))))
  expect_equal(unname(as.matrix(SingleCellExperiment::reducedDim(spe, "cellcharter"))),
               unname(t(matrix(seq_len(15), nrow = 3))))
  expect_equal(as.numeric(SingleCellExperiment::reducedDim(spe, "AUCell")[, 1]),
               seq_len(5) / 5)
  expected_graph <- matrix(0, 5, 5)
  expected_graph[cbind(c(1, 2, 2, 3, 4, 5), c(2, 1, 3, 2, 5, 4))] <-
    c(0.5, 0.5, 1, 1, 0.7, 0.7)
  expect_equal(as.matrix(SingleCellExperiment::colPair(
    spe, "spatial_connectivities", asSparse = TRUE
  )), expected_graph)
  imported_pair <- SingleCellExperiment::colPair(
    spe, "spatial_connectivities", asSparse = FALSE
  )
  expect_identical(S4Vectors::from(imported_pair),
                   c(1L, 2L, 2L, 3L, 4L, 5L))
  expect_identical(S4Vectors::to(imported_pair),
                   c(2L, 1L, 3L, 2L, 5L, 4L))
  expect_identical(as.numeric(S4Vectors::mcols(imported_pair)$value),
                   c(0.5, 0.5, 1, 1, 0.7, 0.7))
  expect_identical(
    S4Vectors::metadata(spe)$transitionR$obsp_report$spatial_connectivities$nnz,
    6L
  )
  expect_setequal(S4Vectors::metadata(spe)$anndata_uns$entries,
                  c("slides", "spatial_neighbors"))
  expect_error(
    .transition_build_spe(
      file.path(bundle, "arrays.h5"), "X", "expression", "x_array", NULL,
      "obsm/spatial", character(), NULL, "slide_id_concat"
    ),
    "sample_column must be imported"
  )

  unlink(source_file)
  expect_true(file.rename(bundle, moved))
  reloaded <- cellcharter_load_spe(moved)
  expect_identical(as.numeric(SummarizedExperiment::assay(reloaded)),
                   as.numeric(matrix(seq_len(15), nrow = 3)))
  expect_identical(dim(SingleCellExperiment::reducedDim(reloaded, "PCA")), c(5L, 2L))
  expect_identical(as.numeric(SummarizedExperiment::colData(reloaded)$risk_group),
                   seq_len(5) + match("risk_group", requested_obs) / 100)
  expect_equal(as.matrix(SingleCellExperiment::colPair(
    reloaded, "spatial_connectivities", asSparse = TRUE
  )), expected_graph)
  expected_distances <- matrix(0, 5, 5)
  expected_distances[cbind(c(1, 2, 2, 3, 4, 5),
                           c(2, 1, 3, 2, 5, 4))] <-
    c(10, 10, 20, 20, 30, 30)
  expect_equal(as.matrix(SingleCellExperiment::colPair(
    reloaded, "spatial_distances", asSparse = TRUE
  )), expected_distances)
  expect_identical(
    normalizePath(S4Vectors::metadata(reloaded)$anndata_uns$filepath),
    normalizePath(file.path(moved, "arrays.h5"))
  )
  expect_equal(as.matrix(.cc_sce_graph_matrix(
    reloaded, "spatial_connectivities"
  ) != 0), expected_graph != 0)
  expect_true(methods::validObject(reloaded, test = TRUE))
})

test_that("imported pair values are separated from graph topology semantics", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(x = matrix(0, nrow = 1L, ncol = 3L))
  )
  SingleCellExperiment::colPair(sce, "graph") <- S4Vectors::SelfHits(
    from = c(1L, 2L), to = c(2L, 3L), nnode = 3L,
    value = c(0, 2), sort.by.query = FALSE
  )
  topology <- .cc_sce_graph_matrix(sce, "graph")
  expect_identical(as.numeric(topology[1, 2]), 0)
  expect_identical(as.numeric(topology[2, 3]), 1)

  SingleCellExperiment::colPair(sce, "negative") <- S4Vectors::SelfHits(
    from = 1L, to = 2L, nnode = 3L,
    value = -1, sort.by.query = FALSE
  )
  expect_error(.cc_sce_graph_matrix(sce, "negative"), "negative edge weights")
})
