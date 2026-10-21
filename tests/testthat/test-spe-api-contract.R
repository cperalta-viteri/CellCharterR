test_that("public experiment APIs consistently use SPE terminology", {
  experiment_apis <- c(
    "cellcharter_graph",
    "cellcharter_aggregate",
    "cellcharter_cluster",
    "cellcharter_cluster_fixed",
    "cellcharter_enrichment",
    "cellcharter_nhood_enrichment",
    "cellcharter_diff_nhood_enrichment",
    "cellcharter_workflow",
    "cellcharter_attach_hdf5_reduced_dim",
    "cellcharter_assign_dominant_program",
    "cellcharter_zonation"
  )

  for (name in experiment_apis) {
    fun <- getExportedValue("CellCharterR", name)
    expect_identical(names(formals(fun))[[1L]], "spe", info = name)
    expect_false("sce" %in% names(formals(fun)), info = name)
  }

  exports <- getNamespaceExports("CellCharterR")
  expect_false(any(grepl("(^|_)sce($|_)", exports, ignore.case = TRUE)))
  expect_true(all(paste0("cellcharter_spe_", c(
    "graph", "aggregate", "cluster", "enrichment",
    "nhood_enrichment", "diff_nhood_enrichment", "workflow"
  )) %in% exports))
  expect_true("read_codex_spe" %in% exports)
})

test_that("named spe calls preserve the native SpatialExperiment", {
  spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
  expect_s4_class(spe, "SpatialExperiment")

  enriched <- cellcharter_enrichment(
    spe = spe,
    group_col = "cell_type_reduced",
    domain_col = "tiny_domain",
    result_name = "named_spe_enrichment",
    progress = FALSE
  )

  expect_s4_class(enriched, "SpatialExperiment")
  expect_identical(dim(enriched), dim(spe))
  expect_true("named_spe_enrichment" %in%
                names(S4Vectors::metadata(enriched)$CellCharterR$results))
})

test_that("SPE-only APIs reject a plain SingleCellExperiment", {
  plain <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(
      0, nrow = 2L, ncol = 3L,
      dimnames = list(c("g1", "g2"), paste0("spot", 1:3))
    ))
  )
  expect_error(
    cellcharter_aucell(
      spe = plain, gene_sets = list(program = "g1"), progress = FALSE
    ),
    "must be a SpatialExperiment", fixed = TRUE
  )
  skip_if_not_installed("ggplot2")
  expect_error(
    cellcharter_plot_enrichment(plain),
    "must be a SpatialExperiment", fixed = TRUE
  )
})
