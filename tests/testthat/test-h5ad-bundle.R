expect_canonical_import_storage <- function(spe) {
  assay <- SummarizedExperiment::assay(spe, "expression")
  expect_s4_class(assay, "DelayedMatrix")
  expect_s4_class(DelayedArray::seed(assay), "HDF5ArraySeed")

  col_data <- SummarizedExperiment::colData(spe)
  row_data <- SummarizedExperiment::rowData(spe)
  for (value in as.list(col_data)) {
    expect_false(inherits(value, "DelayedArray"))
    expect_true(is.atomic(value) || is.factor(value))
  }
  for (value in as.list(row_data)) {
    expect_false(inherits(value, "DelayedArray"))
    expect_true(is.atomic(value) || is.factor(value))
  }
  expect_true(is.matrix(SpatialExperiment::spatialCoords(spe)))
  expect_type(SpatialExperiment::spatialCoords(spe), "integer")
  for (name in SingleCellExperiment::reducedDimNames(spe)) {
    value <- SingleCellExperiment::reducedDim(spe, name)
    expect_true(is.matrix(value))
    expect_true(typeof(value) %in% c("integer", "double"))
    expect_false(inherits(value, "DelayedArray"))
  }
  for (name in SingleCellExperiment::colPairNames(spe)) {
    expect_s4_class(SingleCellExperiment::colPair(spe, name), "SelfHits")
  }
  expect_type(S4Vectors::metadata(spe)$anndata_uns, "list")
  expect_identical(
    S4Vectors::metadata(spe)$transitionR$canonical_storage$assay,
    "HDF5-backed"
  )
}

test_that("H5AD importer exposes its complete conversion contract", {
  expect_identical(
    names(formals(cellcharter_import_h5ad)),
    c(
      "h5ad_file", "output_dir", "assay_path", "assay_name", "obs_columns",
      "var_columns", "spatial_path", "reduced_dims", "aucell_path",
      "sample_column", "chunk_mib", "overwrite", "source_verification",
      "materialize_block_mib", "max_materialized_gib", "progress"
    )
  )
})

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
  rhdf5::h5write(paste0("cell", 1:5), source_file,
                 "obsm/score_aucell/_index")
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

  spe <- cellcharter_import_h5ad(
    source_file, bundle, materialize_block_mib = 1e-6, progress = FALSE
  )
  expect_s4_class(spe, "SpatialExperiment")
  expect_canonical_import_storage(spe)
  expect_identical(dim(spe), c(3L, 5L))
  expect_true(file.exists(file.path(bundle, "arrays.h5")))
  expect_true(file.exists(file.path(bundle, "spe.rds")))
  manifest <- readRDS(file.path(bundle, "manifest.rds"))
  expect_match(manifest$arrays_sha256, "^[0-9a-f]{64}$")
  expect_identical(manifest$arrays_size,
                   unname(file.info(file.path(bundle, "arrays.h5"))$size))
  expect_identical(manifest$arrays_fingerprint,
                   .transition_file_fingerprint(file.path(bundle, "arrays.h5")))
  expect_setequal(SingleCellExperiment::colPairNames(spe),
                  c("spatial_connectivities", "spatial_distances"))
  expect_identical(colnames(SummarizedExperiment::colData(spe)),
                   c("slide_id_concat", requested_obs, "sample_id"))
  expect_identical(colnames(SummarizedExperiment::rowData(spe)), requested_var)
  expect_identical(SummarizedExperiment::colData(spe)$slide_id_concat,
                   factor(c("A", "A", "A", "B", "B"), levels = c("A", "B")))
  expect_identical(
    as.numeric(SummarizedExperiment::rowData(spe)$variance),
    seq_len(3) + 0.01
  )
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
  expect_gt(
    S4Vectors::metadata(spe)$transitionR$materialization_preflight$retained_bytes,
    0
  )
  expect_setequal(S4Vectors::metadata(spe)$anndata_uns$entries,
                  c("slides", "spatial_neighbors"))

  # A rejected overwrite is read-only, including the reconstruction manifest.
  manifest_path <- file.path(bundle, "manifest.rds")
  manifest_before <- digest::digest(file = manifest_path, algo = "sha256")
  spe_before <- digest::digest(file = file.path(bundle, "spe.rds"),
                               algo = "sha256")
  entries_before <- list.files(bundle, all.files = TRUE, no.. = TRUE,
                               recursive = TRUE, include.dirs = TRUE)
  directory_mtime_before <- file.info(bundle)$mtime
  expect_error(
    cellcharter_import_h5ad(
      source_file, bundle, assay_name = "changed_name",
      overwrite = FALSE, progress = FALSE
    ),
    "spe.rds exists"
  )
  expect_identical(digest::digest(file = manifest_path, algo = "sha256"),
                   manifest_before)
  expect_identical(digest::digest(file = file.path(bundle, "spe.rds"),
                                  algo = "sha256"), spe_before)
  expect_identical(
    list.files(bundle, all.files = TRUE, no.. = TRUE,
               recursive = TRUE, include.dirs = TRUE),
    entries_before
  )
  expect_identical(file.info(bundle)$mtime, directory_mtime_before)
  expect_false(any(grepl("[.]conversion-lock", list.files(
    bundle, all.files = TRUE, no.. = TRUE
  ))))
  expect_error(
    .transition_build_spe(
      file.path(bundle, "arrays.h5"), "X", "expression", "x_array", NULL,
      "obsm/spatial", character(), NULL, "slide_id_concat"
    ),
    "sample_column must be imported"
  )
  expect_error(
    .transition_build_spe(
      file.path(bundle, "arrays.h5"), "X", "expression", NULL, NULL,
      "obsm/spatial",
      c(PCA = "obsm/X_pca", cellcharter = "obsm/X_cellcharter"),
      "obsm/score_aucell", "slide_id_concat",
      max_materialized_gib = 1e-12
    ),
    "Estimated native peak materialization.*exceeds"
  )

  # Only the historical format-1 manifest may omit its reconstruction recipe.
  # A newer manifest without import_spec is incomplete/corrupt and must not be
  # silently rebuilt using unrelated defaults.
  current_manifest <- readRDS(file.path(bundle, "manifest.rds"))
  for (version in 2:3) {
    incomplete_manifest <- current_manifest
    incomplete_manifest$format_version <- as.integer(version)
    incomplete_manifest$import_spec <- NULL
    saveRDS(incomplete_manifest, file.path(bundle, "manifest.rds"),
            compress = FALSE)
    expect_error(
      cellcharter_load_spe(bundle, verify = "none", progress = FALSE),
      paste0("format ", version, " manifest lacks the required `import_spec`")
    )
  }
  saveRDS(current_manifest, file.path(bundle, "manifest.rds"),
          compress = FALSE)

  # The first public converter wrote format 1 without its reconstruction
  # recipe. Historical defaults remain automatic, while a custom legacy
  # recipe can be supplied explicitly instead of being guessed.
  legacy_manifest <- readRDS(file.path(bundle, "manifest.rds"))
  legacy_manifest$format_version <- 1L
  legacy_manifest$import_spec <- NULL
  saveRDS(legacy_manifest, file.path(bundle, "manifest.rds"), compress = FALSE)
  expect_message(
    legacy_default <- cellcharter_load_spe(
      bundle, verify = "none", materialize_block_mib = 1e-6,
      progress = FALSE
    ),
    "Legacy format-1 import"
  )
  expect_canonical_import_storage(legacy_default)
  legacy_spec <- list(
    assay_path = "X", assay_name = "legacy_expression",
    obs_columns = NULL, var_columns = NULL,
    spatial_path = "obsm/spatial",
    reduced_dims = c(PCA = "obsm/X_pca",
                     cellcharter = "obsm/X_cellcharter"),
    aucell_path = "obsm/score_aucell", sample_column = "slide_id_concat"
  )
  legacy_custom <- cellcharter_load_spe(
    bundle, verify = "none", legacy_import_spec = legacy_spec,
    materialize_block_mib = 1e-6, progress = FALSE
  )
  expect_identical(SummarizedExperiment::assayNames(legacy_custom),
                   "legacy_expression")

  # Loading a directory produced by the first public converter must be a pure
  # Bioconductor/HDF5 operation. In particular, the optional Torch namespace
  # must not be initialized (or attempt to download Lantern) as a side effect
  # of loading CellCharterR or reconstructing the legacy SPE.
  legacy_without_hash <- legacy_manifest
  legacy_without_hash$arrays_sha256 <- NULL
  legacy_without_hash$arrays_fingerprint <- NULL
  saveRDS(legacy_without_hash, file.path(bundle, "manifest.rds"),
          compress = FALSE)
  child_script <- sprintf(
    paste0(
      "library(CellCharterR); ",
      "stopifnot(!('torch' %%in%% loadedNamespaces())); ",
      "x <- suppressWarnings(cellcharter_load_spe(%s, verify='none', progress=FALSE)); ",
      "stopifnot(methods::validObject(x), identical(dim(x), c(3L, 5L)), ",
      "!('torch' %%in%% loadedNamespaces()))"
    ),
    encodeString(bundle, quote = "\"")
  )
  child_output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "-e", shQuote(child_script)),
    stdout = TRUE, stderr = TRUE
  )
  expect_identical(attr(child_output, "status"), NULL,
                   info = paste(child_output, collapse = "\n"))
  saveRDS(legacy_manifest, file.path(bundle, "manifest.rds"),
          compress = FALSE)

  future_manifest <- legacy_manifest
  future_manifest$format_version <- 99L
  saveRDS(future_manifest, file.path(bundle, "manifest.rds"), compress = FALSE)
  expect_error(
    cellcharter_load_spe(bundle, verify = "none", progress = FALSE),
    "Unsupported SPE import-directory format version"
  )
  saveRDS(legacy_manifest, file.path(bundle, "manifest.rds"), compress = FALSE)

  unlink(source_file)
  expect_true(file.rename(bundle, moved))
  reloaded <- cellcharter_load_spe(
    moved, materialize_block_mib = 1e-6, progress = FALSE
  )
  expect_canonical_import_storage(reloaded)
  expect_identical(as.numeric(SummarizedExperiment::assay(reloaded)),
                   as.numeric(matrix(seq_len(15), nrow = 3)))
  expect_identical(dim(SingleCellExperiment::reducedDim(reloaded, "PCA")), c(5L, 2L))
  expect_identical(as.numeric(SummarizedExperiment::colData(reloaded)$risk_group),
                   seq_len(5) + match("risk_group", requested_obs) / 100)
  expect_identical(
    as.numeric(SummarizedExperiment::rowData(reloaded)$variance),
    seq_len(3) + 0.01
  )
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
  expect_equal(as.matrix(.cc_spe_graph_matrix(
    reloaded, "spatial_connectivities"
  ) != 0), expected_graph != 0)
  expect_true(methods::validObject(reloaded, test = TRUE))
})

test_that("bundle SHA-256 detects same-size corruption without realizing assays", {
  artifact <- tempfile(fileext = ".h5")
  on.exit(unlink(artifact), add = TRUE)
  writeBin(as.raw(seq_len(128L) %% 256L), artifact)
  manifest <- list(
    arrays_size = unname(file.info(artifact)$size),
    arrays_sha256 = .transition_file_sha256(artifact)
  )
  expect_true(.transition_verify_arrays(artifact, manifest, "full"))

  connection <- file(artifact, "r+b")
  seek(connection, 64L)
  writeBin(as.raw(0L), connection)
  close(connection)
  expect_identical(unname(file.info(artifact)$size), manifest$arrays_size)
  expect_error(
    .transition_verify_arrays(artifact, manifest, "full"),
    "SHA-256 differs"
  )
})

test_that("auto verification avoids hashing an unchanged published artifact", {
  artifact <- tempfile(fileext = ".h5")
  on.exit(unlink(artifact), add = TRUE)
  writeBin(as.raw(seq_len(64L)), artifact)
  manifest <- list(
    arrays_size = unname(file.info(artifact)$size),
    arrays_sha256 = .transition_file_sha256(artifact),
    arrays_fingerprint = .transition_file_fingerprint(artifact)
  )
  observed <- .transition_verify_arrays(
    artifact, manifest, "auto",
    hash_fun = function(...) stop("unexpected full hash")
  )
  expect_true(observed)
})

test_that("auto verification hashes after a same-size artifact changes", {
  artifact <- tempfile(fileext = ".h5")
  on.exit(unlink(artifact), add = TRUE)
  writeBin(as.raw(seq_len(128L) %% 256L), artifact)
  manifest <- list(
    arrays_size = unname(file.info(artifact)$size),
    arrays_sha256 = .transition_file_sha256(artifact),
    arrays_fingerprint = .transition_file_fingerprint(artifact)
  )

  connection <- file(artifact, "r+b")
  seek(connection, 64L)
  writeBin(as.raw(0L), connection)
  close(connection)
  Sys.setFileTime(artifact, as.POSIXct(manifest$arrays_fingerprint$mtime,
                                      origin = "1970-01-01", tz = "UTC") + 10)
  expect_identical(unname(file.info(artifact)$size), manifest$arrays_size)
  expect_error(
    .transition_verify_arrays(artifact, manifest, "auto"),
    "SHA-256 differs"
  )
})

test_that("full source verification detects unsampled same-size changes", {
  source_file <- tempfile(fileext = ".h5ad")
  bundle <- tempfile("full-source-verification-")
  on.exit(unlink(c(source_file, bundle), recursive = TRUE, force = TRUE),
          add = TRUE)
  rhdf5::h5createFile(source_file)
  for (group in c("obs", "var", "obsm")) {
    rhdf5::h5createGroup(source_file, group)
  }
  rhdf5::h5write(matrix(seq_len(4), nrow = 2L), source_file, "X")
  # HDF5 permits trailing bytes. Incompressible padding makes room for a
  # mutation that is deliberately outside the three sampled identity windows.
  set.seed(11)
  padding <- as.raw(sample.int(256L, 200000L, replace = TRUE) - 1L)
  append_connection <- file(source_file, "ab")
  writeBin(padding, append_connection)
  close(append_connection)
  dir.create(bundle)

  expect_message(
    clone <- .transition_clone_file(
      source_file, bundle, chunk_mib = 1,
      source_verification = "full"
    ),
    "Computing source H5AD SHA-256"
  )
  expect_true(file.exists(clone))
  original_identity <- .transition_source_identity(source_file)
  mutation_offset <- floor(original_identity$size / 4)

  connection <- file(source_file, "r+b")
  seek(connection, where = mutation_offset, origin = "start", rw = "read")
  old <- readBin(connection, "raw", n = 1L)
  expect_length(old, 1L)
  seek(connection, where = mutation_offset, origin = "start", rw = "write")
  writeBin(as.raw((as.integer(old) + 1L) %% 256L), connection)
  close(connection)
  Sys.setFileTime(
    source_file,
    as.POSIXct(original_identity$mtime, origin = "1970-01-01", tz = "UTC")
  )
  # Some filesystems round sub-second mtimes. Mirror that rounded timestamp in
  # the legacy manifest to construct the precise collision this test targets:
  # size, sampled regions, and recorded mtime agree, while full content differs.
  manifest_path <- file.path(bundle, "manifest.rds")
  manifest <- readRDS(manifest_path)
  manifest$source$mtime <- .transition_source_identity(source_file)$mtime
  saveRDS(manifest, manifest_path, compress = FALSE)
  expect_true(.transition_identity_equal(
    manifest$source, .transition_source_identity(source_file)
  ))
  expect_error(
    .transition_clone_file(
      source_file, bundle, chunk_mib = 1,
      source_verification = "full"
    ),
    "full source SHA-256 differs"
  )
})

test_that("auto persists a trusted fingerprint after verifying a valid copy", {
  artifact <- tempfile(fileext = ".h5")
  manifest_path <- tempfile(fileext = ".rds")
  on.exit(unlink(c(artifact, manifest_path)), add = TRUE)
  writeBin(as.raw(seq_len(64L)), artifact)
  current <- .transition_file_fingerprint(artifact)
  manifest <- list(
    arrays_size = current$size,
    arrays_sha256 = .transition_file_sha256(artifact),
    arrays_fingerprint = list(size = current$size, mtime = current$mtime - 10)
  )
  saveRDS(manifest, manifest_path)

  expect_true(.transition_verify_arrays(
    artifact, manifest, "auto", manifest_path = manifest_path
  ))
  updated <- readRDS(manifest_path)
  expect_identical(updated$arrays_fingerprint,
                   .transition_file_fingerprint(artifact))
  expect_true(.transition_verify_arrays(
    artifact, updated, "auto",
    hash_fun = function(...) stop("unexpected second full hash")
  ))
})

test_that("sample identifiers are complete before graph-isolation auditing", {
  expect_invisible(.transition_validate_sample_ids(c("A", "B"), "slide"))
  expect_error(
    .transition_validate_sample_ids(c("A", NA_character_), "slide"),
    "missing or blank sample identifiers"
  )
  expect_error(
    .transition_validate_sample_ids(c("A", "  "), "slide"),
    "observation position\\(s\\) 2"
  )
})

test_that("index-only AnnData frames retain their declared row count", {
  source_file <- tempfile(fileext = ".h5ad")
  on.exit(unlink(source_file), add = TRUE)
  rhdf5::h5createFile(source_file)
  rhdf5::h5createGroup(source_file, "obs")
  rhdf5::h5write(c("spot1", "spot2", "spot3"), source_file, "obs/_index")

  observed <- .transition_read_frame(
    source_file, "/obs", 3L, "obs/_index", progress = FALSE
  )
  expect_s4_class(observed, "DFrame")
  expect_identical(dim(observed), c(3L, 0L))
  expect_identical(rownames(observed), c("spot1", "spot2", "spot3"))
})

test_that("AnnData dataframe _index attributes select the actual index dataset", {
  source_file <- tempfile(fileext = ".h5ad")
  on.exit(unlink(source_file), add = TRUE)
  rhdf5::h5createFile(source_file)
  rhdf5::h5createGroup(source_file, "obs")
  rhdf5::h5write(c("spot-a", "spot-b", "spot-c"), source_file,
                 "obs/barcode")
  rhdf5::h5write(c(1.5, 2.5, 3.5), source_file, "obs/score")
  file_handle <- rhdf5::H5Fopen(source_file)
  obs_handle <- rhdf5::H5Gopen(file_handle, "obs")
  rhdf5::h5writeAttribute("barcode", obs_handle, "_index")
  rhdf5::h5writeAttribute("score", obs_handle, "column-order")
  rhdf5::H5Gclose(obs_handle)
  rhdf5::H5Fclose(file_handle)

  observed <- .transition_read_frame(
    source_file, "/obs", 3L, "obs/_index", progress = FALSE
  )
  expect_identical(rownames(observed), c("spot-a", "spot-b", "spot-c"))
  expect_identical(colnames(observed), "score")
  expect_equal(as.numeric(observed$score), c(1.5, 2.5, 3.5))
  estimate <- .transition_frame_estimate(source_file, "/obs", 3L)
  expect_gt(estimate, 3 * 8)
  expect_gt(estimate, sum(64 + nchar(c("spot-a", "spot-b", "spot-c"))))
})

test_that("metadata preflight counts long index and annotation strings", {
  source_file <- tempfile(fileext = ".h5ad")
  on.exit(unlink(source_file), add = TRUE)
  rhdf5::h5createFile(source_file)
  rhdf5::h5createGroup(source_file, "obs")
  long_ids <- paste0(c("a", "b", "c"), strrep("x", 1000L))
  long_labels <- paste0(c("d", "e", "f"), strrep("y", 2000L))
  rhdf5::h5write(long_ids, source_file, "obs/_index")
  rhdf5::h5write(long_labels, source_file, "obs/label")

  estimate <- withr::with_options(
    list(CellCharterR.import.character_bytes_per_value = 4096),
    .transition_frame_estimate(source_file, "/obs", 3L)
  )
  raw_payload <- sum(nchar(long_ids, type = "bytes")) +
    sum(nchar(long_labels, type = "bytes"))
  expect_gt(estimate, raw_payload)
  expect_error(
    withr::with_options(
      list(CellCharterR.import.character_bytes_per_value = 0),
      .transition_frame_estimate(source_file, "/obs", 3L)
    ),
    "character_bytes_per_value"
  )
})

test_that("AUCell dataframe indices cannot silently reorder spots", {
  source_file <- tempfile(fileext = ".h5ad")
  on.exit(unlink(source_file), add = TRUE)
  rhdf5::h5createFile(source_file)
  rhdf5::h5createGroup(source_file, "obsm")
  rhdf5::h5createGroup(source_file, "obsm/score_aucell")
  rhdf5::h5write(c("spot2", "spot1", "spot3"), source_file,
                 "obsm/score_aucell/_index")
  rhdf5::h5write(c(0.1, 0.2, 0.3), source_file,
                 "obsm/score_aucell/program")

  info <- .transition_aucell_info(
    source_file, "obsm/score_aucell", 3L
  )
  expect_error(
    .transition_read_aucell(
      source_file, info, 3L, block_mib = 1e-6, progress = FALSE,
      expected_ids = c("spot1", "spot2", "spot3")
    ),
    "differs from `obs/_index` at position 1"
  )
})

test_that("conversion locks preserve active work and quarantine proven-dead owners", {
  output <- tempfile("conversion-lock-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)

  active <- .transition_acquire_lock(output)
  expect_error(.transition_acquire_lock(output), "Another conversion may be using")
  unlink(active, recursive = TRUE)

  dead <- file.path(output, ".conversion-lock")
  dir.create(dead)
  saveRDS(
    list(pid = 2147483000L, host = unname(Sys.info()[["nodename"]]),
         created = Sys.time() - 3600, heartbeat = Sys.time() - 3600),
    file.path(dead, "owner.rds")
  )
  expect_message(replacement <- .transition_acquire_lock(output),
                 "Quarantined stale conversion lock")
  expect_true(dir.exists(replacement))
  expect_length(Sys.glob(file.path(output, ".conversion-lock.stale-*")), 1L)
})

test_that("unsupported sparse assay and matrix obsm layouts fail preflight", {
  source_file <- tempfile(fileext = ".h5ad")
  on.exit(unlink(source_file), add = TRUE)
  rhdf5::h5createFile(source_file)
  rhdf5::h5createGroup(source_file, "X")
  rhdf5::h5createGroup(source_file, "obsm")
  rhdf5::h5write(matrix(1:4, nrow = 2L), source_file, "obsm/score_aucell")
  handle <- rhdf5::H5Fopen(source_file)
  group <- rhdf5::H5Gopen(handle, "X")
  rhdf5::h5writeAttribute("csr_matrix", group, "encoding-type")
  rhdf5::H5Gclose(group)
  rhdf5::H5Fclose(handle)

  expect_error(
    .transition_matrix_preflight(source_file, "X", "assay", "dataset"),
    "Unsupported assay layout.*csr_matrix"
  )
  expect_error(
    .transition_matrix_preflight(source_file, "obsm/score_aucell",
                                 "AUCell obsm", "dataframe"),
    "dense or sparse matrix obsm layouts"
  )
})

test_that("imported pair values are separated from graph topology semantics", {
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(x = matrix(0, nrow = 1L, ncol = 3L)),
    spatialCoords = matrix(seq_len(6), ncol = 2L)
  )
  SingleCellExperiment::colPair(spe, "graph") <- S4Vectors::SelfHits(
    from = c(1L, 2L), to = c(2L, 3L), nnode = 3L,
    value = c(0, 2), sort.by.query = FALSE
  )
  topology <- .cc_spe_graph_matrix(spe, "graph")
  expect_identical(as.numeric(topology[1, 2]), 0)
  expect_identical(as.numeric(topology[2, 3]), 1)

  SingleCellExperiment::colPair(spe, "negative") <- S4Vectors::SelfHits(
    from = 1L, to = 2L, nnode = 3L,
    value = -1, sort.by.query = FALSE
  )
  expect_error(.cc_spe_graph_matrix(spe, "negative"), "negative edge weights")
})
