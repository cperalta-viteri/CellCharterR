.make_state_import_directory <- function() {
  source <- tempfile(fileext = ".h5ad")
  path <- tempfile("cellcharterr-state-")

  rhdf5::h5createFile(source)
  for (group in c("obs", "var", "obsm")) {
    rhdf5::h5createGroup(source, group)
  }
  rhdf5::h5write(matrix(seq_len(12L), nrow = 3L), source, "X")
  rhdf5::h5write(paste0("spot-", seq_len(4L)), source, "obs/_index")
  rhdf5::h5write(rep(c("A", "B"), each = 2L), source,
                 "obs/slide_id_concat")
  rhdf5::h5write(paste0("gene-", seq_len(3L)), source, "var/_index")
  rhdf5::h5write(rbind(seq_len(4L), seq_len(4L) + 10L), source,
                 "obsm/spatial")

  spe <- cellcharter_import_h5ad(
    source,
    path,
    obs_columns = "slide_id_concat",
    var_columns = character(),
    reduced_dims = character(),
    aucell_path = NULL,
    sample_column = "slide_id_concat",
    materialize_block_mib = 1e-6,
    progress = FALSE
  )
  unlink(source)
  list(path = path, spe = spe)
}

.state_file_identity <- function(path) {
  info <- file.info(path)
  inode <- if (nzchar(Sys.which("stat"))) {
    output <- system2("stat", c("-c", "%i", path), stdout = TRUE)
    if (identical(attr(output, "status"), NULL)) output[[1L]] else NA_character_
  } else {
    NA_character_
  }
  list(
    inode = inode,
    size = unname(info$size),
    mtime = unclass(info$mtime),
    sha256 = cellcharter_file_sha256(path)
  )
}

test_that("in-place SPE saves preserve the immutable assay and restore analysis state", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  arrays <- file.path(path, "arrays.h5")
  before <- .state_file_identity(arrays)
  pristine <- cellcharter_load_spe(
    path, verify = "none", restore_state = FALSE,
    materialize_block_mib = 1e-6, progress = FALSE
  )

  # Exercise the mutable compartments used by downstream CellCharterR work.
  SummarizedExperiment::colData(spe)$domain <- factor(c("1", "1", "2", "2"))
  SummarizedExperiment::rowData(spe)$selected <- c(TRUE, FALSE, TRUE)
  SingleCellExperiment::reducedDim(spe, "X_cellcharter_R") <- matrix(
    seq_len(8L), nrow = 4L,
    dimnames = list(NULL, c("feature-1", "feature-2"))
  )
  SingleCellExperiment::colPair(spe, "cellcharter_graph") <- S4Vectors::SelfHits(
    from = c(1L, 2L, 3L), to = c(2L, 1L, 4L), nnode = 4L,
    value = c(1, 1, 0.5), sort.by.query = FALSE
  )
  S4Vectors::metadata(spe)$CellCharterR <- list(
    enrichment = list(dominant_program = list(observed = matrix(1:4, 2L)))
  )

  # Repeated display labels are legitimate for large imported cohorts. The
  # coordinate matrix may retain a separate, unique AnnData observation index.
  colnames(spe) <- rep(c("shared-spot", "other-spot"), 2L)
  coordinates <- SpatialExperiment::spatialCoords(spe)
  rownames(coordinates) <- paste0("coordinate-id-", seq_len(ncol(spe)))
  coordinates[, 1L] <- coordinates[, 1L] + 100
  SpatialExperiment::spatialCoords(spe) <- coordinates

  cellcharter_save_spe(
    spe, path = path, retain_generations = 1L, progress = FALSE
  )

  after <- .state_file_identity(arrays)
  expect_identical(after$size, before$size)
  expect_identical(after$mtime, before$mtime)
  expect_identical(after$sha256, before$sha256)
  if (!is.na(before$inode) && !is.na(after$inode)) {
    expect_identical(after$inode, before$inode)
  }

  # The state overlay must not create a second HDF5 assay container. Comparing
  # file sizes is misleading for tiny fixtures because an RDS header can be
  # larger than this toy HDF5 file.
  h5_files <- list.files(
    path, pattern = "[.]h5$", recursive = TRUE, full.names = TRUE
  )
  expect_identical(normalizePath(h5_files), normalizePath(arrays))

  restored <- cellcharter_load_spe(
    path, verify = "none", restore_state = TRUE,
    materialize_block_mib = 1e-6, progress = FALSE
  )
  expect_identical(as.character(SummarizedExperiment::colData(restored)$domain),
                   c("1", "1", "2", "2"))
  expect_identical(SummarizedExperiment::rowData(restored)$selected,
                   c(TRUE, FALSE, TRUE))
  expect_equal(
    SingleCellExperiment::reducedDim(restored, "X_cellcharter_R"),
    SingleCellExperiment::reducedDim(spe, "X_cellcharter_R")
  )
  expect_equal(
    as.matrix(SingleCellExperiment::colPair(
      restored, "cellcharter_graph", asSparse = TRUE
    )),
    as.matrix(SingleCellExperiment::colPair(
      spe, "cellcharter_graph", asSparse = TRUE
    ))
  )
  expect_equal(
    S4Vectors::metadata(restored)$CellCharterR$enrichment$dominant_program,
    S4Vectors::metadata(spe)$CellCharterR$enrichment$dominant_program
  )
  expect_identical(colnames(restored), colnames(spe))
  expect_identical(
    rownames(SpatialExperiment::spatialCoords(restored)),
    rownames(SpatialExperiment::spatialCoords(spe))
  )
  expect_identical(
    SpatialExperiment::spatialCoords(restored),
    SpatialExperiment::spatialCoords(spe)
  )
  expect_true(methods::validObject(restored))

  # Opting out of restoration must reconstruct the immutable import state.
  base <- cellcharter_load_spe(
    path, verify = "none", restore_state = FALSE,
    materialize_block_mib = 1e-6, progress = FALSE
  )
  expect_equal(
    as.matrix(SummarizedExperiment::assay(base, "expression")),
    as.matrix(SummarizedExperiment::assay(pristine, "expression"))
  )
  expect_identical(colnames(base), colnames(pristine))
  expect_identical(
    SpatialExperiment::spatialCoords(base),
    SpatialExperiment::spatialCoords(pristine)
  )
  expect_false("domain" %in% colnames(SummarizedExperiment::colData(base)))
  expect_false("X_cellcharter_R" %in%
                 SingleCellExperiment::reducedDimNames(base))
  expect_false("cellcharter_graph" %in%
                 SingleCellExperiment::colPairNames(base))
})

test_that("successive SPE saves retain one valid generation by default", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  SummarizedExperiment::colData(spe)$iteration <- 1L
  cellcharter_save_spe(spe, path, retain_generations = 1L, progress = FALSE)
  orphan <- file.path(path, "state", "generations", "orphan-generation")
  dir.create(orphan)
  saveRDS(list(orphan = TRUE), file.path(orphan, "state.rds"))
  SummarizedExperiment::colData(spe)$iteration <- 2L
  cellcharter_save_spe(spe, path, retain_generations = 1L, progress = FALSE)

  generations <- list.dirs(
    file.path(path, "state", "generations"),
    recursive = FALSE, full.names = TRUE
  )
  expect_length(generations, 1L)
  expect_false(dir.exists(orphan))
  restored <- cellcharter_load_spe(
    path, verify = "none", restore_state = TRUE,
    materialize_block_mib = 1e-6, progress = FALSE
  )
  expect_identical(SummarizedExperiment::colData(restored)$iteration,
                   rep(2L, ncol(restored)))
})

test_that("in-place SPE saving rejects assay changes that it cannot persist", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  reordered <- spe[c(3L, 1L, 2L), ]
  expect_error(
    cellcharter_save_spe(reordered, path, progress = FALSE),
    "subset, reordered, or transformed"
  )

  transformed <- spe
  SummarizedExperiment::assay(transformed, "expression") <-
    SummarizedExperiment::assay(transformed, "expression") + 1
  expect_error(
    cellcharter_save_spe(transformed, path, progress = FALSE),
    "subset, reordered, or transformed"
  )

  additional <- spe
  SummarizedExperiment::assay(additional, "extra") <-
    SummarizedExperiment::assay(additional, "expression")
  expect_error(
    cellcharter_save_spe(additional, path, progress = FALSE),
    "exactly the imported assay"
  )
})

test_that("AutoK checkpoints remain separate from in-place SPE state", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  moved <- paste0(path, "-moved")
  spe <- fixture$spe
  checkpoint_dir <- tempfile("separate-autok-")
  dir.create(checkpoint_dir)
  saveRDS(list(model = "kept-separate"), file.path(checkpoint_dir, "fit.rds"))
  on.exit(unlink(c(path, moved, checkpoint_dir), recursive = TRUE, force = TRUE),
          add = TRUE)

  checkpoint_path <- normalizePath(checkpoint_dir)
  SummarizedExperiment::colData(spe)$autok_domain <- factor(c(1L, 1L, 2L, 2L))
  S4Vectors::metadata(spe)$CellCharterR <- list(
    autok = list(autok_domain = structure(list(
      checkpoint_dir = checkpoint_path,
      checkpoint_mode = "none",
      embedded_best_artifacts = list(model = "self-contained-but-omitted")
    ), class = "cellcharter_autok")),
    fixed_k = list(autok_domain = list(
      checkpoint = file.path(checkpoint_path, "fixed_k2.rds"),
      model = "also-omitted"
    )),
    stages = list(`clustering:autok_domain` = list(
      parameters = list(checkpoint_dir = checkpoint_path)
    ), `clustering:orphan_legacy_domain` = list(
      parameters = list(checkpoint_dir = "/deleted/legacy/checkpoint")
    ), enrichment = list(result = "portable")),
    labels_are_in_colData = TRUE
  )
  cellcharter_save_spe(spe, path, progress = FALSE)

  generation <- readRDS(file.path(path, "state", "current.rds"))$generation
  payloads <- list.files(
    file.path(path, "state", "generations", generation),
    recursive = TRUE, all.files = TRUE, no.. = TRUE
  )
  expect_identical(payloads, "state.rds")
  expect_true(file.rename(path, moved))
  unlink(checkpoint_dir, recursive = TRUE, force = TRUE)
  restored <- cellcharter_load_spe(
    moved, verify = "none", progress = FALSE
  )
  expect_identical(as.character(restored$autok_domain), c("1", "1", "2", "2"))
  expect_null(S4Vectors::metadata(restored)$CellCharterR$autok)
  expect_null(S4Vectors::metadata(restored)$CellCharterR$fixed_k)
  expect_null(
    S4Vectors::metadata(restored)$CellCharterR$stages$`clustering:autok_domain`
  )
  expect_null(
    S4Vectors::metadata(restored)$CellCharterR$stages$`clustering:orphan_legacy_domain`
  )
  expect_identical(
    S4Vectors::metadata(restored)$CellCharterR$stages$enrichment$result,
    "portable"
  )
  expect_false(file.exists(checkpoint_path))
})

test_that("in-place SPE state rejects delayed reduced dimensions", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  SingleCellExperiment::reducedDim(spe, "delayed") <-
    DelayedArray::DelayedArray(matrix(seq_len(8L), nrow = 4L))
  expect_error(
    cellcharter_save_spe(spe, path, progress = FALSE),
    "native in-memory reducedDims"
  )
})

test_that("in-place SPE state rejects nested external HDF5 dependencies", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  external <- tempfile("external-state-value-", fileext = ".h5")
  on.exit(unlink(c(path, external), recursive = TRUE, force = TRUE), add = TRUE)

  delayed <- HDF5Array::writeHDF5Array(
    matrix(seq_len(8L), nrow = 4L), filepath = external, name = "external"
  )
  S4Vectors::metadata(spe)$nested <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(external = delayed)
  )
  expect_error(
    cellcharter_save_spe(spe, path, progress = FALSE),
    "delayed|HDF5-backed|external"
  )

  S4Vectors::metadata(spe)$nested <- NULL
  attributed <- 1
  attr(attributed, "external_hdf5") <- delayed
  S4Vectors::metadata(spe)$attributed <- attributed
  expect_error(
    cellcharter_save_spe(spe, path, progress = FALSE),
    "delayed|HDF5-backed|external"
  )

  S4Vectors::metadata(spe)$attributed <- NULL
  SummarizedExperiment::colData(spe)$external <- I(rep(list(delayed), ncol(spe)))
  expect_error(
    cellcharter_save_spe(spe, path, progress = FALSE),
    "delayed|HDF5-backed|external"
  )

  SummarizedExperiment::colData(spe)$external <- NULL
  pair <- S4Vectors::SelfHits(1L, 2L, nnode = ncol(spe), value = 1)
  S4Vectors::metadata(pair)$nested <-
    SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(external = delayed)
    )
  SingleCellExperiment::colPair(spe, "external") <- pair
  expect_error(
    cellcharter_save_spe(spe, path, progress = FALSE),
    "delayed|HDF5-backed|external"
  )
})

test_that("managed SPE state paths cannot escape through symbolic links", {
  skip_on_os("windows")
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  external_state <- tempfile("external-spe-state-")
  on.exit(unlink(c(path, external_state), recursive = TRUE, force = TRUE),
          add = TRUE)

  cellcharter_save_spe(spe, path, progress = FALSE)
  expect_true(file.rename(file.path(path, "state"), external_state))
  expect_true(file.symlink(external_state, file.path(path, "state")))
  expect_error(
    cellcharter_load_spe(path, verify = "none", progress = FALSE),
    "symbolic link"
  )
})

test_that("corrupt SPE state pointers fail closed", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  SummarizedExperiment::colData(spe)$saved <- TRUE
  cellcharter_save_spe(spe, path, retain_generations = 1L, progress = FALSE)
  pointer <- file.path(path, "state", "current.rds")
  expect_true(file.exists(pointer))
  writeBin(charToRaw("not an RDS"), pointer)

  expect_error(
    cellcharter_load_spe(
      path, verify = "none", restore_state = TRUE,
      materialize_block_mib = 1e-6, progress = FALSE
    ),
    "state|current|RDS|read|corrupt|invalid|input format",
    ignore.case = TRUE
  )
  # Explicitly requesting the pristine import remains a recovery path.
  expect_s4_class(
    cellcharter_load_spe(
      path, verify = "none", restore_state = FALSE,
      materialize_block_mib = 1e-6, progress = FALSE
    ),
    "SpatialExperiment"
  )
})

test_that("legacy import manifests upgrade in place without rewriting arrays", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  moved <- paste0(path, "-moved")
  on.exit(unlink(c(path, moved), recursive = TRUE, force = TRUE), add = TRUE)

  arrays <- file.path(path, "arrays.h5")
  before <- .state_file_identity(arrays)
  manifest_path <- file.path(path, "manifest.rds")
  manifest <- readRDS(manifest_path)
  legacy_spec <- manifest$import_spec
  manifest$format_version <- 1L
  manifest$import_spec <- NULL
  manifest$arrays_sha256 <- NULL
  manifest$arrays_fingerprint <- NULL
  saveRDS(manifest, manifest_path)

  SummarizedExperiment::colData(spe)$legacy_domain <- factor(c(1, 1, 2, 2))
  SingleCellExperiment::colPair(spe, "legacy_graph") <- S4Vectors::SelfHits(
    c(1L, 2L), c(2L, 1L), nnode = ncol(spe), value = c(1, 1)
  )
  cellcharter_save_spe(
    spe, path, legacy_import_spec = legacy_spec, progress = FALSE
  )

  after <- .state_file_identity(arrays)
  expect_identical(after, before)
  upgraded <- readRDS(manifest_path)
  expect_identical(upgraded$format_version, 3L)
  expect_identical(upgraded$import_spec, legacy_spec)
  expect_match(upgraded$arrays_sha256, "^[0-9a-f]{64}$")
  expect_true(file.rename(path, moved))
  restored <- cellcharter_load_spe(moved, verify = "full", progress = FALSE)
  expect_identical(as.character(restored$legacy_domain), c("1", "1", "2", "2"))
  expect_equal(
    as.matrix(SingleCellExperiment::colPair(restored, "legacy_graph", asSparse = TRUE)),
    as.matrix(SingleCellExperiment::colPair(spe, "legacy_graph", asSparse = TRUE))
  )
})

test_that("all loads honor an active writer lock", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  lock <- CellCharterR:::.transition_acquire_lock(path)
  expect_error(
    cellcharter_load_spe(
      path, verify = "none", restore_state = FALSE, progress = FALSE
    ),
    "Another conversion|using this bundle|lock"
  )
  CellCharterR:::.transition_release_lock(lock)
})

test_that("retained state generations support explicit verified rollback", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  spe$iteration <- 1L
  cellcharter_save_spe(spe, path, retain_generations = 2L, progress = FALSE)
  spe$iteration <- 2L
  cellcharter_save_spe(spe, path, retain_generations = 2L, progress = FALSE)
  pointer <- readRDS(file.path(path, "state", "current.rds"))
  expect_length(pointer$history, 2L)
  previous <- pointer$history[[2L]]
  writeBin(
    charToRaw("damaged"),
    file.path(path, "state", "generations", pointer$generation, "state.rds")
  )
  expect_error(
    cellcharter_load_spe(path, verify = "none", progress = FALSE),
    "checksum|size|invalid|corrupt",
    ignore.case = TRUE
  )
  restored <- cellcharter_load_spe(
    path, verify = "none", state_generation = previous, progress = FALSE
  )
  expect_identical(restored$iteration, rep(1L, ncol(restored)))
})

test_that("retained generation records cannot escape their managed directory", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  cellcharter_save_spe(spe, path, retain_generations = 2L, progress = FALSE)
  pointer_path <- file.path(path, "state", "current.rds")
  pointer <- readRDS(pointer_path)
  pointer$generations[[pointer$generation]]$generation <- "../../../outside-state"
  saveRDS(pointer, pointer_path)
  expect_error(
    cellcharter_load_spe(path, verify = "none", progress = FALSE),
    "invalid|does not match|managed child",
    ignore.case = TRUE
  )
})

test_that("internal reduced-dimension row names round-trip exactly", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)

  internal <- matrix(
    seq_len(8L), nrow = ncol(spe),
    dimnames = list(paste0("internal-", seq_len(ncol(spe))), c("a", "b"))
  )
  SingleCellExperiment::reducedDim(spe, "internal", withDimnames = FALSE) <- internal
  cellcharter_save_spe(spe, path, progress = FALSE)
  restored <- cellcharter_load_spe(path, verify = "none", progress = FALSE)
  expect_identical(
    SingleCellExperiment::reducedDim(restored, "internal", withDimnames = FALSE),
    internal
  )
})

test_that("managed manifest and dangling pointer symlinks fail closed", {
  skip_on_os("windows")
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  fixture <- .make_state_import_directory()
  path <- fixture$path
  spe <- fixture$spe
  outside <- tempfile("outside-manifest-", fileext = ".rds")
  on.exit(unlink(c(path, outside), recursive = TRUE, force = TRUE), add = TRUE)

  manifest_path <- file.path(path, "manifest.rds")
  expect_true(file.rename(manifest_path, outside))
  expect_true(file.symlink(outside, manifest_path))
  expect_error(
    cellcharter_load_spe(path, verify = "none", progress = FALSE),
    "symbolic link"
  )
  unlink(manifest_path)
  expect_true(file.rename(outside, manifest_path))
  cellcharter_save_spe(spe, path, progress = FALSE)
  pointer <- file.path(path, "state", "current.rds")
  unlink(pointer)
  expect_true(file.symlink(file.path(path, "missing-pointer.rds"), pointer))
  expect_error(
    cellcharter_load_spe(path, verify = "none", progress = FALSE),
    "symbolic link"
  )
})

test_that("real mouse-spleen state remains movable and cleans temporary files", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("HDF5Array")
  source <- system.file(
    "extdata", "codex_mouse_spleen_tiny_raw.h5ad",
    package = "CellCharterR", mustWork = TRUE
  )
  path <- tempfile("real-mouse-spleen-state-")
  moved <- paste0(path, "-moved")
  on.exit(unlink(c(path, moved), recursive = TRUE, force = TRUE), add = TRUE)

  spe <- cellcharter_import_h5ad(
    source, path, reduced_dims = character(), aucell_path = NULL,
    sample_column = "sample", progress = FALSE
  )
  reference <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
  expect_identical(dim(spe), c(29L, 1000L))
  for (name in SingleCellExperiment::reducedDimNames(reference)) {
    SingleCellExperiment::reducedDim(spe, name, withDimnames = FALSE) <-
      SingleCellExperiment::reducedDim(reference, name, withDimnames = FALSE)
  }
  SingleCellExperiment::colPair(spe, "cellcharter_graph") <-
    SingleCellExperiment::colPair(reference, "cellcharter_graph")
  spe$r_domain_k4 <- reference$r_domain_k4
  S4Vectors::metadata(spe)$real_mouse_spleen_acceptance <- list(spots = ncol(spe))

  arrays <- file.path(path, "arrays.h5")
  before <- cellcharter_file_sha256(arrays)
  generations <- file.path(path, "state", "generations")
  stale <- file.path(generations, ".staging-interrupted")
  dir.create(stale, recursive = TRUE)
  writeLines("incomplete", file.path(stale, "state.tmp"))
  cellcharter_save_spe(spe, path, progress = FALSE)
  expect_false(dir.exists(stale))
  artifacts <- list.files(
    path, recursive = TRUE, all.files = TRUE, include.dirs = TRUE,
    no.. = TRUE
  )
  expect_false(any(grepl("(^|/)(\\.staging-|\\.tmp-)", artifacts)))
  expect_false(file.exists(file.path(path, ".conversion-lock")))
  expect_true(file.rename(path, moved))
  restored <- cellcharter_load_spe(moved, verify = "full", progress = FALSE)
  expect_identical(cellcharter_file_sha256(file.path(moved, "arrays.h5")), before)
  expect_identical(as.character(restored$r_domain_k4), as.character(reference$r_domain_k4))
  expect_identical(
    SingleCellExperiment::reducedDim(restored, "X_cellcharter", withDimnames = FALSE),
    SingleCellExperiment::reducedDim(reference, "X_cellcharter", withDimnames = FALSE)
  )
  expect_equal(
    as.matrix(SingleCellExperiment::colPair(restored, "cellcharter_graph", asSparse = TRUE)),
    as.matrix(SingleCellExperiment::colPair(reference, "cellcharter_graph", asSparse = TRUE))
  )
  expect_true(methods::validObject(restored))
  moved_artifacts <- list.files(
    moved, recursive = TRUE, all.files = TRUE, include.dirs = TRUE,
    no.. = TRUE
  )
  expect_false(any(grepl("(^|/)(\\.staging-|\\.tmp-)", moved_artifacts)))
  expect_false(file.exists(file.path(moved, ".conversion-lock")))
})
