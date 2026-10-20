make_aucell_spe <- function(x, delayed = FALSE, sparse = FALSE,
                            sample_id = rep("sample01", ncol(x))) {
  if (delayed) x <- DelayedArray::DelayedArray(x)
  if (sparse) x <- Matrix::Matrix(x, sparse = TRUE)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(expression = x),
    spatialCoords = cbind(seq_len(ncol(x)), 0),
    sample_id = sample_id
  )
  rownames(spe) <- rownames(x)
  colnames(spe) <- colnames(x)
  spe
}

test_that("AUCell uses the strict cutoff and normalized denominator", {
  x <- cbind(
    s1 = c(5, 4, 3, 2, 1),
    s2 = c(1, 2, 3, 4, 5)
  )
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  observed <- cellcharter_aucell(
    make_aucell_spe(x),
    list(A = c("g1", "g3"), B = "g4", C = c("g1", "g2", "g3", "g4")),
    auc_max_rank = 4, tie_method = "ordinal", threads = 2L,
    progress = FALSE
  )
  scores <- SingleCellExperiment::reducedDim(observed, "AUCell")

  # T=4 retains ranks 1:3. A's denominator is 2*4-(1+2)=5;
  # C is capped at three available ranks, with denominator 6.
  expect_equal(
    unname(scores),
    rbind(c(4 / 5, 0, 1), c(1 / 5, 2 / 3, 1 / 2)),
    tolerance = 1e-15
  )
  expect_identical(dimnames(scores), list(c("s1", "s2"), c("A", "B", "C")))
})

test_that("AUCell matches a frozen AUCell-formula tie-free oracle", {
  x <- cbind(
    s1 = 12:1,
    s2 = c(3, 10, 1, 8, 5, 12, 2, 11, 4, 9, 6, 7),
    s3 = c(8, 2, 12, 4, 10, 6, 1, 11, 5, 9, 3, 7)
  )
  rownames(x) <- paste0("g", 1:12)
  gene_sets <- list(
    small = c("g1", "g3"),
    equal = c("g2", "g4", "g6", "g8", "g10"),
    large = c(paste0("g", 1:9), "g1"),
    partial = c("g5", "absent")
  )
  expected <- rbind(
    c(8 / 9, 0.4, 1, 0.2),
    c(0, 1, 13 / 15, 0),
    c(2 / 3, 0.4, 13 / 15, 0.6)
  )
  observed <- cellcharter_aucell(
    make_aucell_spe(x), gene_sets, auc_max_rank = 5.6,
    tie_method = "ordinal", threads = 3L, progress = FALSE
  )
  expect_equal(
    unname(SingleCellExperiment::reducedDim(observed, "AUCell")),
    expected, tolerance = 1e-12
  )
})

test_that("ordinal ties follow assay feature order", {
  x <- matrix(0, nrow = 5L, ncol = 1L,
              dimnames = list(paste0("g", 1:5), "spot"))
  observed <- cellcharter_aucell(
    make_aucell_spe(x), list(first = "g1", fourth = "g4"),
    auc_max_rank = 4, tie_method = "ordinal", progress = FALSE
  )
  expect_equal(
    unname(SingleCellExperiment::reducedDim(observed, "AUCell")),
    matrix(c(1, 0), nrow = 1L), tolerance = 1e-15
  )
})

test_that("random ties are invariant to threads and chunk boundaries", {
  set.seed(9)
  x <- matrix(sample.int(4L, 1000L * 100L, replace = TRUE), 1000L, 100L,
              dimnames = list(paste0("g", 1:1000), paste0("s", 1:100)))
  gene_sets <- list(A = paste0("g", 1:100), B = paste0("g", 401:550))
  one <- cellcharter_aucell(
    make_aucell_spe(x), gene_sets, auc_max_rank = 51,
    threads = 1L, memory_mib = 65, seed = 83L, progress = FALSE
  )
  many <- cellcharter_aucell(
    make_aucell_spe(x), gene_sets, auc_max_rank = 51,
    threads = 8L, memory_mib = 128, seed = 83L, progress = FALSE
  )
  expect_identical(
    SingleCellExperiment::reducedDim(one, "AUCell"),
    SingleCellExperiment::reducedDim(many, "AUCell")
  )

  changed <- cellcharter_aucell(
    make_aucell_spe(x), gene_sets, auc_max_rank = 51,
    threads = 4L, memory_mib = 65, seed = 84L, progress = FALSE
  )
  ordinal <- cellcharter_aucell(
    make_aucell_spe(x), gene_sets, auc_max_rank = 51,
    tie_method = "ordinal", threads = 4L, memory_mib = 65,
    progress = FALSE
  )
  one_scores <- SingleCellExperiment::reducedDim(one, "AUCell")
  expect_false(identical(
    one_scores, SingleCellExperiment::reducedDim(changed, "AUCell")
  ))
  expect_false(identical(
    one_scores, SingleCellExperiment::reducedDim(ordinal, "AUCell")
  ))
  expect_true(all(one_scores >= 0 & one_scores <= 1))
})

test_that("dense, sparse, delayed, and HDF5-backed assays agree", {
  x <- matrix(seq_len(120), 20L, 6L,
              dimnames = list(paste0("g", 1:20), paste0("s", 1:6)))
  gene_sets <- list(A = paste0("g", c(2, 7, 11)), B = paste0("g", c(5, 15)))
  score <- function(spe) {
    spe <- cellcharter_aucell(
      spe, gene_sets, auc_max_rank = 8, tie_method = "ordinal",
      threads = 2L, progress = FALSE
    )
    SingleCellExperiment::reducedDim(spe, "AUCell")
  }
  expected <- score(make_aucell_spe(x))
  expect_identical(score(make_aucell_spe(x, delayed = TRUE)), expected)
  expect_identical(score(make_aucell_spe(x, sparse = TRUE)), expected)

  path <- tempfile(fileext = ".h5")
  h5 <- HDF5Array::writeHDF5Array(x, filepath = path, name = "expression",
                                  with.dimnames = TRUE)
  expect_identical(score(make_aucell_spe(h5)), expected)
})

test_that("HDF5 scoring uses multiple bounded chunks", {
  x <- matrix(seq_len(1000L * 100L), 1000L, 100L,
              dimnames = list(paste0("g", 1:1000), paste0("s", 1:100)))
  gene_sets <- list(A = paste0("g", 1:100), B = paste0("g", 501:650))
  chunk_width <- CellCharterR:::.cc_aucell_chunk_width(1000L, 2L, 1L, 65)
  expect_lt(chunk_width, ncol(x))
  dense <- cellcharter_aucell(
    make_aucell_spe(x), gene_sets, auc_max_rank = 51,
    tie_method = "ordinal", memory_mib = 65, threads = 1L,
    progress = FALSE
  )
  path <- tempfile(fileext = ".h5")
  h5 <- HDF5Array::writeHDF5Array(x, filepath = path, name = "expression",
                                  with.dimnames = TRUE)
  backed <- cellcharter_aucell(
    make_aucell_spe(h5), gene_sets, auc_max_rank = 51,
    tie_method = "ordinal", memory_mib = 65, threads = 1L,
    progress = FALSE
  )
  expect_identical(
    SingleCellExperiment::reducedDim(backed, "AUCell"),
    SingleCellExperiment::reducedDim(dense, "AUCell")
  )
})

test_that("AUCell attaches only the completed reduced dimension", {
  x <- matrix(seq_len(50), 10L, 5L,
              dimnames = list(paste0("g", 1:10), paste0("s", 1:5)))
  spe <- make_aucell_spe(x)
  S4Vectors::metadata(spe)$sentinel <- "unchanged"
  first <- cellcharter_aucell(
    spe, list(A = c("g1", "g2")), auc_max_rank = 5,
    tie_method = "ordinal", progress = FALSE
  )
  second <- cellcharter_aucell(
    first, list(B = c("g3", "g4")), auc_max_rank = 5,
    tie_method = "ordinal", progress = FALSE
  )
  expect_identical(S4Vectors::metadata(second), S4Vectors::metadata(spe))
  expect_identical(colnames(SingleCellExperiment::reducedDim(second, "AUCell")), "B")
  expect_false("AUCell" %in% SingleCellExperiment::reducedDimNames(spe))
})

test_that("late chunk failure leaves the input unchanged and creates no files", {
  x <- matrix(seq_len(1000L * 100L), 1000L, 100L,
              dimnames = list(paste0("g", 1:1000), paste0("s", 1:100)))
  x[1L, 100L] <- Inf
  spe <- make_aucell_spe(x)
  sentinel <- matrix(seq_len(200), 100L, 2L,
                     dimnames = list(colnames(spe), c("old1", "old2")))
  SingleCellExperiment::reducedDim(spe, "AUCell") <- sentinel
  before <- spe
  work <- tempfile("aucell-no-files-")
  dir.create(work)
  old <- setwd(work)
  on.exit(setwd(old), add = TRUE)
  expect_error(
    cellcharter_aucell(
      spe, list(A = paste0("g", 1:100)), auc_max_rank = 51,
      memory_mib = 65, threads = 1L, progress = FALSE
    ),
    "non-finite"
  )
  expect_identical(spe, before)
  expect_identical(list.files(work, all.files = TRUE, no.. = TRUE), character())
})

test_that("scores preserve repeated barcodes and feed downstream APIs", {
  x <- matrix(
    c(9, 8, 1, 0, 0, 1, 8, 9,
      8, 9, 0, 1, 1, 0, 9, 8),
    nrow = 4L,
    dimnames = list(paste0("g", 1:4), NULL)
  )
  spe <- make_aucell_spe(
    x, sample_id = rep(c("patient1", "patient2"), each = 2L)
  )
  colnames(spe) <- rep(c("barcode1", "barcode2"), 2L)
  rownames(SpatialExperiment::spatialCoords(spe)) <- colnames(spe)
  scored <- cellcharter_aucell(
    spe, list(A = c("g1", "g2"), B = c("g3", "g4")),
    auc_max_rank = 3, tie_method = "ordinal", progress = FALSE
  )
  expect_identical(
    rownames(SingleCellExperiment::reducedDim(scored, "AUCell")),
    colnames(spe)
  )
  assigned <- cellcharter_assign_dominant_program(
    scored, scores = "AUCell", block = NULL, programs = c("A", "B"),
    output = "dominant", progress = FALSE
  )
  expect_identical(length(assigned$dominant), ncol(spe))
})

test_that("AUCell rejects unsafe inputs before scoring", {
  x <- matrix(seq_len(50), 10L, 5L,
              dimnames = list(paste0("g", 1:10), paste0("s", 1:5)))
  spe <- make_aucell_spe(x)
  expect_error(cellcharter_aucell(spe, list(A = "missing"), progress = FALSE),
               "no features")
  expect_error(cellcharter_aucell(
    spe, list(A = c("g1", paste0("missing", 1:5))),
    auc_max_rank = 2, progress = FALSE
  ), "80%")
  expect_error(cellcharter_aucell(
    spe, list(A = ""), progress = FALSE
  ), "character vector")
  expect_error(cellcharter_aucell(
    spe, list(A = "g1"), auc_max_rank = 1, progress = FALSE
  ), "between 2")
  expect_error(cellcharter_aucell(
    spe, list(A = "g1"), norm_auc = FALSE, progress = FALSE
  ), "normalized")
  expect_error(cellcharter_aucell(
    spe, list(A = "g1"), threads = 33, progress = FALSE
  ), "cannot exceed 32")

  negative <- x
  negative[] <- -negative
  expect_no_error(cellcharter_aucell(
    make_aucell_spe(negative), list(A = "g1"), auc_max_rank = 2,
    tie_method = "ordinal", progress = FALSE
  ))
})
