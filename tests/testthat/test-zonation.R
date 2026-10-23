make_zonation_spe <- function(use_alt_exp = FALSE) {
  blocks <- rep(paste0("p", 1:8), each = 9L)
  zones <- rep(rep(c("Pre", "Early", "Late"), each = 3L), 8L)
  base <- matrix(-1, nrow = length(blocks), ncol = 5L)
  colnames(base) <- c("Pre-EMT", "Early-EMT", "Late-EMT", "A", "B")
  base[cbind(seq_len(nrow(base)), match(zones, c("Pre", "Early", "Late")))] <- 3
  patient <- match(blocks, unique(blocks))
  base[, "A"] <- ifelse(zones == "Early", 3, -0.5) + patient / 20
  base[, "B"] <- ifelse(zones == "Late", 3, -0.5) - patient / 30
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(0, 1L, length(blocks))),
    colData = S4Vectors::DataFrame(
      patient = blocks,
      risk = rep(rep(c("High", "Low"), each = 4L), each = 9L),
      site = rep(rep(c("AR", "BM"), 4L), each = 9L)
    ),
    spatialCoords = cbind(seq_along(blocks), 0),
    sample_id = blocks
  )
  colnames(spe) <- paste0("spot", seq_along(blocks))
  rownames(SpatialExperiment::spatialCoords(spe)) <- colnames(spe)
  rownames(base) <- colnames(spe)
  if (use_alt_exp) {
    alternative <- SingleCellExperiment::SingleCellExperiment(
      assays = list(score = t(base))
    )
    SingleCellExperiment::altExp(spe, "AUCell") <- alternative
  } else {
    SingleCellExperiment::reducedDim(spe, "AUCell") <- DelayedArray::DelayedArray(base)
  }
  spe
}

test_that("zonation public APIs name their experiment argument spe", {
  expect_identical(names(formals(cellcharter_assign_dominant_program))[[1L]], "spe")
  expect_identical(names(formals(cellcharter_zonation))[[1L]], "spe")
})

test_that("dominant programs are assigned block-wise without running zonation", {
  spe <- make_zonation_spe()
  observed <- cellcharter_assign_dominant_program(
    spe = spe, block = "patient",
    programs = c("Pre-EMT", "Early-EMT", "Late-EMT"),
    output = "dominant_emt", progress = FALSE
  )
  expect_identical(
    as.character(observed$dominant_emt[1:9]),
    rep(c("Pre-EMT", "Early-EMT", "Late-EMT"), each = 3L)
  )
  record <- S4Vectors::metadata(observed)$CellCharterR$stages$
    `dominant_program:dominant_emt`
  expect_identical(record$parameters$block, "patient")
  expect_identical(record$parameters$score_source, "reducedDim:AUCell")
  expect_error(
    cellcharter_assign_dominant_program(
      observed, block = "patient", output = "dominant_emt", progress = FALSE
    ),
    "already contains"
  )

  alternative <- cellcharter_assign_dominant_program(
    make_zonation_spe(TRUE), block = "patient",
    programs = c("Pre-EMT", "Early-EMT", "Late-EMT"),
    output = "dominant_emt", progress = FALSE
  )
  expect_identical(observed$dominant_emt, alternative$dominant_emt)
})

test_that("zonation is deterministic, block-aware, and stores audit data", {
  spe <- make_zonation_spe()
  set.seed(812L)
  before <- .Random.seed
  observed <- cellcharter_zonation(
    spe = spe, block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 4L, seed = 19L
  )
  expect_identical(.Random.seed, before)
  record <- S4Vectors::metadata(observed)$CellCharterR$results$zonation
  expect_identical(record$value$schema, "cellcharter_zonation_v2")
  expect_identical(record$value$block_mode, "blocked")
  expect_identical(record$value$response_programs, c("A", "B"))
  expect_equal(nrow(record$value$block_results), 8L * 5L * 3L)
  expect_true(all(record$value$zone_audit$retained))
  expect_identical(record$parameters$min_spots_per_zone, 2L)
  expect_true("cellcharter_zone" %in% colnames(SummarizedExperiment::colData(observed)))

  repeated <- cellcharter_zonation(
    spe, block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 4L, seed = 19L
  )
  expect_identical(
    S4Vectors::metadata(observed)$CellCharterR$results$zonation$value$block_results,
    S4Vectors::metadata(repeated)$CellCharterR$results$zonation$value$block_results
  )
})

test_that("zone definitions and measured response programs are independent", {
  observed <- cellcharter_zonation(
    make_zonation_spe(), block = "patient",
    zone_programs = c("Pre-EMT", "Early-EMT", "Late-EMT"),
    composition_programs = c("A", "B"),
    min_spots_per_zone = 2L, sample_size = 2L, n_iterations = 3L
  )
  result <- S4Vectors::metadata(observed)$CellCharterR$results$zonation$value
  expect_identical(result$composition_programs, c("A", "B"))
  expect_identical(result$response_programs, c("A", "B"))
  expect_identical(unique(result$block_results$program), c("A", "B"))
  expect_equal(nrow(result$block_results), 8L * 2L * 3L)
})

test_that("block is optional for pooled descriptive zonation", {
  observed <- cellcharter_zonation(
    make_zonation_spe(), block = NULL,
    min_spots_per_zone = 2L, sample_size = 2L, n_iterations = 3L
  )
  record <- S4Vectors::metadata(observed)$CellCharterR$results$zonation
  expect_null(record$parameters$block)
  expect_identical(record$value$block_mode, "pooled")
  expect_identical(unique(record$value$block_results$block), ".all")
  expect_false("block" %in% names(record$input_fingerprints))

  dominant <- cellcharter_assign_dominant_program(
    make_zonation_spe(), block = NULL, programs = c("A", "B"),
    progress = FALSE
  )
  expect_true("dominant_program" %in% colnames(SummarizedExperiment::colData(dominant)))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(
      cellcharter_plot_zonation(observed, compare = "none"), "ggplot"
    )
    expect_error(
      cellcharter_plot_zonation(observed, compare = "groups", group = "risk"),
      "require zonation computed with a `block`"
    )
  }
})

test_that("zonation accepts sample-local Visium barcodes", {
  spe <- make_zonation_spe()
  colnames(spe) <- rep(paste0("barcode", 1:9), 8L)
  rownames(SpatialExperiment::spatialCoords(spe)) <- colnames(spe)
  expect_no_error(
    observed <- cellcharter_zonation(
      spe, block = "sample_id", composition_programs = c("A", "B"),
      min_spots_per_zone = 2L, sample_size = 2L, n_iterations = 2L
    )
  )
  expect_identical(
    S4Vectors::metadata(observed)$CellCharterR$results$zonation$value$block_mode,
    "blocked"
  )
})

test_that("altExp and delayed reducedDim score sources agree", {
  reduced <- cellcharter_zonation(
    make_zonation_spe(FALSE), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 3L
  )
  alternative <- cellcharter_zonation(
    make_zonation_spe(TRUE), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 3L
  )
  expect_equal(
    S4Vectors::metadata(reduced)$CellCharterR$results$zonation$value$block_results,
    S4Vectors::metadata(alternative)$CellCharterR$results$zonation$value$block_results
  )
})

test_that("authors-style numerical mechanics match a hand-derived fixture", {
  composition <- CellCharterR:::.cc_zonation_positive_composition(
    rbind(c(1, 1, -1), c(2, 0, 2), c(0, 0, 0)),
    threshold = 0.1
  )
  expect_equal(
    composition,
    rbind(c(0.5, 0.5, 0), c(0.5, 0, 0.5), c(0, 0, 0))
  )
  rarefied <- CellCharterR:::.cc_zonation_rarefy(
    list(Pre = c(A = 2, B = 0), Early = c(A = 0, B = 2)),
    programs = c("A", "B"), sample_size = 2L,
    n_iterations = 3L, seed = 7L
  )
  expect_identical(rarefied$depth, 2L)
  expect_equal(
    rarefied$value,
    matrix(c(2, 0, 0, 2), nrow = 2L,
           dimnames = list(c("A", "B"), c("Pre", "Early")))
  )

  observed <- cellcharter_zonation(
    make_zonation_spe(), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 4L, seed = 19L
  )
  value <- S4Vectors::metadata(observed)$CellCharterR$results$zonation$value
  one <- value$block_results[value$block_results$block == "p1", ]
  expect_equal(
    one$enrichment_z[one$program == "A"],
    c(-1 / sqrt(3), 2 / sqrt(3), -1 / sqrt(3)),
    tolerance = 1e-10
  )
  expect_equal(
    one$enrichment_z[one$program == "B"],
    c(-1 / sqrt(3), -1 / sqrt(3), 2 / sqrt(3)),
    tolerance = 1e-10
  )
  expect_identical(
    as.character(SummarizedExperiment::colData(observed)$cellcharter_zone[1:9]),
    rep(c("Pre-EMT", "Early-EMT", "Late-EMT"), each = 3L)
  )
})

test_that("zone cutoff is strict and invalid inputs fail early", {
  spe <- make_zonation_spe()
  observed <- cellcharter_zonation(
    spe, block = "patient", min_spots_per_zone = 3L,
    sample_size = 2L, n_iterations = 2L
  )
  audit <- S4Vectors::metadata(observed)$CellCharterR$results$zonation$value$zone_audit
  expect_false(any(audit$retained))
  expect_error(cellcharter_zonation(spe, block = "missing"), "lacks requested")
  expect_error(
    cellcharter_zonation(spe, block = "patient", zone_programs = c("Pre-EMT", "missing")),
    "lack zone"
  )
})

test_that("zonation plot supports group and within-group strata comparisons", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_zonation(
    make_zonation_spe(), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 3L
  )
  group_plot <- cellcharter_plot_zonation(
    spe, group = "risk", group_levels = c("High", "Low"),
    show_individual = FALSE
  )
  expect_s3_class(group_plot, "ggplot")
  expect_equal(nrow(attr(group_plot, "cellcharter_zonation_tests")), 6L)
  expect_identical(
    levels(attr(group_plot, "cellcharter_zonation_data")$program), c("A", "B")
  )

  strata_plot <- suppressWarnings(cellcharter_plot_zonation(
    spe, compare = "strata", group = "risk", within_group = "High",
    stratify_by = "site", strata = c("AR", "BM"),
    show_individual = TRUE
  ))
  expect_s3_class(strata_plot, "ggplot")
  strata_data <- attr(strata_plot, "cellcharter_zonation_data")
  expect_true(all(strata_data$risk == "High"))
  expect_identical(levels(strata_data$curve_group), c("AR", "BM"))
  expect_equal(nrow(attr(strata_plot, "cellcharter_zonation_tests")), 6L)
})

test_that("zonation plot can facet by group while styling programs by linetype", {
  skip_if_not_installed("ggplot2")
  public_formals <- names(formals(cellcharter_plot_zonation))
  legacy_tail <- match("colours", public_formals):match("title", public_formals)
  expect_identical(
    public_formals[legacy_tail],
    c("colours", "facet", "facet_ncol", "title")
  )
  spe <- cellcharter_zonation(
    make_zonation_spe(), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 3L
  )
  plot <- cellcharter_plot_zonation(
    spe,
    programs = c("A", "B"),
    compare = "groups",
    group = "risk",
    group_levels = c("Low", "High"),
    facet_by = "group",
    facet_ncol = 1L,
    linetypes = c(A = "solid", B = "dashed"),
    show_individual = FALSE,
    show_significance = TRUE
  )

  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)
  expect_identical(
    as.character(built$layout$layout$curve_group),
    c("Low", "High")
  )
  line_layer <- built$data[[2L]]
  expect_setequal(unique(line_layer$linetype), c("solid", "dashed"))
  text_layer <- built$data[[5L]]
  expect_true(all(text_layer$PANEL == 1L))

  expect_error(
    cellcharter_plot_zonation(spe, compare = "groups", facet_by = "group"),
    "requires `group`"
  )
  expect_error(
    cellcharter_plot_zonation(
      spe, compare = "strata", group = "risk", stratify_by = "site",
      facet_by = "group"
    ),
    "requires `compare = 'groups'` or `compare = 'none'`",
    fixed = TRUE
  )
  expect_error(
    cellcharter_plot_zonation(
      spe, programs = c("A", "B"), compare = "groups", group = "risk",
      facet_by = "group", linetypes = c(A = "solid")
    ),
    "named for every displayed program"
  )
})

test_that("zonation plot overlays and compares programs within each zone", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_zonation(
    make_zonation_spe(), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 3L
  )
  programs <- c("Pre-EMT", "Early-EMT", "Late-EMT")
  plot <- cellcharter_plot_zonation(
    spe, programs = programs, compare = "programs", facet = FALSE,
    show_individual = FALSE
  )
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$facet, "FacetNull")
  tests <- attr(plot, "cellcharter_zonation_tests")
  expect_equal(nrow(tests), 3L)
  expect_identical(as.character(tests$zone), programs)
  expect_true(all(tests$n_blocks == 8L))
  expect_true(all(is.finite(tests$q_value)))
  expect_error(
    cellcharter_plot_zonation(
      spe, programs = programs, compare = "programs", facet = TRUE
    ),
    "requires `facet_by = 'none'`"
  )
})

test_that("stored zonation results reject stale experiments", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_zonation(
    make_zonation_spe(), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 2L
  )
  stale <- spe[, rev(seq_len(ncol(spe)))]
  expect_error(cellcharter_plot_zonation(stale, group = "risk"), "stale")
})

test_that("missing comparison levels produce an explicit untested result", {
  skip_if_not_installed("ggplot2")
  spe <- cellcharter_zonation(
    make_zonation_spe(), block = "patient", min_spots_per_zone = 2L,
    sample_size = 2L, n_iterations = 2L
  )
  record <- S4Vectors::metadata(spe)$CellCharterR$results$zonation
  remove <- record$value$block_results$program == "A" &
    record$value$block_results$zone == "Late-EMT" &
    record$value$block_results$block %in% paste0("p", 5:8)
  record$value$block_results$enrichment_z[remove] <- NA_real_
  S4Vectors::metadata(spe)$CellCharterR$results$zonation <- record
  plot <- cellcharter_plot_zonation(
    spe, group = "risk", group_levels = c("High", "Low"),
    show_individual = FALSE
  )
  tests <- attr(plot, "cellcharter_zonation_tests")
  row <- tests[tests$program == "A" & tests$zone == "Late-EMT", ]
  expect_true(is.na(row$p_value))
  expect_true(is.na(row$q_value))
  expect_identical(row$reason, "insufficient comparison groups")
  expect_identical(row$label, "not tested")
  expect_identical(row$n_Low, 0L)
})
