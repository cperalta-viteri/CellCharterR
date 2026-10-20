.cc_zonation_scores <- function(sce, scores, assay_name) {
  .cc_validate_sce_identity(sce, reduced_dims = NULL)
  if (!is.character(scores) || length(scores) != 1L || is.na(scores) || !nzchar(scores)) {
    stop("`scores` must name one altExp or reducedDim.", call. = FALSE)
  }
  if (scores %in% SingleCellExperiment::altExpNames(sce)) {
    alternative <- SingleCellExperiment::altExp(sce, scores)
    if (!assay_name %in% SummarizedExperiment::assayNames(alternative)) {
      stop("altExp `", scores, "` lacks assay `", assay_name, "`.", call. = FALSE)
    }
    value <- SummarizedExperiment::assay(alternative, assay_name)
    if (ncol(value) != ncol(sce)) {
      stop("AUCell altExp assay must have one column per spot.", call. = FALSE)
    }
    programs <- rownames(value)
    orientation <- "program_by_spot"
    subset <- function(rows, selected) {
      t(as.matrix(value[selected, rows, drop = FALSE]))
    }
    source <- paste0("altExp:", scores, "/assay:", assay_name)
  } else if (scores %in% SingleCellExperiment::reducedDimNames(sce)) {
    value <- SingleCellExperiment::reducedDim(sce, scores, withDimnames = TRUE)
    if (nrow(value) != ncol(sce)) {
      stop("AUCell reducedDim must have one row per spot.", call. = FALSE)
    }
    programs <- colnames(value)
    orientation <- "spot_by_program"
    subset <- function(rows, selected) {
      as.matrix(value[rows, selected, drop = FALSE])
    }
    source <- paste0("reducedDim:", scores)
  } else {
    stop("SCE lacks AUCell score source `", scores, "`.", call. = FALSE)
  }
  numeric_value <- is.numeric(value) ||
    (inherits(value, "DelayedArray") &&
       DelayedArray::type(value) %in% c("integer", "double"))
  if (!numeric_value || is.null(programs) || anyNA(programs) ||
      any(!nzchar(programs)) || anyDuplicated(programs)) {
    stop("AUCell scores must be numeric with unique, non-empty program names.",
         call. = FALSE)
  }
  list(programs = programs, subset = subset, source = source,
       orientation = orientation)
}

.cc_zonation_score_sha256 <- function(sce, scores, assay_name) {
  .cc_validate_sce_identity(sce, reduced_dims = NULL)
  if (scores %in% SingleCellExperiment::altExpNames(sce)) {
    alternative <- SingleCellExperiment::altExp(sce, scores)
    if (!assay_name %in% SummarizedExperiment::assayNames(alternative)) {
      stop("altExp `", scores, "` lacks assay `", assay_name, "`.", call. = FALSE)
    }
    value <- SummarizedExperiment::assay(alternative, assay_name)
  } else if (scores %in% SingleCellExperiment::reducedDimNames(sce)) {
    value <- SingleCellExperiment::reducedDim(sce, scores, withDimnames = TRUE)
  } else {
    stop("SCE lacks AUCell score source `", scores, "`.", call. = FALSE)
  }
  .autok_representation_sha256(value)
}

.cc_zonation_block_values <- function(sce, block) {
  if (is.null(block)) {
    return(rep.int(".all", ncol(sce)))
  }
  if (!is.character(block) || length(block) != 1L || is.na(block) ||
      !nzchar(block)) {
    stop("`block` must be `NULL` or one non-empty colData column name.",
         call. = FALSE)
  }
  values <- as.character(.cc_sce_coldata(sce, block, "block column"))
  if (anyNA(values) || any(!nzchar(trimws(values)))) {
    stop("`block` contains missing or blank values.", call. = FALSE)
  }
  values
}

.cc_zonation_positive_composition <- function(x, threshold) {
  x[x <= 0] <- 0
  totals <- rowSums(x)
  x <- x / totals
  x[!is.finite(x)] <- 0
  x[x < threshold] <- 0
  totals <- rowSums(x)
  x <- x / totals
  x[!is.finite(x)] <- 0
  x
}

.cc_zonation_rarefy <- function(totals, programs, sample_size,
                                n_iterations, seed) {
  answer <- matrix(
    NA_real_, nrow = length(programs), ncol = length(totals),
    dimnames = list(programs, names(totals))
  )
  if (!length(totals)) return(list(value = answer, depth = NA_integer_))
  depth <- min(sample_size, as.integer(min(vapply(totals, sum, numeric(1)))))
  if (!is.finite(depth) || depth <= 0L) {
    return(list(value = answer, depth = as.integer(depth)))
  }
  pools <- lapply(totals, function(counts) {
    rep.int(programs, times = as.integer(counts))
  })
  trials <- withr::with_seed(seed, lapply(seq_len(n_iterations), function(iteration) {
    lapply(pools, function(pool) {
      tabulate(match(sample(pool, depth, replace = FALSE), programs),
               nbins = length(programs))
    })
  }))
  for (position in seq_along(totals)) {
    values <- vapply(trials, `[[`, integer(length(programs)), position)
    answer[, position] <- rowMeans(values)
  }
  list(value = answer, depth = depth)
}

#' Assign a block-wise dominant program
#'
#' Standardizes every selected score program independently within each
#' patient/block, or globally when `block = NULL`, and assigns each spot to the
#' program with the largest scaled score. Processing occurs one block at a time,
#' so a disk-backed score matrix is never realized in full.
#'
#' @param sce A `SpatialExperiment` or compatible `SingleCellExperiment`.
#' @param scores Name of an AUCell `altExp` (programs by spots) or `reducedDim`
#'   (spots by programs).
#' @param assay_name Assay used when `scores` identifies an `altExp`.
#' @param block Optional patient/subject column in `colData(sce)`. `NULL`
#'   treats the complete object as one block and performs global scaling.
#' @param programs Programs allowed to compete for the assignment. `NULL` uses
#'   every program in `scores`.
#' @param output Name of the `colData` column receiving the assignment.
#' @param require_positive If `TRUE`, leave a spot unassigned when its largest
#'   block-standardized score is not positive. This matches zonation assignment.
#' @param overwrite Whether an existing `output` column may be replaced.
#' @param progress Whether to display a block progress bar.
#' @return The updated experiment with a factor in `colData(sce)[[output]]` and
#'   provenance recorded under `metadata(sce)$CellCharterR$stages`.
#' @export
cellcharter_assign_dominant_program <- function(
    sce, scores = "AUCell", assay_name = "score", block = NULL,
    programs = NULL, output = "dominant_program", require_positive = FALSE,
    overwrite = FALSE, progress = interactive()) {
  block_values <- .cc_zonation_block_values(sce, block)
  source <- .cc_zonation_scores(sce, scores, assay_name)
  if (is.null(programs)) programs <- source$programs
  if (!is.character(programs) || !length(programs) || anyNA(programs) ||
      any(!nzchar(programs)) || anyDuplicated(programs) ||
      !all(programs %in% source$programs)) {
    stop("`programs` must be unique names present in the score source.",
         call. = FALSE)
  }
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(output)) {
    stop("`output` must be one non-empty column name.", call. = FALSE)
  }
  if (!is.logical(require_positive) || length(require_positive) != 1L ||
      is.na(require_positive) || !is.logical(overwrite) ||
      length(overwrite) != 1L || is.na(overwrite) || !is.logical(progress) ||
      length(progress) != 1L || is.na(progress)) {
    stop("`require_positive`, `overwrite`, and `progress` must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (output %in% colnames(SummarizedExperiment::colData(sce)) && !overwrite) {
    stop("colData already contains `", output,
         "`; set `overwrite = TRUE` to replace it.", call. = FALSE)
  }

  blocks <- unique(block_values)
  assigned <- rep.int(NA_character_, ncol(sce))
  bar <- if (progress) utils::txtProgressBar(0L, length(blocks), style = 3L) else NULL
  if (!is.null(bar)) on.exit(close(bar), add = TRUE)
  for (i in seq_along(blocks)) {
    rows <- which(block_values == blocks[[i]])
    raw <- source$subset(rows, programs)
    if (any(!is.finite(raw))) {
      stop("Scores contain non-finite values in block `", blocks[[i]], "`.",
           call. = FALSE)
    }
    centers <- colMeans(raw)
    spreads <- apply(raw, 2L, stats::sd)
    spreads[!is.finite(spreads) | spreads == 0] <- 1
    scaled <- sweep(sweep(raw, 2L, centers, "-"), 2L, spreads, "/")
    winner <- max.col(scaled, ties.method = "first")
    if (require_positive) {
      maximum <- scaled[cbind(seq_len(nrow(scaled)), winner)]
      winner[maximum <= 0] <- NA_integer_
    }
    assigned[rows] <- ifelse(is.na(winner), NA_character_, programs[winner])
    if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
  }

  SummarizedExperiment::colData(sce)[[output]] <- factor(
    assigned, levels = programs
  )
  .cc_sce_record_stage(
    sce, paste0("dominant_program:", output),
    list(scores = scores, score_source = source$source,
         score_orientation = source$orientation, assay_name = assay_name,
         block = block, programs = programs, output = output,
         require_positive = require_positive)
  )
}

#' Compute authors-style score-defined zonation
#'
#' Adapts the equal-depth zonal enrichment mechanics used by Simkin et al. to
#' ordered, mutually exclusive zones defined from AUCell programs. Scores are
#' standardized separately within each patient/block, or globally when
#' `block = NULL`. The strongest positive
#' zone program assigns each spot; spots without a positive zone score remain
#' unassigned. This analysis does not use spatial distance or a graph. Stored
#' score, block, and parameter fingerprints prevent plotting stale results.
#'
#' @param sce A `SpatialExperiment` or compatible `SingleCellExperiment`.
#' @param scores Name of an AUCell `altExp` (programs by spots) or legacy
#'   `reducedDim` (spots by programs).
#' @param assay_name Assay name used when `scores` identifies an `altExp`.
#' @param block Optional patient/subject column in `colData(sce)`. `NULL`
#'   treats the complete object as one block and supports descriptive plots,
#'   but not block-level group inference.
#' @param zone_programs Ordered programs defining mutually exclusive zones.
#' @param composition_programs Response programs quantified across the ordered
#'   zones. They need not include `zone_programs`. `NULL` uses every available
#'   program and treats non-zone programs as the default plotting responses.
#' @param result_name Storage name under
#'   `metadata(sce)$CellCharterR$results`.
#' @param zone_col Optional `colData` column receiving spot-level zone labels.
#' @param min_spots_per_zone Zones with this many or fewer spots are excluded.
#' @param sample_size Equal rarefaction depth used when available.
#' @param n_iterations Number of rarefaction iterations.
#' @param score_threshold Minimum within-spot composition retained after
#'   positive-score normalization.
#' @param seed Deterministic rarefaction seed.
#' @return The updated experiment with compact block-level results stored in
#'   `metadata` and, when requested, spot-level zones in `colData`.
#' @export
cellcharter_zonation <- function(
    sce, scores = "AUCell", assay_name = "score", block = NULL,
    zone_programs = c("Pre-EMT", "Early-EMT", "Late-EMT"),
    composition_programs = NULL, result_name = "zonation",
    zone_col = "cellcharter_zone", min_spots_per_zone = 70L,
    sample_size = 200L, n_iterations = 100L, score_threshold = 0.1,
    seed = 5410L) {
  block_values <- .cc_zonation_block_values(sce, block)
  source <- .cc_zonation_scores(sce, scores, assay_name)
  if (!is.character(zone_programs) || length(zone_programs) < 2L ||
      anyNA(zone_programs) || any(!nzchar(zone_programs)) ||
      anyDuplicated(zone_programs)) {
    stop("`zone_programs` must contain at least two unique program names.",
         call. = FALSE)
  }
  if (!all(zone_programs %in% source$programs)) {
    stop("AUCell scores lack zone program(s): ",
         paste(setdiff(zone_programs, source$programs), collapse = ", "),
         call. = FALSE)
  }
  composition_was_default <- is.null(composition_programs)
  if (composition_was_default) composition_programs <- source$programs
  if (!is.character(composition_programs) || !length(composition_programs) ||
      anyNA(composition_programs) || any(!nzchar(composition_programs)) ||
      anyDuplicated(composition_programs) ||
      !all(composition_programs %in% source$programs)) {
    stop("`composition_programs` must be unique names present in the AUCell source.",
         call. = FALSE)
  }
  integer_args <- list(min_spots_per_zone = min_spots_per_zone,
                       sample_size = sample_size, n_iterations = n_iterations,
                       seed = seed)
  for (name in names(integer_args)) {
    value <- integer_args[[name]]
    minimum <- if (name == "min_spots_per_zone") 0 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < minimum || value != as.integer(value)) {
      stop("`", name, "` must be one integer >= ", minimum, ".", call. = FALSE)
    }
  }
  if (!is.numeric(score_threshold) || length(score_threshold) != 1L ||
      !is.finite(score_threshold) || score_threshold < 0 || score_threshold >= 1) {
    stop("`score_threshold` must be in [0, 1).", call. = FALSE)
  }
  if (!is.character(result_name) || length(result_name) != 1L ||
      is.na(result_name) || !nzchar(result_name)) {
    stop("`result_name` must be one non-empty string.", call. = FALSE)
  }
  if (!is.null(zone_col) && (!is.character(zone_col) || length(zone_col) != 1L ||
      is.na(zone_col) || !nzchar(zone_col))) {
    stop("`zone_col` must be `NULL` or one non-empty string.", call. = FALSE)
  }

  block_levels <- unique(block_values)
  assigned <- rep.int(NA_character_, ncol(sce))
  output <- vector("list", length(block_levels))
  zone_audit <- vector("list", length(block_levels))
  for (block_index in seq_along(block_levels)) {
    current_block <- block_levels[[block_index]]
    rows <- which(block_values == current_block)
    analysis_programs <- unique(c(zone_programs, composition_programs))
    raw <- source$subset(rows, analysis_programs)
    if (any(!is.finite(raw))) {
      stop("AUCell scores contain non-finite values in block `", current_block,
           "`.", call. = FALSE)
    }
    centers <- colMeans(raw)
    spreads <- apply(raw, 2L, stats::sd)
    spreads[!is.finite(spreads) | spreads == 0] <- 1
    scaled <- sweep(raw, 2L, centers, FUN = "-")
    scaled <- sweep(scaled, 2L, spreads, FUN = "/")

    zone_scores <- scaled[, zone_programs, drop = FALSE]
    zone_index <- max.col(zone_scores, ties.method = "first")
    maxima <- zone_scores[cbind(seq_len(nrow(zone_scores)), zone_index)]
    zone_index[maxima <= 0] <- NA_integer_
    assigned[rows] <- ifelse(is.na(zone_index), NA_character_,
                             zone_programs[zone_index])

    composition <- .cc_zonation_positive_composition(
      scaled[, composition_programs, drop = FALSE], score_threshold
    )
    zone_sizes <- vapply(seq_along(zone_programs), function(zone) {
      sum(zone_index == zone, na.rm = TRUE)
    }, integer(1))
    retained <- which(zone_sizes > as.integer(min_spots_per_zone))
    totals <- lapply(retained, function(zone) {
      in_zone <- !is.na(zone_index) & zone_index == zone
      round(colSums(composition[in_zone, , drop = FALSE]))
    })
    names(totals) <- zone_programs[retained]
    rarefied <- .cc_zonation_rarefy(
      totals, composition_programs, as.integer(sample_size),
      as.integer(n_iterations), as.integer(seed)
    )
    mean_sampled <- matrix(
      NA_real_, nrow = length(composition_programs), ncol = length(zone_programs),
      dimnames = list(composition_programs, zone_programs)
    )
    if (length(retained)) mean_sampled[, retained] <- rarefied$value

    long <- expand.grid(
      program = composition_programs,
      zone = zone_programs,
      stringsAsFactors = FALSE
    )
    long$mean_sampled_count <- mapply(
      function(program, zone) mean_sampled[program, zone],
      long$program, long$zone
    )
    long$enrichment_z <- NA_real_
    for (program in composition_programs) {
      positions <- which(long$program == program & is.finite(long$mean_sampled_count))
      values <- long$mean_sampled_count[positions]
      if (length(values) >= 2L && stats::sd(values) > 0) {
        long$enrichment_z[positions] <- as.numeric(scale(values))
      }
    }
    long$block <- current_block
    long$zone_spots <- zone_sizes[match(long$zone, zone_programs)]
    long$retained_zone <- match(long$zone, zone_programs) %in% retained
    long$sampled_state_units <- ifelse(long$retained_zone, rarefied$depth, NA_integer_)
    output[[block_index]] <- long
    zone_audit[[block_index]] <- data.frame(
      block = current_block, zone = zone_programs, spots = zone_sizes,
      retained = seq_along(zone_programs) %in% retained,
      sampled_state_units = ifelse(seq_along(zone_programs) %in% retained,
                                   rarefied$depth, NA_integer_),
      stringsAsFactors = FALSE
    )
  }
  result <- list(
    schema = "cellcharter_zonation_v2",
    block_mode = if (is.null(block)) "pooled" else "blocked",
    block_results = do.call(rbind, output),
    zone_audit = do.call(rbind, zone_audit),
    zone_programs = zone_programs,
    composition_programs = composition_programs,
    response_programs = if (composition_was_default) {
      setdiff(composition_programs, zone_programs)
    } else {
      composition_programs
    },
    reference = list(
      method = "Simkin et al. authors-style equal-depth zonal enrichment adaptation",
      commit = "d53c7332270fb2f41a9ca14647ae390314f24ddd"
    )
  )
  if (!is.null(zone_col)) {
    SummarizedExperiment::colData(sce)[[zone_col]] <- factor(
      assigned, levels = zone_programs, ordered = TRUE
    )
  }
  .cc_sce_store_result(
    sce, result_name, result,
    list(
      scores = scores, score_source = source$source,
      score_orientation = source$orientation, assay_name = assay_name,
      block = block, zone_programs = zone_programs,
      composition_programs = composition_programs, zone_col = zone_col,
      min_spots_per_zone = as.integer(min_spots_per_zone),
      sample_size = as.integer(sample_size),
      n_iterations = as.integer(n_iterations),
      score_threshold = score_threshold, seed = as.integer(seed)
    ),
    c(list(
      scores = list(kind = "zonation_scores", scores = scores,
                    assay_name = assay_name)
    ), if (is.null(block)) list() else list(
      block = list(kind = "colData", name = block)
    ))
  )
}

.cc_zonation_result <- function(x, result_name) {
  if (!inherits(x, "SingleCellExperiment")) {
    stop("`x` must be a SingleCellExperiment or SpatialExperiment.", call. = FALSE)
  }
  record <- .cc_sce_require_result(
    x, result_name, context = "Stored zonation result"
  )
  if (!record$value$schema %in% c(
      "cellcharter_zonation_v1", "cellcharter_zonation_v2")) {
    stop("SCE lacks stored zonation result `", result_name,
         "`; run `cellcharter_zonation()` first.", call. = FALSE)
  }
  record
}

.cc_zonation_block_annotations <- function(sce, block_col, columns) {
  if (is.null(block_col)) {
    if (length(columns)) {
      stop(
        "Group or stratum comparisons require zonation computed with a `block` column.",
        call. = FALSE
      )
    }
    return(data.frame(block = ".all", stringsAsFactors = FALSE))
  }
  block <- as.character(.cc_sce_coldata(sce, block_col, "stored block column"))
  answer <- data.frame(block = unique(block), stringsAsFactors = FALSE)
  for (column in columns) {
    values <- as.character(.cc_sce_coldata(sce, column, paste0("annotation column `", column, "`")))
    split_values <- split(values, block)
    conflicts <- names(split_values)[vapply(split_values, function(z) length(unique(z)) != 1L,
                                            logical(1))]
    if (length(conflicts)) {
      stop("Annotation `", column, "` is not constant within block(s): ",
           paste(utils::head(conflicts, 5L), collapse = ", "), call. = FALSE)
    }
    mapping <- vapply(split_values, `[[`, character(1), 1L)
    answer[[column]] <- unname(mapping[answer$block])
  }
  answer
}

#' Plot score-defined zonation and patient-level comparisons
#'
#' Plots compact block-level results from [cellcharter_zonation()]. It can
#' compare two groups (for example High versus Low risk), compare metadata
#' strata within one group (for example oral sites among High-risk patients),
#' compare the displayed programs within each zone, or display descriptive
#' curves without inference.
#'
#' @param x An experiment containing a stored zonation result.
#' @param result_name Stored result name.
#' @param programs Programs to display. `NULL` displays every non-zone program.
#' @param compare One of `"groups"`, `"strata"`, `"programs"`, or `"none"`.
#'   `"programs"` performs a block-paired Friedman test among the displayed
#'   programs separately within each zone and requires `facet = FALSE`.
#' @param group `colData` column defining groups.
#' @param group_levels Optional ordered group levels to retain.
#' @param within_group For `compare = "strata"`, one group value to retain.
#' @param stratify_by Optional `colData` column defining strata.
#' @param strata Optional ordered strata to retain.
#' @param show_individual Draw faint block/patient trajectories.
#' @param show_significance Add adjusted q-values to the plot.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param colours Optional named colour vector for comparison levels. With
#'   `facet = FALSE` and `compare = "none"` or `"programs"`, name colours by
#'   program instead.
#' @param facet Whether to draw one facet per program. Set to `FALSE` to overlay
#'   all requested programs in one panel.
#' @param facet_ncol Number of facet columns; ignored when `facet = FALSE`.
#' @param title Optional title.
#' @return A `ggplot` object. Block-level data, summaries, and tests are
#'   attached as `cellcharter_zonation_data`, `cellcharter_zonation_summary`,
#'   and `cellcharter_zonation_tests` attributes.
#' @export
cellcharter_plot_zonation <- function(
    x, result_name = "zonation", programs = NULL,
    compare = c("groups", "strata", "programs", "none"), group = NULL,
    group_levels = NULL, within_group = NULL,
    stratify_by = NULL, strata = NULL, show_individual = TRUE,
    show_significance = TRUE, p_adjust_method = "BH", colours = NULL,
    facet = TRUE, facet_ncol = 3L, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install `ggplot2` to plot zonation.", call. = FALSE)
  }
  compare <- match.arg(compare)
  if (!is.logical(facet) || length(facet) != 1L || is.na(facet)) {
    stop("`facet` must be TRUE or FALSE.", call. = FALSE)
  }
  record <- .cc_zonation_result(x, result_name)
  result <- record$value
  parameters <- record$parameters
  data <- result$block_results
  if (is.null(programs)) programs <- result$response_programs
  if (!is.character(programs) || !length(programs) ||
      !all(programs %in% result$composition_programs)) {
    stop("`programs` must name programs stored in the zonation result.", call. = FALSE)
  }
  if (compare == "programs" && facet) {
    stop("`compare = 'programs'` requires `facet = FALSE` so the programs share one panel.",
         call. = FALSE)
  }
  if (compare == "programs" && length(programs) < 2L) {
    stop("`compare = 'programs'` requires at least two displayed programs.",
         call. = FALSE)
  }
  columns <- unique(c(group, stratify_by))
  columns <- columns[!is.na(columns) & nzchar(columns)]
  annotations <- .cc_zonation_block_annotations(x, parameters$block, columns)
  data <- merge(data, annotations, by = "block", all.x = TRUE, sort = FALSE)
  data <- data[data$program %in% programs & is.finite(data$enrichment_z), , drop = FALSE]
  if (!is.null(group_levels)) {
    if (is.null(group)) stop("`group_levels` requires `group`.", call. = FALSE)
    data <- data[data[[group]] %in% group_levels, , drop = FALSE]
  }
  if (!is.null(strata)) {
    if (is.null(stratify_by)) stop("`strata` requires `stratify_by`.", call. = FALSE)
    data <- data[data[[stratify_by]] %in% strata, , drop = FALSE]
  }
  if (!is.null(within_group)) {
    if (is.null(group) || length(within_group) != 1L) {
      stop("`within_group` requires `group` and one group value.", call. = FALSE)
    }
    data <- data[data[[group]] == within_group, , drop = FALSE]
  }
  if (!nrow(data)) stop("No finite zonation values remain after filtering.", call. = FALSE)

  if (compare == "groups") {
    if (is.null(group)) stop("`compare = 'groups'` requires `group`.", call. = FALSE)
    data$curve_group <- as.character(data[[group]])
    levels_used <- if (is.null(group_levels)) unique(data$curve_group) else group_levels
    levels_used <- levels_used[levels_used %in% data$curve_group]
    if (length(levels_used) != 2L) {
      stop("Group comparison requires exactly two retained group levels.", call. = FALSE)
    }
    test_kind <- "Wilcoxon rank-sum"
  } else if (compare == "strata") {
    if (is.null(stratify_by)) stop("`compare = 'strata'` requires `stratify_by`.", call. = FALSE)
    data$curve_group <- as.character(data[[stratify_by]])
    levels_used <- if (is.null(strata)) unique(data$curve_group) else strata
    levels_used <- levels_used[levels_used %in% data$curve_group]
    if (length(levels_used) < 2L) {
      stop("Strata comparison requires at least two retained strata.", call. = FALSE)
    }
    test_kind <- "Kruskal-Wallis"
  } else if (compare == "programs") {
    data$curve_group <- "All blocks"
    levels_used <- "All blocks"
    test_kind <- "Block-paired Friedman"
  } else {
    if (!is.null(group)) {
      data$curve_group <- as.character(data[[group]])
      levels_used <- if (is.null(group_levels)) unique(data$curve_group) else group_levels
    } else {
      data$curve_group <- "All blocks"
      levels_used <- "All blocks"
    }
    test_kind <- "none"
  }
  data <- data[data$curve_group %in% levels_used, , drop = FALSE]
  data$curve_group <- factor(data$curve_group, levels = levels_used)
  data$zone <- factor(data$zone, levels = result$zone_programs, ordered = TRUE)
  data$program <- factor(data$program, levels = programs)

  summary_keys <- interaction(data$program, data$curve_group, data$zone, drop = TRUE)
  summary <- do.call(rbind, lapply(split(data, summary_keys), function(z) {
    data.frame(
      program = z$program[[1]], curve_group = z$curve_group[[1]], zone = z$zone[[1]],
      mean = mean(z$enrichment_z),
      se = if (nrow(z) > 1L) stats::sd(z$enrichment_z) / sqrt(nrow(z)) else NA_real_,
      blocks = length(unique(z$block)), stringsAsFactors = FALSE
    )
  }))
  summary$curve_group <- factor(summary$curve_group, levels = levels_used)
  summary$zone <- factor(summary$zone, levels = result$zone_programs, ordered = TRUE)
  summary$program <- factor(summary$program, levels = programs)

  overlay_by_program <- !facet && compare %in% c("none", "programs")
  colour_levels <- if (overlay_by_program) programs else levels_used
  data$plot_colour <- factor(
    if (overlay_by_program) as.character(data$program) else as.character(data$curve_group),
    levels = colour_levels
  )
  summary$plot_colour <- factor(
    if (overlay_by_program) as.character(summary$program) else as.character(summary$curve_group),
    levels = colour_levels
  )
  data$plot_group <- interaction(
    data$block, data$curve_group, data$program, drop = TRUE
  )
  summary$plot_group <- interaction(
    summary$curve_group, summary$program, drop = TRUE
  )

  tests <- NULL
  if (compare == "programs") {
    tests <- do.call(rbind, lapply(split(data, data$zone), function(z) {
      # Retain only blocks with one finite value for every displayed program;
      # these are the repeated-measure units in the Friedman test.
      availability <- table(z$block, factor(z$program, levels = programs))
      complete_blocks <- rownames(availability)[rowSums(availability == 1L) == length(programs)]
      paired <- z[z$block %in% complete_blocks, , drop = FALSE]
      if (length(complete_blocks) < 3L) {
        warning("Zone `", as.character(z$zone[[1L]]),
                "` has fewer than three complete blocks for the program comparison; ",
                "the zone was not tested.", call. = FALSE)
        statistic <- NA_real_
        p_value <- NA_real_
        reason <- "fewer than three complete blocks"
      } else {
        test <- stats::friedman.test(
          enrichment_z ~ program | block,
          data = paired
        )
        statistic <- unname(test$statistic)
        p_value <- test$p.value
        reason <- NA_character_
      }
      data.frame(
        zone = as.character(z$zone[[1L]]), statistic = statistic,
        p_value = p_value, reason = reason,
        n_blocks = length(complete_blocks), stringsAsFactors = FALSE
      )
    }))
  } else if (compare != "none") {
    test_keys <- interaction(data$program, data$zone, drop = TRUE)
    tests <- do.call(rbind, lapply(split(data, test_keys), function(z) {
      counts <- table(factor(z$curve_group, levels = levels_used))
      present <- counts > 0L
      sufficient_groups <- if (compare == "groups") {
        all(present)
      } else {
        sum(present) >= 2L
      }
      if (any(counts[present] < 3L)) {
        warning("At least one comparison level has fewer than three blocks; ",
                "interpret `", as.character(z$program[[1]]), "` / `",
                as.character(z$zone[[1]]), "` cautiously.", call. = FALSE)
      }
      if (!sufficient_groups) {
        statistic <- NA_real_
        p_value <- NA_real_
        reason <- "insufficient comparison groups"
      } else {
        test <- if (compare == "groups") {
          stats::wilcox.test(
            z$enrichment_z[z$curve_group == levels_used[[1L]]],
            z$enrichment_z[z$curve_group == levels_used[[2L]]],
            exact = FALSE
          )
        } else {
          stats::kruskal.test(z$enrichment_z, droplevels(z$curve_group))
        }
        statistic <- unname(test$statistic)
        p_value <- test$p.value
        reason <- NA_character_
      }
      row <- data.frame(
        program = as.character(z$program[[1]]), zone = as.character(z$zone[[1]]),
        statistic = statistic, p_value = p_value, reason = reason,
        stringsAsFactors = FALSE
      )
      for (level in levels_used) row[[paste0("n_", make.names(level))]] <- unname(counts[[level]])
      row
    }))
  }
  if (!is.null(tests)) {
    tests$q_value <- stats::p.adjust(tests$p_value, method = p_adjust_method)
    tests$significance <- cut(
      tests$q_value, c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "ns"), right = TRUE
    )
    tests$label <- sprintf("q=%.2g %s", tests$q_value, tests$significance)
    tests$label[!is.finite(tests$q_value)] <- "not tested"
    tests$zone <- factor(tests$zone, levels = result$zone_programs, ordered = TRUE)
    tests$annotation_y <- max(data$enrichment_z) + 0.15 * diff(range(data$enrichment_z))
    if (compare != "programs") {
      tests$program <- factor(tests$program, levels = programs)
    }
    if (!facet && compare != "programs" && length(programs) > 1L) {
      step <- max(0.12, 0.1 * diff(range(data$enrichment_z)))
      tests$annotation_y <- tests$annotation_y +
        (match(as.character(tests$program), programs) - 1L) * step
      tests$label <- paste0(as.character(tests$program), ": ", tests$label)
    }
  }

  if (is.null(colours)) {
    palette <- c("#7F0000", "#08306B", "#009E73", "#CC79A7", "#6A3D9A",
                 "#E69F00", "#56B4E9", "#000000")
    if (length(colour_levels) > length(palette)) {
      palette <- grDevices::hcl.colors(length(colour_levels), "Dark 3")
    }
    colours <- stats::setNames(palette[seq_along(colour_levels)], colour_levels)
  } else if (is.null(names(colours)) || !all(colour_levels %in% names(colours))) {
    stop(
      "`colours` must be named for every ",
      if (overlay_by_program) "displayed program." else "retained comparison level.",
      call. = FALSE
    )
  }
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = zone, y = enrichment_z, colour = plot_colour)
  ) + ggplot2::geom_hline(yintercept = 0, colour = "grey70", linetype = 2)
  if (isTRUE(show_individual)) {
    plot <- plot + ggplot2::geom_line(
      ggplot2::aes(group = plot_group),
      alpha = 0.15, linewidth = 0.35, na.rm = TRUE
    )
  }
  summary_line_mapping <- if (!facet && !overlay_by_program &&
                              compare != "none" && length(programs) > 1L) {
    ggplot2::aes(
      x = zone, y = mean, colour = plot_colour,
      group = plot_group, linetype = program
    )
  } else {
    ggplot2::aes(
      x = zone, y = mean, colour = plot_colour, group = plot_group
    )
  }
  plot <- plot +
    ggplot2::geom_line(
      data = summary,
      summary_line_mapping,
      inherit.aes = FALSE, linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = summary,
      ggplot2::aes(x = zone, y = mean, colour = plot_colour),
      inherit.aes = FALSE, size = 2.4
    ) +
    ggplot2::geom_errorbar(
      data = summary,
      ggplot2::aes(x = zone, ymin = mean - se, ymax = mean + se,
                   colour = plot_colour),
      inherit.aes = FALSE, width = 0.07, linewidth = 0.5, na.rm = TRUE
    )
  if (isTRUE(show_significance) && !is.null(tests)) {
    plot <- plot + ggplot2::geom_text(
      data = tests,
      ggplot2::aes(x = zone, y = annotation_y, label = label),
      inherit.aes = FALSE, colour = "black", size = 2.8, fontface = "bold"
    )
  }
  if (facet) {
    plot <- plot + ggplot2::facet_wrap(~program, ncol = as.integer(facet_ncol))
  }
  plot <- plot +
    ggplot2::scale_colour_manual(values = colours[colour_levels], drop = FALSE) +
    ggplot2::labs(
      x = "Ordered score-defined zone",
      y = if (is.null(parameters$block)) {
        "Global program enrichment z-score"
      } else {
        "Within-block program enrichment z-score"
      },
      colour = if (overlay_by_program) {
        "Program"
      } else if (compare == "strata") {
        stratify_by
      } else if (!is.null(group)) {
        group
      } else {
        NULL
      },
      linetype = if (!facet && !overlay_by_program &&
                    compare != "none" && length(programs) > 1L) {
        "Program"
      } else {
        NULL
      },
      title = title,
      subtitle = if (compare == "none") NULL else if (compare == "programs") {
        paste0(test_kind, " test within each zone; ", p_adjust_method,
               " adjustment across displayed zones")
      } else {
        paste0(test_kind, "; ", p_adjust_method,
               " adjustment across all displayed program x zone tests")
      },
      caption = paste0(
        if (is.null(parameters$block)) {
          "Authors-style equal-depth zonation; all spots form one descriptive block. "
        } else {
          "Authors-style equal-depth zonation; blocks are independent units. "
        },
        "Zones with <=", parameters$min_spots_per_zone, " spots were excluded."
      )
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(), legend.position = "top",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 9),
      plot.title.position = "plot"
    )
  attr(plot, "cellcharter_zonation_data") <- data
  attr(plot, "cellcharter_zonation_summary") <- summary
  attr(plot, "cellcharter_zonation_tests") <- tests
  plot
}
