.cc_zonation_scores <- function(spe, scores, assay_name) {
  .cc_validate_spe_identity(spe, reduced_dims = NULL)
  if (!is.character(scores) || length(scores) != 1L || is.na(scores) || !nzchar(scores)) {
    stop("`scores` must name one altExp or reducedDim.", call. = FALSE)
  }
  if (scores %in% SingleCellExperiment::altExpNames(spe)) {
    alternative <- SingleCellExperiment::altExp(spe, scores)
    if (!assay_name %in% SummarizedExperiment::assayNames(alternative)) {
      stop("altExp `", scores, "` lacks assay `", assay_name, "`.", call. = FALSE)
    }
    value <- SummarizedExperiment::assay(alternative, assay_name)
    if (ncol(value) != ncol(spe)) {
      stop("AUCell altExp assay must have one column per spot.", call. = FALSE)
    }
    programs <- rownames(value)
    orientation <- "program_by_spot"
    subset <- function(rows, selected) {
      t(as.matrix(value[selected, rows, drop = FALSE]))
    }
    source <- paste0("altExp:", scores, "/assay:", assay_name)
  } else if (scores %in% SingleCellExperiment::reducedDimNames(spe)) {
    value <- SingleCellExperiment::reducedDim(spe, scores, withDimnames = TRUE)
    if (nrow(value) != ncol(spe)) {
      stop("AUCell reducedDim must have one row per spot.", call. = FALSE)
    }
    programs <- colnames(value)
    orientation <- "spot_by_program"
    subset <- function(rows, selected) {
      as.matrix(value[rows, selected, drop = FALSE])
    }
    source <- paste0("reducedDim:", scores)
  } else {
    stop("SpatialExperiment lacks AUCell score source `", scores, "`.", call. = FALSE)
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

.cc_zonation_score_sha256 <- function(spe, scores, assay_name) {
  .cc_validate_spe_identity(spe, reduced_dims = NULL)
  if (scores %in% SingleCellExperiment::altExpNames(spe)) {
    alternative <- SingleCellExperiment::altExp(spe, scores)
    if (!assay_name %in% SummarizedExperiment::assayNames(alternative)) {
      stop("altExp `", scores, "` lacks assay `", assay_name, "`.", call. = FALSE)
    }
    value <- SummarizedExperiment::assay(alternative, assay_name)
  } else if (scores %in% SingleCellExperiment::reducedDimNames(spe)) {
    value <- SingleCellExperiment::reducedDim(spe, scores, withDimnames = TRUE)
  } else {
    stop("SpatialExperiment lacks AUCell score source `", scores, "`.", call. = FALSE)
  }
  .autok_representation_sha256(value)
}

.cc_zonation_block_values <- function(spe, block) {
  if (is.null(block)) {
    return(rep.int(".all", ncol(spe)))
  }
  if (!is.character(block) || length(block) != 1L || is.na(block) ||
      !nzchar(block)) {
    stop("`block` must be `NULL` or one non-empty colData column name.",
         call. = FALSE)
  }
  values <- as.character(.cc_spe_coldata(spe, block, "block column"))
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
#' @param spe A `SpatialExperiment` containing program scores for every spot.
#' @param scores One non-empty string naming either an `altExp` whose rows are
#'   programs and columns are spots, or a `reducedDim` whose rows are spots and
#'   columns are programs. Program names must be unique and non-empty.
#' @param assay_name One non-empty string naming the assay to read when
#'   `scores` identifies an `altExp`; ignored for a `reducedDim` source.
#' @param block `NULL` for a single global standardization, or one `colData(spe)`
#'   column containing non-missing patient, subject, slide, or sample IDs.
#'   Scores are centered and scaled independently within every block; this is
#'   a scaling unit, not a statistical correction for block effects.
#' @param programs Unique program names allowed to compete for assignment, in
#'   tie-breaking and output-factor order. `NULL` uses every program available
#'   from `scores`. Ties are resolved in favor of the first listed program.
#' @param output One non-empty name for the factor written to
#'   `colData(spe)[[output]]`. The stage record is written under
#'   `metadata(spe)$CellCharterR$stages[[paste0("dominant_program:", output)]]`.
#' @param require_positive If `TRUE`, leave a spot unassigned when its largest
#'   block-standardized score is zero or negative. If `FALSE`, every spot is
#'   assigned, including spots for which all standardized scores are negative.
#' @param overwrite A single logical indicating whether an existing
#'   `colData(spe)[[output]]` column may be replaced. Nothing is overwritten
#'   when this is `FALSE`.
#' @param progress A single logical indicating whether to display progress over
#'   blocks. It does not affect the result.
#' @return The updated `SpatialExperiment`. The dominant-program factor is in
#'   `colData(spe)[[output]]`; parameters, spot count, and spot-order identity
#'   are in the corresponding `metadata(spe)$CellCharterR$stages` record. No
#'   score matrix or external file is created.
#' @export
cellcharter_assign_dominant_program <- function(
    spe, scores = "AUCell", assay_name = "score", block = NULL,
    programs = NULL, output = "dominant_program", require_positive = FALSE,
    overwrite = FALSE, progress = interactive()) {
  block_values <- .cc_zonation_block_values(spe, block)
  source <- .cc_zonation_scores(spe, scores, assay_name)
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
  if (output %in% colnames(SummarizedExperiment::colData(spe)) && !overwrite) {
    stop("colData already contains `", output,
         "`; set `overwrite = TRUE` to replace it.", call. = FALSE)
  }

  blocks <- unique(block_values)
  assigned <- rep.int(NA_character_, ncol(spe))
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

  SummarizedExperiment::colData(spe)[[output]] <- factor(
    assigned, levels = programs
  )
  .cc_spe_record_stage(
    spe, paste0("dominant_program:", output),
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
#' @param spe A `SpatialExperiment` containing program scores for every spot.
#' @param scores One non-empty string naming either an `altExp` whose rows are
#'   programs and columns are spots, or a `reducedDim` whose rows are spots and
#'   columns are programs.
#' @param assay_name One non-empty string naming the assay to read when
#'   `scores` identifies an `altExp`; ignored for a `reducedDim` source.
#' @param block `NULL` to treat all spots as one descriptive block, or one
#'   `colData(spe)` column containing complete patient, subject, slide, or
#'   sample IDs. Score standardization, rarefaction, and inferential units are
#'   block-specific when supplied. Group inference in
#'   [cellcharter_plot_zonation()] requires a non-`NULL` block.
#' @param zone_programs At least two unique score-program names in biological
#'   order. Within every block, each program is z-scored and each spot is
#'   assigned to its largest positive zone score; ties follow this vector's
#'   order and spots with no positive maximum remain unassigned.
#' @param composition_programs Unique response-program names quantified across
#'   the ordered zones. They may overlap `zone_programs`. `NULL` quantifies all
#'   available programs and makes non-zone programs the default plot set.
#' @param result_name One non-empty storage key. The complete freshness-checked
#'   record is written to
#'   `metadata(spe)$CellCharterR$results[[result_name]]`; scientific outputs
#'   such as `block_results` and `zone_audit` are in that record's `$value`.
#'   Reusing a key replaces its previous record.
#' @param zone_col `NULL` to avoid adding spot-level labels, or one non-empty
#'   name for the ordered zone factor written to `colData(spe)[[zone_col]]`.
#' @param min_spots_per_zone A non-negative whole number. A block-zone is
#'   retained only when its spot count is strictly greater than this value.
#' @param sample_size A positive whole number giving the target rarefaction
#'   depth. Each block uses the smaller of this value and its smallest retained
#'   zone total.
#' @param n_iterations A positive whole number of rarefaction draws whose
#'   program counts are averaged for every retained block-zone.
#' @param score_threshold A finite number in `[0, 1)`. After negative z-scores
#'   are set to zero and positive scores are normalized within a spot,
#'   proportions below this value are discarded and the remainder renormalized.
#' @param seed A positive whole-number seed for deterministic rarefaction.
#' @return The updated `SpatialExperiment`. The result record is stored in
#'   `metadata(spe)$CellCharterR$results[[result_name]]`; when `zone_col` is not
#'   `NULL`, the ordered spot-level zone factor is also stored in
#'   `colData(spe)[[zone_col]]`. No graph, distance, or external file is made.
#' @export
cellcharter_zonation <- function(
    spe, scores = "AUCell", assay_name = "score", block = NULL,
    zone_programs = c("Pre-EMT", "Early-EMT", "Late-EMT"),
    composition_programs = NULL, result_name = "zonation",
    zone_col = "cellcharter_zone", min_spots_per_zone = 70L,
    sample_size = 200L, n_iterations = 100L, score_threshold = 0.1,
    seed = 5410L) {
  block_values <- .cc_zonation_block_values(spe, block)
  source <- .cc_zonation_scores(spe, scores, assay_name)
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
  assigned <- rep.int(NA_character_, ncol(spe))
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
    SummarizedExperiment::colData(spe)[[zone_col]] <- factor(
      assigned, levels = zone_programs, ordered = TRUE
    )
  }
  .cc_spe_store_result(
    spe, result_name, result,
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
  if (!inherits(x, "SpatialExperiment")) {
    stop("`x` must be a SpatialExperiment.", call. = FALSE)
  }
  record <- .cc_spe_require_result(
    x, result_name, context = "Stored zonation result"
  )
  if (!record$value$schema %in% c(
      "cellcharter_zonation_v1", "cellcharter_zonation_v2")) {
    stop("SpatialExperiment lacks stored zonation result `", result_name,
         "`; run `cellcharter_zonation()` first.", call. = FALSE)
  }
  record
}

.cc_zonation_block_annotations <- function(spe, block_col, columns) {
  if (is.null(block_col)) {
    if (length(columns)) {
      stop(
        "Group or stratum comparisons require zonation computed with a `block` column.",
        call. = FALSE
      )
    }
    return(data.frame(block = ".all", stringsAsFactors = FALSE))
  }
  block <- as.character(.cc_spe_coldata(spe, block_col, "stored block column"))
  answer <- data.frame(block = unique(block), stringsAsFactors = FALSE)
  for (column in columns) {
    values <- as.character(.cc_spe_coldata(spe, column, paste0("annotation column `", column, "`")))
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
#' @param x A `SpatialExperiment` containing a result previously stored by
#'   [cellcharter_zonation()]. Plotting is read-only and does not modify `x`.
#' @param result_name One non-empty key under
#'   `metadata(x)$CellCharterR$results`. The result's score, block, and
#'   parameter fingerprints are validated before plotting.
#' @param programs Unique stored composition-program names to display, in plot
#'   order. `NULL` uses the result's response programs: normally every
#'   quantified non-zone program when `composition_programs` was not specified.
#' @param compare One of `"groups"`, `"strata"`, `"programs"`, or `"none"`.
#'   `"groups"` performs two-group Wilcoxon rank-sum tests for each displayed
#'   program-zone combination; `"strata"` performs Kruskal-Wallis tests across
#'   at least two retained strata; `"programs"` performs a block-paired
#'   Friedman test across displayed programs separately in each zone and
#'   requires `facet_by = "none"` (or the legacy `facet = FALSE`);
#'   `"none"` performs no test.
#' @param group `NULL`, or one `colData(x)` column containing a group label that
#'   is constant within every stored block. It supplies the curves and test
#'   groups for `compare = "groups"`, filters for `within_group`, and can supply
#'   descriptive curves for `compare = "none"`.
#' @param group_levels `NULL` to use retained groups in encounter order, or an
#'   ordered vector of group values to keep. A group comparison requires
#'   exactly two retained levels, and their order controls display order.
#' @param within_group `NULL`, or one value from `group` to retain before a
#'   strata comparison or descriptive plot. Supplying it requires `group`.
#' @param stratify_by `NULL`, or one `colData(x)` column containing a stratum
#'   label that is constant within every block. It is required for
#'   `compare = "strata"`.
#' @param strata `NULL` to use retained strata in encounter order, or an
#'   ordered vector of stratum values to keep. Supplying it requires
#'   `stratify_by`; at least two retained strata are needed for inference.
#' @param show_individual A single logical indicating whether faint block-level
#'   trajectories are drawn behind the mean curves.
#' @param show_significance A single logical indicating whether adjusted
#'   q-values and significance symbols are annotated. Tests are still computed
#'   and attached when this is `FALSE`.
#' @param p_adjust_method One method accepted by [stats::p.adjust()]. For group
#'   and strata comparisons, correction is over all displayed program-by-zone
#'   tests together; for program comparisons, it is over displayed zones.
#' @param colours Optional named colour vector for comparison levels. With
#'   `facet_by = "none"` (or the legacy `facet = FALSE`) and
#'   `compare = "none"` or `"programs"`, name colours by program instead.
#'   Every displayed colour level must be represented.
#' @param linetypes `NULL` to use ggplot2's default program linetypes, or a
#'   named vector assigning a linetype to every displayed program (for example,
#'   `c("Pre-EMT" = "solid", "Early-EMT" = "dashed")`). It is used when
#'   multiple programs share a panel, including `facet_by = "group"`.
#' @param facet Backward-compatible logical faceting switch. `TRUE` draws one
#'   facet per program and `FALSE` draws no facets. When `facet_by` is supplied,
#'   `facet_by` takes precedence.
#' @param facet_by `NULL` to retain the behavior selected by `facet`, or one of
#'   `"program"`, `"group"`, or `"none"`. `"group"` draws one facet per
#'   retained value of `group` while distinguishing programs by linetype; it
#'   requires `group` and `compare = "groups"` or `"none"`. Cross-group
#'   significance labels are drawn once in the first group facet. `"program"`
#'   draws one facet per program, and `"none"` overlays the requested curves.
#' @param facet_ncol A positive whole number of facet columns; ignored when
#'   faceting is disabled.
#' @param title `NULL` for no title, or a character value used as the plot title.
#' @return A read-only `ggplot` object. The plotted block-level observations,
#'   mean/standard-error summaries, and inferential table (including raw
#'   `p_value`, adjusted `q_value`, sample counts, and skipped-test reasons) are
#'   attached as `cellcharter_zonation_data`, `cellcharter_zonation_summary`,
#'   and `cellcharter_zonation_tests` attributes.
#' @export
cellcharter_plot_zonation <- function(
    x, result_name = "zonation", programs = NULL,
    compare = c("groups", "strata", "programs", "none"), group = NULL,
    group_levels = NULL, within_group = NULL,
    stratify_by = NULL, strata = NULL, show_individual = TRUE,
    show_significance = TRUE, p_adjust_method = "BH", colours = NULL,
    facet = TRUE, facet_ncol = 3L, title = NULL,
    facet_by = NULL, linetypes = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install `ggplot2` to plot zonation.", call. = FALSE)
  }
  compare <- match.arg(compare)
  if (!is.logical(facet) || length(facet) != 1L || is.na(facet)) {
    stop("`facet` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(facet_by)) {
    facet_by <- if (facet) "program" else "none"
  } else {
    facet_by <- match.arg(facet_by, c("program", "group", "none"))
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
  if (compare == "programs" && facet_by != "none") {
    stop("`compare = 'programs'` requires `facet_by = 'none'` so the programs share one panel.",
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
  if (facet_by == "group" && is.null(group)) {
    stop("`facet_by = 'group'` requires `group`.", call. = FALSE)
  }
  if (facet_by == "group" && !compare %in% c("groups", "none")) {
    stop(
      "`facet_by = 'group'` requires `compare = 'groups'` or `compare = 'none'`.",
      call. = FALSE
    )
  }

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

  facet_by_program <- facet_by == "program"
  facet_by_group <- facet_by == "group"
  overlay_by_program <- facet_by == "none" && compare %in% c("none", "programs")
  uses_program_linetype <- !facet_by_program && !overlay_by_program &&
    length(programs) > 1L
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
    if (facet_by == "none" && compare != "programs" && length(programs) > 1L) {
      step <- max(0.12, 0.1 * diff(range(data$enrichment_z)))
      tests$annotation_y <- tests$annotation_y +
        (match(as.character(tests$program), programs) - 1L) * step
      tests$label <- paste0(as.character(tests$program), ": ", tests$label)
    }
    if (facet_by_group) {
      # The inferential result compares groups and therefore belongs to neither
      # group alone. Draw it once in the first facet instead of repeating the
      # same label in every panel.
      tests$curve_group <- factor(levels_used[[1L]], levels = levels_used)
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
  if (!is.null(linetypes)) {
    if (!uses_program_linetype) {
      stop(
        "`linetypes` requires multiple programs sharing a panel or group facet.",
        call. = FALSE
      )
    }
    if (is.null(names(linetypes)) || !all(programs %in% names(linetypes))) {
      stop("`linetypes` must be named for every displayed program.", call. = FALSE)
    }
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
  summary_line_mapping <- if (uses_program_linetype) {
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
  if (facet_by_program) {
    plot <- plot + ggplot2::facet_wrap(~program, ncol = as.integer(facet_ncol))
  } else if (facet_by_group) {
    plot <- plot + ggplot2::facet_wrap(~curve_group, ncol = as.integer(facet_ncol))
  }
  if (!is.null(linetypes)) {
    plot <- plot + ggplot2::scale_linetype_manual(
      values = linetypes[programs], drop = FALSE
    )
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
      linetype = if (uses_program_linetype) {
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
