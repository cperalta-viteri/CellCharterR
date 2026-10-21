#' Plot CellCharter cell-type enrichment
#'
#' Draws a CellCharter-style dot plot that focuses on positive log2
#' observed-to-expected enrichment. Non-positive values remain in the returned
#' statistical result but are omitted from the plot. The input can be a
#' clustered `SpatialExperiment`, a result returned by the enrichment kernel,
#' or an enrichment matrix. For experiment input, stored scientific-input and
#' parameter fingerprints are checked before plotting.
#'
#' @param x A `SpatialExperiment` containing a stored enrichment result, a
#'   result list with a numeric log2-enrichment matrix, or a numeric matrix that
#'   the caller asserts already contains log2 fold changes. Stored and direct
#'   result lists explicitly marked `log2_fold_change = FALSE` are rejected so
#'   raw ratios cannot be mislabeled as positive log2 fold changes. Plotting is
#'   read-only and never changes `x`.
#' @param result_name One non-empty key under
#'   `metadata(x)$CellCharterR$results`. It is read only when `x` is a
#'   `SpatialExperiment`; stored input and parameter fingerprints are validated
#'   before the result is plotted.
#' @param title `NULL` for no title, or a character value passed to the plot's
#'   title label.
#' @param point_range Two ordered positive finite numbers giving the minimum
#'   and maximum rendered point sizes.
#' @param size_threshold One positive finite log2 fold-change at which point
#'   size is capped. The CellCharter CODEX tutorial uses `0.58`.
#' @param colour_threshold One positive finite log2 fold-change at which the
#'   dark-red colour intensity is capped. It also separates weakly positive
#'   from enriched values for hierarchical ordering. The CellCharter CODEX
#'   tutorial uses `0.58`.
#' @param cluster_domains A single logical. When `TRUE`, domains are ordered by
#'   complete-linkage hierarchical clustering of ternary enrichment states:
#'   negative, weakly positive, and at least `colour_threshold`.
#' @param cluster_programs A single logical applying the same complete-linkage
#'   ordering to programs rather than retaining the input row order.
#' @return A read-only, customizable `ggplot` object. Non-positive and
#'   non-finite enrichment entries are omitted from the displayed points but
#'   remain unchanged in the stored result.
#' @examples
#' spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' cellcharter_plot_enrichment(spe)
#' @export
cellcharter_plot_enrichment <- function(
    x, result_name = "enrichment", title = NULL, point_range = c(1, 8),
    size_threshold = 0.58, colour_threshold = 0.58,
    cluster_domains = TRUE, cluster_programs = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install `ggplot2` to plot CellCharter enrichment.", call. = FALSE)
  }
  if (inherits(x, "SingleCellExperiment") &&
      !inherits(x, "SpatialExperiment")) {
    stop("Experiment input `x` must be a SpatialExperiment.", call. = FALSE)
  }
  if (inherits(x, "SpatialExperiment")) {
    record <- .cc_spe_require_result(
      x, result_name, context = "Stored enrichment result"
    )
    if (is.null(record$value$enrichment)) {
      stop("SpatialExperiment lacks stored enrichment result `", result_name,
           "`; run `cellcharter_enrichment()` first.", call. = FALSE)
    }
    values <- record$value$enrichment
    log2_fold_change <- record$parameters$log2_fold_change %||%
      record$value$params$log2_fold_change
  } else if (is.list(x) && !is.null(x$enrichment)) {
    values <- x$enrichment
    log2_fold_change <- x$params$log2_fold_change
  } else {
    values <- x
    log2_fold_change <- TRUE
  }
  if (!is.null(log2_fold_change) && !isTRUE(log2_fold_change)) {
    stop(
      "`cellcharter_plot_enrichment()` requires log2 fold-change enrichment; ",
      "recompute with `log2_fold_change = TRUE`.", call. = FALSE
    )
  }
  if (!is.matrix(values) || !is.numeric(values) ||
      any(dim(values) == 0L) || is.null(rownames(values)) ||
      is.null(colnames(values))) {
    stop("`x` must resolve to a non-empty numeric enrichment matrix with row and column names.",
         call. = FALSE)
  }
  if (!is.numeric(point_range) || length(point_range) != 2L ||
      any(!is.finite(point_range)) || any(point_range <= 0) ||
      point_range[[1L]] > point_range[[2L]]) {
    stop("`point_range` must contain two ordered positive finite numbers.", call. = FALSE)
  }
  if (!is.numeric(size_threshold) || length(size_threshold) != 1L ||
      !is.finite(size_threshold) || size_threshold <= 0) {
    stop("`size_threshold` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(colour_threshold) || length(colour_threshold) != 1L ||
      !is.finite(colour_threshold) || colour_threshold <= 0) {
    stop("`colour_threshold` must be one positive finite number.", call. = FALSE)
  }
  if (!is.logical(cluster_domains) || length(cluster_domains) != 1L ||
      is.na(cluster_domains)) {
    stop("`cluster_domains` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.logical(cluster_programs) || length(cluster_programs) != 1L ||
      is.na(cluster_programs)) {
    stop("`cluster_programs` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  finite <- values[is.finite(values)]
  if (!length(finite)) {
    stop("The enrichment matrix has no finite values.", call. = FALSE)
  }

  # Match CellCharter's ordering states: depleted, weakly positive, and
  # enriched rather than allowing unbounded fold changes to dominate distance.
  colour_states <- ifelse(
    values < 0, -1,
    ifelse(values < colour_threshold, 0, 1)
  )
  nonfinite <- !is.finite(colour_states)
  colour_states[nonfinite] <- ifelse(
    is.infinite(values[nonfinite]) & values[nonfinite] > 0, 1, -1
  )

  domain_order <- colnames(values)
  if (cluster_domains && ncol(values) > 1L) {
    domain_order <- domain_order[stats::hclust(
      stats::dist(t(colour_states)), method = "complete"
    )$order]
  }
  program_order <- rownames(values)
  if (cluster_programs && nrow(values) > 1L) {
    program_order <- program_order[stats::hclust(
      stats::dist(colour_states), method = "complete"
    )$order]
  }

  data <- as.data.frame(as.table(values), stringsAsFactors = FALSE)
  names(data) <- c("group", "domain", "enrichment")
  data$group <- factor(data$group, levels = program_order)
  data$domain <- factor(data$domain, levels = rev(domain_order))
  # CellCharter's signed size normalization is clipped at zero by Scanpy.
  # Keeping only positive finite values makes that behavior explicit and also
  # prevents absent group-domain combinations (-Inf) from appearing as dots.
  data <- data[is.finite(data$enrichment) & data$enrichment > 0, , drop = FALSE]
  data$point_size <- pmin(data$enrichment, size_threshold)
  data$colour_value <- pmin(data$enrichment, colour_threshold)

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = group, y = domain,
      colour = colour_value, size = point_size
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(
      low = "#FDE0DD", high = "#67000D",
      limits = c(0, colour_threshold), name = "Positive log2 FC"
    ) +
    ggplot2::scale_size_continuous(
      limits = c(0, size_threshold), name = "Positive log2 FC",
      range = point_range
    ) +
    ggplot2::labs(x = "Cell type", y = "CellCharter domain", title = title) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

utils::globalVariables(c(
  "domain", "group", "enrichment", "point_size", "colour_value"
))
