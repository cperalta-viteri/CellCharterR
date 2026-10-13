#' Plot CellCharter cell-type enrichment
#'
#' Draws a dot plot in which colour represents signed log2 observed-to-expected
#' enrichment and point area represents its magnitude. The input can be a
#' clustered `SpatialExperiment`, a result returned by the enrichment kernel,
#' or an enrichment matrix.
#'
#' @param x A `SingleCellExperiment`/`SpatialExperiment` containing a stored
#'   enrichment result, a list with an `enrichment` matrix, or the matrix itself.
#' @param result_name Result name under
#'   `metadata(x)$CellCharterR$results`; used only when `x` is an experiment.
#' @param title Optional plot title.
#' @param point_range Two positive numbers controlling the point-size range.
#' @return A `ggplot` object.
#' @examples
#' spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' cellcharter_plot_enrichment(spe)
#' @export
cellcharter_plot_enrichment <- function(
    x, result_name = "enrichment", title = NULL, point_range = c(1, 8)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install `ggplot2` to plot CellCharter enrichment.", call. = FALSE)
  }
  if (inherits(x, "SingleCellExperiment")) {
    record <- S4Vectors::metadata(x)$CellCharterR$results[[result_name]]
    if (is.null(record$value$enrichment)) {
      stop("SCE lacks stored enrichment result `", result_name,
           "`; run `cellcharter_enrichment()` first.", call. = FALSE)
    }
    values <- record$value$enrichment
  } else if (is.list(x) && !is.null(x$enrichment)) {
    values <- x$enrichment
  } else {
    values <- x
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

  plot_values <- values
  finite <- values[is.finite(values)]
  if (!length(finite)) stop("The enrichment matrix has no finite values.", call. = FALSE)
  plot_values[is.infinite(values) & values < 0] <- min(finite)
  plot_values[is.infinite(values) & values > 0] <- max(finite)

  data <- as.data.frame(as.table(plot_values), stringsAsFactors = FALSE)
  names(data) <- c("group", "domain", "enrichment")
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = domain, y = group,
      colour = enrichment, size = abs(enrichment)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "log2 FC"
    ) +
    ggplot2::scale_size_continuous(name = "|log2 FC|", range = point_range) +
    ggplot2::labs(x = "CellCharter domain", y = "Group", title = title) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

utils::globalVariables(c("domain", "group", "enrichment"))
