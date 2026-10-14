.cc_diff_network_result <- function(x, result_name) {
  if (inherits(x, "SingleCellExperiment")) {
    record <- S4Vectors::metadata(x)$CellCharterR$results[[result_name]]
    if (is.null(record) || is.null(record$value)) {
      stop(
        "Experiment lacks stored differential neighborhood-enrichment result `",
        result_name, "`; run `cellcharter_diff_nhood_enrichment()` first.",
        call. = FALSE
      )
    }
    identity <- .cc_validate_sce_identity(
      x, context = paste0("Differential result `", result_name, "`")
    )
    if (!identical(record$n_cells, ncol(x)) ||
        !identical(record$cell_order_sha256, identity$cell_order_sha256)) {
      stop("Stored differential result `", result_name,
           "` is stale after cell reordering or subsetting; recompute it.",
           call. = FALSE)
    }
    return(record$value)
  }
  if (is.list(x) && !is.null(x$value) && is.list(x$value)) x <- x$value
  if (!is.list(x)) {
    stop("`x` must be an experiment or differential neighborhood-enrichment result.",
         call. = FALSE)
  }
  x
}

.cc_diff_network_validate_matrix <- function(effect, pvalue) {
  valid_matrix <- function(x) {
    is.matrix(x) && is.numeric(x) && length(dim(x)) == 2L &&
      nrow(x) > 0L && nrow(x) == ncol(x)
  }
  if (!valid_matrix(effect) || !valid_matrix(pvalue)) {
    stop("Differential `effect` and `pvalue` must be non-empty square numeric matrices.",
         call. = FALSE)
  }
  if (!identical(dim(effect), dim(pvalue)) ||
      !identical(dimnames(effect), dimnames(pvalue))) {
    stop("Differential `effect` and `pvalue` must have identical dimensions and dimnames.",
         call. = FALSE)
  }
  domains <- rownames(effect)
  if (is.null(domains) || is.null(colnames(effect)) ||
      !identical(domains, colnames(effect)) || anyNA(domains) ||
      any(!nzchar(domains)) || anyDuplicated(domains)) {
    stop("Differential matrices require identical, unique, non-empty row and column names.",
         call. = FALSE)
  }
  invisible(domains)
}

.cc_diff_nhood_network_data <- function(
    x, result_name = "diff_nhood_enrichment", alpha = 0.05,
    p_adjust_method = "none", effect_threshold = 0,
    significant_only = FALSE, top_n = NULL) {
  result <- .cc_diff_network_result(x, result_name)
  if (is.null(result$effect) || is.null(result$pvalue) || is.null(result$params)) {
    stop("Malformed differential result: `effect`, `pvalue`, and `params` are required.",
         call. = FALSE)
  }
  domains <- .cc_diff_network_validate_matrix(result$effect, result$pvalue)
  groups <- as.character(result$params$condition_level_order)
  if (length(groups) != 2L || anyNA(groups) || any(!nzchar(groups)) ||
      anyDuplicated(groups)) {
    stop("Differential result lacks two valid inferential condition groups.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha > 1) {
    stop("`alpha` must be one finite number in (0, 1].", call. = FALSE)
  }
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) || !p_adjust_method %in% stats::p.adjust.methods) {
    stop("`p_adjust_method` must be one value from `p.adjust.methods`.",
         call. = FALSE)
  }
  if (!is.numeric(effect_threshold) || length(effect_threshold) != 1L ||
      is.na(effect_threshold) || !is.finite(effect_threshold) ||
      effect_threshold < 0) {
    stop("`effect_threshold` must be one non-negative finite number.",
         call. = FALSE)
  }
  if (!is.logical(significant_only) || length(significant_only) != 1L ||
      is.na(significant_only)) {
    stop("`significant_only` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.null(top_n) &&
      (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
       !is.finite(top_n) || top_n < 1 || top_n != floor(top_n))) {
    stop("`top_n` must be `NULL` or one positive whole number.", call. = FALSE)
  }

  effect <- result$effect
  pvalue <- result$pvalue
  indices <- which(row(effect) != col(effect) & is.finite(effect) & effect != 0,
                   arr.ind = TRUE)
  if (nrow(indices) &&
      any(!is.finite(pvalue[indices]) | pvalue[indices] < 0 | pvalue[indices] > 1)) {
    stop("Eligible directed effects require finite p-values between zero and one.",
         call. = FALSE)
  }
  adjusted <- if (nrow(indices)) {
    stats::p.adjust(pvalue[indices], method = p_adjust_method)
  } else numeric()
  edges <- data.frame(
    source_domain = domains[indices[, 1L]],
    neighbor_domain = domains[indices[, 2L]],
    signed_effect = effect[indices],
    display_effect = abs(effect[indices]),
    favoured_condition = ifelse(effect[indices] > 0, groups[[1L]], groups[[2L]]),
    pvalue = pvalue[indices],
    adjusted_pvalue = adjusted,
    significant = adjusted <= alpha,
    source_order = indices[, 1L], neighbor_order = indices[, 2L],
    stringsAsFactors = FALSE
  )
  edges <- edges[edges$display_effect >= effect_threshold, , drop = FALSE]
  if (isTRUE(significant_only)) edges <- edges[edges$significant, , drop = FALSE]
  if (nrow(edges)) {
    edges$favoured_condition <- factor(edges$favoured_condition, levels = groups)
    edges <- edges[order(
      edges$favoured_condition, -edges$display_effect,
      edges$source_order, edges$neighbor_order
    ), , drop = FALSE]
    if (!is.null(top_n)) {
      rank_in_panel <- stats::ave(
        seq_len(nrow(edges)), edges$favoured_condition, FUN = seq_along
      )
      edges <- edges[rank_in_panel <= as.integer(top_n), , drop = FALSE]
    }
    reciprocal_key <- paste(edges$neighbor_domain, edges$source_domain, sep = "\r")
    edge_key <- paste(edges$source_domain, edges$neighbor_domain, sep = "\r")
    edges$curvature <- ifelse(reciprocal_key %in% edge_key, 0.18, 0)
  } else {
    edges$favoured_condition <- factor(character(), levels = groups)
    edges$curvature <- numeric()
  }
  edges$source_domain <- factor(edges$source_domain, levels = domains)
  edges$neighbor_domain <- factor(edges$neighbor_domain, levels = domains)
  rownames(edges) <- NULL
  attr(edges, "condition_groups") <- groups
  attr(edges, "p_adjust_method") <- p_adjust_method
  attr(edges, "alpha") <- alpha
  attr(edges, "empirical_tail") <- result$params$empirical_tail %||%
    "directional one-sided empirical"
  edges
}

.cc_diff_network_colours <- function(values, keys, argument, default) {
  if (is.null(values)) values <- default
  if (!is.character(values) || anyNA(values) || !length(values)) {
    stop("`", argument, "` must contain valid R colour strings.", call. = FALSE)
  }
  if (length(values) == 1L) values <- rep(values, length(keys))
  if (!is.null(names(values))) {
    if (!all(keys %in% names(values))) {
      stop("Named `", argument, "` must provide every required name.", call. = FALSE)
    }
    values <- values[keys]
  } else if (length(values) != length(keys)) {
    stop("Unnamed `", argument, "` must have length one or match the required names.",
         call. = FALSE)
  }
  tryCatch(grDevices::col2rgb(values), error = function(e) {
    stop("`", argument, "` contains an invalid R colour.", call. = FALSE)
  })
  stats::setNames(unname(values), keys)
}

#' Plot differential neighborhood enrichment as directed networks
#'
#' Uses a stored result from [cellcharter_diff_nhood_enrichment()] to draw two
#' panels with the same node positions. Positive effects appear in the first
#' inferential condition and negative effects in the second. The permutation
#' analysis is never recomputed. Matrix rows are source/central domains and
#' columns are surrounding/neighbor domains.
#'
#' @param x A `SingleCellExperiment`/`SpatialExperiment` containing a stored
#'   differential result, its stored result record, or the raw result list.
#' @param result_name Stored result name, used only for an experiment input.
#' @param alpha Significance threshold in `(0, 1]`.
#' @param p_adjust_method Method passed to [stats::p.adjust()]. The default
#'   `"none"` matches the reference CellCharter display.
#' @param effect_threshold Minimum absolute effect retained for display.
#' @param significant_only Display only statistically supported edges.
#' @param top_n Optional maximum number of strongest edges per condition panel.
#' @param edge_colours One colour for every arrow, or two colours optionally
#'   named by condition.
#' @param node_colours One node colour, or one colour per domain optionally
#'   named by domain.
#' @param show_edge_values Label arrows with absolute displayed effects.
#' @return A customizable faceted `ggplot` object. Its transformed edge and
#'   node tables are available in the `cellcharter_network_data` and
#'   `cellcharter_node_coordinates` attributes.
#' @examples
#' spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' cellcharter_plot_graph_diff(spe, top_n = 10)
#' @export
cellcharter_plot_graph_diff <- function(
    x, result_name = "diff_nhood_enrichment", alpha = 0.05,
    p_adjust_method = "none", effect_threshold = 0,
    significant_only = FALSE, top_n = NULL, edge_colours = "black",
    node_colours = NULL, show_edge_values = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install `ggplot2` to plot differential neighborhood networks.",
         call. = FALSE)
  }
  if (!is.logical(show_edge_values) || length(show_edge_values) != 1L ||
      is.na(show_edge_values)) {
    stop("`show_edge_values` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  edges <- .cc_diff_nhood_network_data(
    x, result_name, alpha, p_adjust_method, effect_threshold,
    significant_only, top_n
  )
  groups <- attr(edges, "condition_groups")
  domains <- levels(edges$source_domain)
  edge_colours <- .cc_diff_network_colours(
    edge_colours, groups, "edge_colours", "black"
  )
  node_colours <- .cc_diff_network_colours(
    node_colours, domains, "node_colours", "grey85"
  )
  n_domains <- length(domains)
  angles <- pi / 2 - 2 * pi * (seq_len(n_domains) - 1L) / n_domains
  nodes <- data.frame(
    domain = factor(domains, levels = domains),
    x = cos(angles), y = sin(angles),
    label_x = 1.22 * cos(angles), label_y = 1.22 * sin(angles)
  )
  nodes <- merge(
    data.frame(favoured_condition = factor(groups, levels = groups)),
    nodes, all = TRUE
  )
  position <- stats::setNames(seq_len(n_domains), domains)
  if (nrow(edges)) {
    source <- position[as.character(edges$source_domain)]
    neighbor <- position[as.character(edges$neighbor_domain)]
    dx <- cos(angles[neighbor]) - cos(angles[source])
    dy <- sin(angles[neighbor]) - sin(angles[source])
    distance <- sqrt(dx^2 + dy^2)
    inset <- 0.11
    edges$x <- cos(angles[source]) + inset * dx / distance
    edges$y <- sin(angles[source]) + inset * dy / distance
    edges$xend <- cos(angles[neighbor]) - inset * dx / distance
    edges$yend <- sin(angles[neighbor]) - inset * dy / distance
    edges$label_x <- (edges$x + edges$xend) / 2
    edges$label_y <- (edges$y + edges$yend) / 2
    edges$significance <- factor(
      ifelse(edges$significant, "Supported", "Not statistically supported"),
      levels = c("Supported", "Not statistically supported")
    )
  } else {
    edges$x <- edges$y <- edges$xend <- edges$yend <- numeric()
    edges$label_x <- edges$label_y <- numeric()
    edges$significance <- factor(
      character(), levels = c("Supported", "Not statistically supported")
    )
  }
  if (nrow(edges) > 100L) {
    warning("The differential network contains more than 100 arrows; consider `top_n`, `effect_threshold`, or `significant_only = TRUE` for readability.",
            call. = FALSE)
  }

  basis <- if (identical(p_adjust_method, "none")) "raw" else
    paste0(p_adjust_method, "-adjusted")
  significance_title <- paste0(
    basis, " directional empirical p ", "\u2264 ", format(alpha)
  )
  plot <- ggplot2::ggplot()
  for (curve in sort(unique(edges$curvature))) {
    layer_data <- edges[edges$curvature == curve, , drop = FALSE]
    plot <- plot + ggplot2::geom_curve(
      data = layer_data,
      ggplot2::aes(
        x = x, y = y, xend = xend, yend = yend,
        colour = favoured_condition, linetype = significance
      ),
      curvature = curve, linewidth = 0.7, alpha = 0.85,
      arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed")
    )
  }
  plot <- plot +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, fill = domain),
      shape = 21, size = 7, colour = "grey15", stroke = 0.8
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = label_x, y = label_y, label = domain),
      size = 4.2, lineheight = 0.9
    )
  if (isTRUE(show_edge_values) && nrow(edges)) {
    plot <- plot + ggplot2::geom_text(
      data = edges,
      ggplot2::aes(x = label_x, y = label_y, label = sprintf("%.3g", display_effect)),
      size = 3.8, colour = "grey10"
    )
  }
  panel_counts <- table(factor(edges$favoured_condition, levels = groups))
  empty_groups <- groups[panel_counts == 0L]
  empty <- data.frame(
    favoured_condition = factor(empty_groups, levels = groups),
    x = rep(0, length(empty_groups)), y = rep(0, length(empty_groups)),
    label = rep("No edges pass the selected filters", length(empty_groups))
  )
  if (nrow(empty)) {
    plot <- plot + ggplot2::geom_text(
      data = empty, ggplot2::aes(x = x, y = y, label = label),
      size = 4.2, colour = "grey35"
    )
  }
  plot <- plot +
    ggplot2::facet_wrap(~favoured_condition, nrow = 1, drop = FALSE) +
    ggplot2::scale_colour_manual(values = edge_colours, name = "Favoured condition") +
    ggplot2::scale_fill_manual(values = node_colours, guide = "none") +
    ggplot2::scale_linetype_manual(
      values = c("Supported" = "solid", "Not statistically supported" = "22"),
      name = significance_title, drop = FALSE
    ) +
    ggplot2::coord_equal(xlim = c(-1.55, 1.55), ylim = c(-1.45, 1.45), clip = "off") +
    ggplot2::labs(
      title = "Differential neighborhood enrichment",
      subtitle = paste0(
        groups[[1L]], " - ", groups[[2L]],
        "; arrows show source/central ", "\u2192", " surrounding/neighbor"
      ),
      caption = paste0(
        "Solid arrows satisfy ", significance_title,
        "; dashed arrows are not statistically supported."
      )
    ) +
    ggplot2::theme_void(base_size = 16) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 16),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(18, 35, 18, 35),
      legend.position = "bottom"
    )
  attr(plot, "cellcharter_network_data") <- edges
  attr(plot, "cellcharter_node_coordinates") <- nodes
  plot
}

utils::globalVariables(c(
  "source_domain", "neighbor_domain", "signed_effect", "display_effect",
  "favoured_condition", "adjusted_pvalue", "significant", "curvature",
  "x", "y", "xend", "yend", "label_x", "label_y", "domain",
  "significance", "label"
))
