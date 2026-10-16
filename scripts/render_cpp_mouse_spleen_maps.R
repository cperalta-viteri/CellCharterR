#!/usr/bin/env Rscript

# Render the two README spatial maps from the completed compiled-C++ labels.
# This script only reads saved report inputs; it never refits a model.

input_root <- "results/CellCharterR_mouse_spleen_cpp_v1/report_artifacts"
output_root <- "docs/images"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

map_to_spe <- function(data, k) {
  samples <- c("BALBc-1", "MRL-4")
  keep <- as.character(data$sample) %in% samples
  data <- data[keep, , drop = FALSE]
  SpatialExperiment::SpatialExperiment(
    assays = list(display = matrix(0, nrow = 1L, ncol = nrow(data))),
    colData = S4Vectors::DataFrame(
      sample_id = factor(
        as.character(data$sample), levels = samples,
        labels = c("BALBc-1", "MRL-4 (early)")
      ),
      domain = factor(data$domain, levels = seq.int(0L, k - 1L))
    ),
    spatialCoords = as.matrix(data[, c("x", "y")]),
    sample_id = factor(
      as.character(data$sample), levels = samples,
      labels = c("BALBc-1", "MRL-4 (early)")
    )
  )
}

render_map <- function(relative_input, k, title, output_name) {
  input <- file.path(input_root, relative_input, "map_input.rds")
  if (!file.exists(input)) stop("Missing saved map input: ", input)
  map <- readRDS(input)
  spe <- map_to_spe(map, k)
  plot <- ggspavis::plotCoords(
    spe,
    annotate = "domain",
    sample_id = "sample_id",
    in_tissue = NULL,
    point_size = 0.1,
    pal = grDevices::hcl.colors(k, "Dynamic")
  ) +
    ggplot2::labs(colour = "Domain", title = title) +
    ggplot2::theme_minimal(base_size = 16)
  ggplot2::ggsave(
    file.path(output_root, output_name), plot,
    width = 12, height = 5.5, units = "in", dpi = 150, bg = "white"
  )
}

render_map(
  "primary_AutoK_K10", 10L,
  "Compiled C++ AutoK primary spatial domains: K = 10",
  "mouse-spleen-cpp-k10-balbc1-mrl4.png"
)
render_map(
  "tutorial_standalone_K11", 11L,
  "Tutorial-matched sensitivity spatial domains: K = 11",
  "mouse-spleen-k11-balbc1-mrl4.png"
)

message("Rendered compiled-C++ mouse-spleen maps in ", output_root)
