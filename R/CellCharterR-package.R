#' CellCharterR: disk-backed spatial-domain analysis
#'
#' CellCharterR provides a `SpatialExperiment`-native implementation of the
#' post-representation CellCharter workflow. Large assays may remain
#' disk-backed while spatial coordinates, graph edges, reduced representations,
#' annotations, and analysis results stay in their standard SPE compartments.
#' Start with [cellcharter_workflow()] for a compact
#' example or use the staged graph, aggregation, clustering, and enrichment
#' adapters when durable checkpoints are required.
#'
#' CellCharterR is an independent R implementation of methods described by the
#' CellCharter project. It is not affiliated with, endorsed by, or maintained
#' by the original CellCharter authors or the Computational Systems Oncology
#' group. Questions and bug reports concerning CellCharterR should be directed
#' to <https://github.com/cperalta-viteri/CellCharterR/issues>, not to the
#' upstream CellCharter project.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "annotation_y", "block", "curve_group", "enrichment_z", "plot_colour",
  "plot_group", "program", "se", "zone"
))
