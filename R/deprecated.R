#' Deprecated SingleCellExperiment-prefixed API
#'
#' CellCharterR now uses neutral public names because its primary workflow is
#' SpatialExperiment-native. These wrappers remain temporarily so existing
#' scripts receive an actionable deprecation warning.
#'
#' @param ... Arguments forwarded to the replacement function.
#' @name deprecated-sce-api
#' @keywords internal
NULL

#' @rdname deprecated-sce-api
#' @export
cellcharter_sce_graph <- function(...) {
  .Deprecated("cellcharter_graph")
  cellcharter_graph(...)
}

#' @rdname deprecated-sce-api
#' @export
cellcharter_sce_aggregate <- function(...) {
  .Deprecated("cellcharter_aggregate")
  cellcharter_aggregate(...)
}

#' @rdname deprecated-sce-api
#' @export
cellcharter_sce_cluster <- function(...) {
  .Deprecated("cellcharter_cluster")
  cellcharter_cluster(...)
}

#' @rdname deprecated-sce-api
#' @export
cellcharter_sce_enrichment <- function(...) {
  .Deprecated("cellcharter_enrichment")
  cellcharter_enrichment(...)
}

#' @rdname deprecated-sce-api
#' @export
cellcharter_sce_nhood_enrichment <- function(...) {
  .Deprecated("cellcharter_nhood_enrichment")
  cellcharter_nhood_enrichment(...)
}

#' @rdname deprecated-sce-api
#' @export
cellcharter_sce_diff_nhood_enrichment <- function(...) {
  .Deprecated("cellcharter_diff_nhood_enrichment")
  cellcharter_diff_nhood_enrichment(...)
}

#' @rdname deprecated-sce-api
#' @export
cellcharter_sce_workflow <- function(...) {
  .Deprecated("cellcharter_workflow")
  cellcharter_workflow(...)
}
