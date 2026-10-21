# SpatialExperiment-explicit aliases. The concise cellcharter_* entry points
# remain canonical; each alias is documented on its canonical help page so F1
# exposes the complete parameter, return-value, and storage contract.

#' @rdname cellcharter_graph
#' @param ... Additional arguments accepted by [cellcharter_graph()].
#' @export
cellcharter_spe_graph <- function(spe, ...) {
  cellcharter_graph(spe, ...)
}

#' @rdname cellcharter_aggregate
#' @param ... Additional arguments accepted by [cellcharter_aggregate()].
#' @export
cellcharter_spe_aggregate <- function(spe, ...) {
  cellcharter_aggregate(spe, ...)
}

#' @rdname cellcharter_cluster
#' @param ... Additional arguments accepted by [cellcharter_cluster()].
#' @export
cellcharter_spe_cluster <- function(spe, ...) {
  cellcharter_cluster(spe, ...)
}

#' @rdname cellcharter_enrichment
#' @param ... Additional arguments accepted by [cellcharter_enrichment()].
#' @export
cellcharter_spe_enrichment <- function(spe, ...) {
  cellcharter_enrichment(spe, ...)
}

#' @rdname cellcharter_nhood_enrichment
#' @param ... Additional arguments accepted by
#'   [cellcharter_nhood_enrichment()].
#' @export
cellcharter_spe_nhood_enrichment <- function(spe, ...) {
  cellcharter_nhood_enrichment(spe, ...)
}

#' @rdname cellcharter_diff_nhood_enrichment
#' @param ... Additional arguments accepted by
#'   [cellcharter_diff_nhood_enrichment()].
#' @export
cellcharter_spe_diff_nhood_enrichment <- function(spe, ...) {
  cellcharter_diff_nhood_enrichment(spe, ...)
}

#' @rdname cellcharter_workflow
#' @param ... Additional arguments accepted by [cellcharter_workflow()].
#' @export
cellcharter_spe_workflow <- function(spe, ...) {
  cellcharter_workflow(spe, ...)
}
