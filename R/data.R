#' Construct the portable mouse-spleen teaching SpatialExperiment
#'
#' Reconstructs the 1,000-segmented-cell CODEX teaching object from
#' version-neutral matrices
#' and data frames. This avoids serialized Bioconductor class incompatibilities
#' between R 4.4/Bioconductor 3.20 and newer Bioconductor releases.
#' The subset retains the protein assay, fixed trVAE representation, spatial
#' coordinates, annotations, K=4 labels, and K=11 sensitivity labels. It is for
#' examples and teaching, not for reproducing the full biological analysis.
#'
#' @param with_results If `TRUE`, also restore the validated teaching workflow's
#'   graph, AutoK summary, labels, and downstream statistical results. The
#'   default returns the clean input object used to demonstrate the workflow.
#' @return A 29-feature by 1,000-segmented-cell `SpatialExperiment`.
#' @source CODEX mouse-spleen data, dataset DOI
#'   \doi{10.17632/zjnpwh8m5b.1}; original CODEX publication
#'   \doi{10.1016/j.cell.2018.07.010}.
#' @export
#' @examples
#' spe <- cellcharter_mouse_spleen_tiny()
cellcharter_mouse_spleen_tiny <- function(with_results = FALSE) {
  if (!is.logical(with_results) || length(with_results) != 1L || is.na(with_results)) {
    stop("`with_results` must be TRUE or FALSE.", call. = FALSE)
  }
  filename <- if (with_results) {
    "codex_mouse_spleen_tiny_workflow_components.rds"
  } else {
    "codex_mouse_spleen_tiny_components.rds"
  }
  path <- system.file(
    "extdata", filename,
    package = "CellCharterR", mustWork = TRUE
  )
  components <- readRDS(path)
  row_data <- S4Vectors::DataFrame(components$row_data)
  col_data <- S4Vectors::DataFrame(components$col_data)
  # data.frame() sanitizes this original AnnData annotation name.
  colnames(col_data)[colnames(col_data) == "i.niche"] <- "i-niche"
  rownames(row_data) <- components$row_names
  rownames(col_data) <- components$col_names
  assays <- S4Vectors::SimpleList(components$assays)
  for (name in names(assays)) {
    rownames(assays[[name]]) <- components$row_names
    colnames(assays[[name]]) <- components$col_names
  }
  spatial <- components$spatial
  rownames(spatial) <- components$col_names
  spe <- SpatialExperiment::SpatialExperiment(
    assays = assays,
    rowData = row_data,
    colData = col_data,
    spatialCoords = spatial,
    sample_id = as.character(col_data$sample_id)
  )
  for (name in names(components$reduced_dims)) {
    value <- components$reduced_dims[[name]]
    rownames(value) <- components$col_names
    SingleCellExperiment::reducedDim(spe, name) <- value
  }
  for (name in names(components$col_pairs)) {
    pair <- components$col_pairs[[name]]
    hits <- S4Vectors::SelfHits(pair$from, pair$to, nnode = ncol(spe))
    extra <- setdiff(colnames(pair), c("from", "to"))
    if (length(extra)) S4Vectors::mcols(hits) <- S4Vectors::DataFrame(pair[extra])
    SingleCellExperiment::colPair(spe, name) <- hits
  }
  S4Vectors::metadata(spe) <- components$metadata
  spe
}
