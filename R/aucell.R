# AUCell-formula scoring for disk-backed spatial experiments.

.cc_aucell_scalar_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop("`", name, "` must be one integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.cc_aucell_validate_name <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be one non-empty string.", call. = FALSE)
  }
  x
}

.cc_aucell_prepare_gene_sets <- function(gene_sets, features) {
  if (!is.list(gene_sets) || !length(gene_sets)) {
    stop("`gene_sets` must be a non-empty named list of character vectors.",
         call. = FALSE)
  }
  program_names <- names(gene_sets)
  if (is.null(program_names) || anyNA(program_names) ||
      any(!nzchar(program_names)) || anyDuplicated(program_names)) {
    stop("`gene_sets` must have unique, non-empty program names.", call. = FALSE)
  }
  valid <- vapply(gene_sets, function(x) {
    is.character(x) && !anyNA(x) && all(nzchar(x))
  }, logical(1L))
  if (!all(valid)) {
    stop("Every element of `gene_sets` must be a character vector without missing values.",
         call. = FALSE)
  }

  gene_sets <- lapply(gene_sets, unique)
  feature_index <- stats::setNames(seq_along(features), features)
  indices <- lapply(gene_sets, function(x) {
    unname(feature_index[intersect(x, features)])
  })
  available <- lengths(indices)
  if (any(available == 0L)) {
    stop("Gene set(s) have no features in the assay: ",
         paste(program_names[available == 0L], collapse = ", "), ".",
         call. = FALSE)
  }
  missing_fraction <- 1 - available / pmax(lengths(gene_sets), 1L)
  mostly_missing <- missing_fraction >= 0.8
  if (any(mostly_missing)) {
    stop("At least 80% of the genes are absent for gene set(s): ",
         paste(program_names[mostly_missing], collapse = ", "),
         ". Revise these signatures before scoring.", call. = FALSE)
  }
  names(indices) <- program_names
  indices
}

.cc_aucell_chunk_width <- function(n_features, n_programs, threads, memory_mib) {
  budget <- memory_mib * 1024^2
  # Each worker keeps one feature-sized value/key/index buffer. The R side may
  # transiently hold the extracted matrix, its dense coercion, and a native
  # numeric coercion, so reserve three doubles per assay value plus 64 MiB.
  # Use 32 bytes per ranking item as a portable upper budget; the current
  # x86_64 C++ layout is smaller, but padding is implementation-defined.
  scratch <- threads * n_features * 32
  safety <- 64 * 1024^2
  available <- budget - scratch - safety
  if (available <= 3 * n_features * 8) {
    stop("`memory_mib` is too small for one assay column and native worker scratch; ",
         "reduce `threads` or increase `memory_mib`.", call. = FALSE)
  }
  max(1L, min(.Machine$integer.max,
              floor(available / (3 * n_features * 8 + n_programs * 8))))
}

.cc_aucell_schedule <- function(n_spots, chunk_width) {
  seq.int(1L, n_spots, by = chunk_width)
}

#' Calculate AUCell scores for a disk-backed SpatialExperiment
#'
#' Ranks assay features independently within each spot and calculates the
#' normalized AUCell area for every supplied gene program. Assay data are read
#' in bounded chunks, while the final spot-by-program score matrix is held in
#' memory and stored in a reduced dimension.
#'
#' Equal expression values use deterministic counter-based random ordering by
#' default. This has the random-tie semantics of AUCell and is invariant to
#' chunking and thread scheduling, but is not a bit-identical reproduction of
#' AUCell's R random-number stream. Use `tie_method = "ordinal"` to resolve
#' ties by assay feature order, as in the decoupler ordinal implementation.
#' Feature identifiers and program names must be unique and non-empty. Marker
#' genes are deduplicated before scoring. A program with no available markers,
#' or with at least 80 percent of its markers absent, is rejected rather than
#' assigned a potentially misleading score. Finite negative assay values are
#' valid because the calculation depends only on ranks.
#'
#' @param spe A [SpatialExperiment::SpatialExperiment] with unique, non-empty
#'   row names and a non-empty numeric assay. Rows are features and columns are
#'   spatial spots; dense, sparse, and disk-backed assays are accepted.
#' @param gene_sets A non-empty named list of character vectors containing
#'   assay feature identifiers. Program names must be unique and non-empty.
#'   Duplicate markers within a program are removed, absent markers are
#'   ignored, programs with no matches are rejected, and programs with at
#'   least 80 percent of markers absent are rejected.
#' @param assay_name One non-empty string naming the feature-by-spot assay to
#'   rank independently within each spot.
#' @param auc_max_rank One finite number whose rounded value must be between 2
#'   and the assay feature count. The AUCell sum includes ranks strictly below
#'   this cutoff; the default is five percent of all features.
#' @param norm_auc A single logical that currently must be `TRUE`. Scores are
#'   divided by the maximum attainable recovery area for the available marker
#'   count and chosen cutoff.
#' @param tie_method Either `"random"` for deterministic counter-based random
#'   ordering of equal values or `"ordinal"` for assay-row-order ties. Random
#'   results are invariant to chunking and thread scheduling for a fixed seed,
#'   but do not reproduce AUCell's R random-number stream bit for bit.
#' @param output One non-empty name for the reduced dimension receiving the in-memory
#'   spot-by-program matrix. An existing reduced dimension of this name is
#'   replaced only after successful completion.
#' @param memory_mib Working-memory budget in MiB for assay chunks, coercion
#'   workspace, and native worker scratch. It excludes the final score matrix,
#'   the already loaded `SpatialExperiment`, package/runtime overhead, and
#'   transient allocator overhead. The final matrix alone requires
#'   `8 * ncol(spe) * length(gene_sets)` bytes.
#' @param threads A whole number from 1 through 32 giving the number of native
#'   CPU worker threads used within each assay chunk. It does not change the
#'   deterministic result.
#' @param seed Non-negative integer seed for deterministic random tie ordering.
#'   It has no effect when `tie_method = "ordinal"`.
#' @param progress A single logical indicating whether to display progress over
#'   assay chunks. It does not affect the result.
#'
#' @return The updated `SpatialExperiment`, with a dense, in-memory,
#'   spot-by-program double matrix in `reducedDim(spe, output)`. Row names are
#'   spot identifiers and column names are gene-set names. No `colData`,
#'   `metadata`, checkpoint, or external file is created.
#' @examples
#' expression <- matrix(
#'   c(8, 7, 2, 1, 0, 0, 0, 1, 2, 7, 8, 6),
#'   nrow = 6,
#'   dimnames = list(paste0("gene", 1:6), c("spot1", "spot2"))
#' )
#' spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(expression = expression),
#'   spatialCoords = cbind(x = 1:2, y = 0)
#' )
#' scored <- cellcharter_aucell(
#'   spe,
#'   gene_sets = list(program_a = c("gene1", "gene2"),
#'                    program_b = c("gene5", "gene6")),
#'   auc_max_rank = 4,
#'   threads = 1,
#'   progress = FALSE
#' )
#' SingleCellExperiment::reducedDim(scored, "AUCell")
#' @export
cellcharter_aucell <- function(
    spe,
    gene_sets,
    assay_name = "expression",
    auc_max_rank = ceiling(0.05 * nrow(spe)),
    norm_auc = TRUE,
    tie_method = c("random", "ordinal"),
    output = "AUCell",
    memory_mib = 1024,
    threads = 8L,
    seed = 5410L,
    progress = interactive()) {
  if (!methods::is(spe, "SpatialExperiment")) {
    stop("`spe` must be a SpatialExperiment.", call. = FALSE)
  }
  assay_name <- .cc_aucell_validate_name(assay_name, "assay_name")
  output <- .cc_aucell_validate_name(output, "output")
  if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
    stop("`spe` has no assay named `", assay_name, "`.", call. = FALSE)
  }
  if (!is.logical(norm_auc) || length(norm_auc) != 1L || is.na(norm_auc) ||
      !norm_auc) {
    stop("Only normalized AUCell scores are currently supported; set `norm_auc = TRUE`.",
         call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  tie_method <- match.arg(tie_method)
  threads <- .cc_aucell_scalar_integer(threads, "threads")
  if (threads > 32L) {
    stop("`threads` cannot exceed 32.", call. = FALSE)
  }
  seed <- .cc_aucell_scalar_integer(seed, "seed", minimum = 0L)
  if (!is.numeric(memory_mib) || length(memory_mib) != 1L ||
      is.na(memory_mib) || !is.finite(memory_mib) || memory_mib <= 0) {
    stop("`memory_mib` must be one positive finite number.", call. = FALSE)
  }

  assay <- SummarizedExperiment::assay(spe, assay_name, withDimnames = FALSE)
  matrix_type <- if (inherits(assay, "Matrix")) {
    DelayedArray::type(DelayedArray::DelayedArray(assay))
  } else {
    NULL
  }
  numeric_assay <- is.numeric(assay) ||
    (inherits(assay, "DelayedArray") &&
       DelayedArray::type(assay) %in% c("integer", "double")) ||
    (!is.null(matrix_type) && matrix_type %in% c("integer", "double"))
  if (!numeric_assay || length(dim(assay)) != 2L ||
      nrow(assay) == 0L || ncol(assay) == 0L) {
    stop("The selected assay must be a non-empty numeric matrix-like object.",
         call. = FALSE)
  }
  features <- rownames(spe)
  if (is.null(features) || anyNA(features) || any(!nzchar(features)) ||
      anyDuplicated(features)) {
    stop("Assay feature identifiers must be unique and non-empty.", call. = FALSE)
  }
  indices <- .cc_aucell_prepare_gene_sets(gene_sets, features)

  if (!is.numeric(auc_max_rank) || length(auc_max_rank) != 1L ||
      is.na(auc_max_rank) || !is.finite(auc_max_rank)) {
    stop("`auc_max_rank` must be one finite number.", call. = FALSE)
  }
  cutoff <- round(auc_max_rank)
  if (cutoff < 2L || cutoff > nrow(assay)) {
    stop("`round(auc_max_rank)` must be between 2 and the assay feature count.",
         call. = FALSE)
  }

  chunk_width <- .cc_aucell_chunk_width(
    nrow(assay), length(indices), threads, memory_mib
  )
  chunk_starts <- .cc_aucell_schedule(ncol(assay), chunk_width)
  scores <- matrix(NA_real_, nrow = ncol(assay), ncol = length(indices),
                   dimnames = list(colnames(spe), names(indices)))
  progress_bar <- if (progress) {
    utils::txtProgressBar(min = 0L, max = length(chunk_starts), style = 3L)
  } else {
    NULL
  }
  if (!is.null(progress_bar)) {
    on.exit(if (!is.null(progress_bar)) close(progress_bar), add = TRUE)
  }

  for (chunk_number in seq_along(chunk_starts)) {
    start <- chunk_starts[[chunk_number]]
    spot_index <- seq.int(start, min(ncol(assay), start + chunk_width - 1L))
    chunk <- as.matrix(assay[, spot_index, drop = FALSE])
    if (!is.double(chunk)) storage.mode(chunk) <- "double"
    if (any(!is.finite(chunk))) {
      stop("The assay contains non-finite values in processing chunk ",
           chunk_number, ".", call. = FALSE)
    }
    chunk_scores <- .Call(
      "C_cpp_aucell_score",
      chunk,
      unname(indices),
      as.integer(cutoff),
      as.double(seed),
      as.integer(spot_index),
      threads,
      tie_method,
      PACKAGE = "CellCharterR"
    )
    if (!is.numeric(chunk_scores) || !is.matrix(chunk_scores) ||
        !identical(dim(chunk_scores), c(length(spot_index), length(indices))) ||
        any(!is.finite(chunk_scores))) {
      stop("The native AUCell engine returned an invalid score matrix.",
           call. = FALSE)
    }
    scores[spot_index, ] <- chunk_scores
    if (!is.null(progress_bar)) utils::setTxtProgressBar(progress_bar, chunk_number)
  }
  if (!is.null(progress_bar)) {
    close(progress_bar)
    progress_bar <- NULL
  }

  SingleCellExperiment::reducedDim(spe, output) <- scores
  spe
}
