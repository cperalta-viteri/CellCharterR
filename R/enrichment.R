#' CellCharter observed-versus-expected group enrichment
#'
#' Reproduces `cellcharter.gr.enrichment()` from CellCharter v0.2.0. For
#' cell-type enrichment, `group` contains cell types and `label` contains
#' spatial-domain labels. The observed value is the within-domain cell-type
#' proportion and the expected value is the cohort-wide cell-type proportion.
#'
#' @param group Group annotation, typically cell type.
#' @param label Label annotation, typically spatial domain.
#' @param log2_fold_change Return log2 fold changes when `TRUE`; otherwise
#'   return fold changes.
#' @return A list containing `observed`, `expected`, `enrichment`, and `params`.
#' @export
cellcharter_enrichment <- function(group, label, log2_fold_change = TRUE) {
  group <- .cc_complete_factor(group, "group")
  label <- .cc_complete_factor(label, "label")
  if (length(group) != length(label)) {
    stop("`group` and `label` must have the same length.", call. = FALSE)
  }
  if (length(group) == 0L) {
    stop("`group` and `label` must not be empty.", call. = FALSE)
  }
  if (!is.logical(log2_fold_change) || length(log2_fold_change) != 1L ||
      is.na(log2_fold_change)) {
    stop("`log2_fold_change` must be TRUE or FALSE.", call. = FALSE)
  }

  # CellCharter first calculates group proportions within every label and then
  # transposes the table so groups are rows and labels are columns.
  counts <- unclass(table(label, group))
  label_totals <- rowSums(counts)
  observed <- t(sweep(counts, 1L, label_totals, FUN = "/"))

  expected_proportion <- as.numeric(table(group)) / length(group)
  names(expected_proportion) <- levels(group)
  expected <- matrix(
    expected_proportion,
    nrow = length(expected_proportion), ncol = nlevels(label),
    dimnames = list(levels(group), levels(label))
  )
  dimnames(observed) <- dimnames(expected)

  enrichment <- observed / expected
  if (log2_fold_change) enrichment <- log2(enrichment)

  # pandas.DataFrame.fillna(enrichment.min()) fills missing values using the
  # minimum of the corresponding column. Infinite fold changes are retained.
  for (column in seq_len(ncol(enrichment))) {
    missing <- is.nan(enrichment[, column])
    available <- enrichment[!missing, column]
    if (any(missing) && length(available) > 0L) {
      enrichment[missing, column] <- min(available)
    }
  }

  list(
    observed = observed,
    expected = expected,
    enrichment = enrichment,
    params = list(
      log2_fold_change = log2_fold_change,
      reference = "CellCharter 0.2.0"
    )
  )
}

#' Asymmetric CellCharter neighborhood enrichment
#'
#' Computes directed inter-domain neighborhood proportions and the analytical
#' expectation used by CellCharter v0.2.0. With `only_inter = TRUE`, same-domain
#' links are removed before calculating source-domain proportions. The pinned
#' analytical expectation is the target domain's cell fraction and is not
#' renormalized after same-domain links are removed.
#'
#' @param adjacency Symmetric sparse cell-by-cell connectivity matrix.
#' @param cluster Spatial-domain annotation in adjacency row order.
#' @param log2_fold_change Use `log2(observed / expected)` instead of the
#'   v0.2.0 default `observed - expected`.
#' @param only_inter Remove same-domain links and return an `NA` diagonal.
#' @return A list containing `observed`, `expected`, `enrichment`, and `params`.
#' @export
cellcharter_nhood_enrichment <- function(adjacency, cluster,
                                         log2_fold_change = FALSE,
                                         only_inter = TRUE) {
  checked <- .cc_nhood_inputs(adjacency, cluster)
  .cc_nhood_from_entries(
    entries = checked$entries,
    cluster = checked$cluster,
    log2_fold_change = log2_fold_change,
    only_inter = only_inter
  )
}

#' Differential CellCharter neighborhood enrichment
#'
#' Computes condition-specific asymmetric neighborhood enrichments and their
#' difference, then estimates directional empirical p-values by permuting the
#' condition assignment at the library level. Per-library sufficient statistics
#' are computed once, avoiding a full graph scan for every permutation.
#'
#' `permutation_mode = "cellcharter_v0.2.0"` intentionally reproduces the
#' pinned Python implementation: it draws a vector of library names with
#' replacement, uses the first `n_condition1` draws as the selected set, and
#' calculates unselected minus selected enrichment. `"shuffle"` instead makes
#' the conventional size-preserving, without-replacement assignment and
#' calculates condition 1 minus condition 2.
#'
#' @param adjacency Symmetric sparse cell-by-cell connectivity matrix.
#' @param cluster Spatial-domain annotation in adjacency row order.
#' @param condition Condition annotation in adjacency row order.
#' @param library Physical sample/library annotation in adjacency row order.
#'   Every library must belong to exactly one condition, and graph edges may
#'   not join different libraries.
#' @param condition_groups Two conditions in the desired subtraction order.
#' @param n_perms Positive number of library-label permutations.
#' @param seed Deterministic R random seed.
#' @param log2_fold_change Forwarded to neighborhood enrichment.
#' @param only_inter Remove same-domain links before analysis.
#' @param permutation_mode Pinned CellCharter parity mode or a conventional
#'   without-replacement shuffle.
#' @return A list with condition-specific `observed`, `expected`, and
#'   `enrichment` matrices, the differential `effect`, `permutation_counts`,
#'   empirical `pvalue`, and complete `params`.
#' @export
cellcharter_diff_nhood_enrichment <- function(
    adjacency, cluster, condition, library, condition_groups = NULL,
    n_perms = 100L, seed = 0L, log2_fold_change = FALSE,
    only_inter = TRUE,
    permutation_mode = c("cellcharter_v0.2.0", "shuffle")) {
  permutation_mode <- match.arg(permutation_mode)
  checked <- .cc_nhood_inputs(adjacency, cluster)
  cluster <- checked$cluster
  condition <- .cc_complete_factor(condition, "condition")
  library <- .cc_complete_factor(library, "library")
  n <- nrow(adjacency)
  if (length(condition) != n || length(library) != n) {
    stop("`condition` and `library` must have one value per adjacency row.", call. = FALSE)
  }
  if (!is.numeric(n_perms) || length(n_perms) != 1L || is.na(n_perms) ||
      !is.finite(n_perms) || n_perms < 1 || n_perms != floor(n_perms)) {
    stop("`n_perms` must be one positive whole number.", call. = FALSE)
  }
  n_perms <- as.integer(n_perms)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    stop("`seed` must be one finite whole number.", call. = FALSE)
  }
  if (!is.logical(log2_fold_change) || length(log2_fold_change) != 1L ||
      is.na(log2_fold_change) || !is.logical(only_inter) ||
      length(only_inter) != 1L || is.na(only_inter)) {
    stop("`log2_fold_change` and `only_inter` must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(condition_groups)) {
    present <- levels(condition)[tabulate(as.integer(condition), nlevels(condition)) > 0L]
    condition_groups <- present
  }
  if (length(condition_groups) != 2L || anyNA(condition_groups) ||
      anyDuplicated(condition_groups) ||
      !all(condition_groups %in% levels(condition))) {
    stop("`condition_groups` must identify two distinct observed conditions.", call. = FALSE)
  }
  condition_groups <- as.character(condition_groups)
  keep <- as.character(condition) %in% condition_groups
  if (!all(condition_groups %in% as.character(condition[keep]))) {
    stop("Both `condition_groups` must contain observations.", call. = FALSE)
  }

  # A library is the permutation unit and therefore cannot be split across
  # conditions. This also catches accidental reuse of sample identifiers.
  library_condition <- split(
    as.character(condition[keep]), droplevels(library[keep])
  )
  if (any(vapply(library_condition, function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("Each library must belong to exactly one condition.", call. = FALSE)
  }

  retained_rows <- which(keep)
  adjacency_subset <- adjacency[retained_rows, retained_rows, drop = FALSE]
  subset_checked <- .cc_nhood_inputs(adjacency_subset, cluster[retained_rows])
  stats <- .cc_nhood_library_stats(
    entries = subset_checked$entries,
    cluster = subset_checked$cluster,
    library = droplevels(library[retained_rows]),
    only_inter = only_inter
  )

  library_to_condition <- vapply(
    split(
      as.character(condition[retained_rows]),
      droplevels(library[retained_rows])
    ),
    function(x) unique(x)[[1L]], character(1)
  )
  libraries1 <- names(library_to_condition)[library_to_condition == condition_groups[[1L]]]
  libraries2 <- names(library_to_condition)[library_to_condition == condition_groups[[2L]]]
  if (length(libraries1) == 0L || length(libraries2) == 0L) {
    stop("Each requested condition must contain at least one library.", call. = FALSE)
  }

  result1 <- .cc_nhood_from_library_stats(
    stats, libraries1, log2_fold_change, only_inter
  )
  result2 <- .cc_nhood_from_library_stats(
    stats, libraries2, log2_fold_change, only_inter
  )
  effect <- result1$enrichment - result2$enrichment

  all_libraries <- c(libraries1, libraries2)
  n_condition1 <- length(libraries1)
  permutation_counts <- matrix(
    0L, nrow = nrow(effect), ncol = ncol(effect), dimnames = dimnames(effect)
  )

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) old_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  base::set.seed(seed)

  for (permutation in seq_len(n_perms)) {
    if (identical(permutation_mode, "cellcharter_v0.2.0")) {
      # Exact v0.2.0 semantics, including replacement and 0-minus-1 direction.
      draws <- sample(all_libraries, length(all_libraries), replace = TRUE)
      selected <- unique(draws[seq_len(n_condition1)])
      unselected <- setdiff(all_libraries, selected)
      permuted_effect <-
        .cc_nhood_from_library_stats(stats, unselected, log2_fold_change, only_inter)$enrichment -
        .cc_nhood_from_library_stats(stats, selected, log2_fold_change, only_inter)$enrichment
    } else {
      selected <- sample(all_libraries, n_condition1, replace = FALSE)
      unselected <- setdiff(all_libraries, selected)
      permuted_effect <-
        .cc_nhood_from_library_stats(stats, selected, log2_fold_change, only_inter)$enrichment -
        .cc_nhood_from_library_stats(stats, unselected, log2_fold_change, only_inter)$enrichment
    }

    positive <- !is.na(effect) & effect > 0
    negative <- !is.na(effect) & effect < 0
    # NumPy comparisons against NaN are FALSE. Because v0.2.0 computes
    # 1 - mean(null < observed) (or > for negative observations), a NaN null
    # value contributes to the complementary tail count used here.
    permutation_counts[positive] <- permutation_counts[positive] +
      as.integer(
        is.na(permuted_effect[positive]) |
          permuted_effect[positive] >= effect[positive]
      )
    permutation_counts[negative] <- permutation_counts[negative] +
      as.integer(
        is.na(permuted_effect[negative]) |
          permuted_effect[negative] <= effect[negative]
      )
  }

  # This is the directional empirical tail definition in v0.2.0: no +1
  # correction, p=0 for an observed zero, and p=0 on the masked diagonal.
  pvalue <- permutation_counts / n_perms
  pvalue[is.na(effect) | effect == 0] <- 0

  condition_results <- stats::setNames(list(result1, result2), condition_groups)
  list(
    observed = lapply(condition_results, `[[`, "observed"),
    expected = lapply(condition_results, `[[`, "expected"),
    enrichment = lapply(condition_results, `[[`, "enrichment"),
    effect = effect,
    permutation_counts = permutation_counts,
    pvalue = pvalue,
    params = list(
      condition_groups = condition_groups,
      library_unit = "sample",
      library_counts = stats::setNames(c(length(libraries1), length(libraries2)), condition_groups),
      n_perms = n_perms,
      seed = seed,
      log2_fold_change = log2_fold_change,
      only_inter = only_inter,
      permutation_mode = permutation_mode,
      pvalue_correction = "none",
      reference = "CellCharter 0.2.0"
    )
  )
}

.cc_complete_factor <- function(x, name) {
  if (is.null(x) || is.data.frame(x) || is.matrix(x)) {
    stop(sprintf("`%s` must be an atomic vector or factor.", name), call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("`%s` must not contain missing values.", name), call. = FALSE)
  }
  if (is.factor(x)) return(x)
  factor(as.character(x), levels = unique(as.character(x)))
}

.cc_nhood_inputs <- function(adjacency, cluster) {
  if (!inherits(adjacency, "sparseMatrix") || nrow(adjacency) != ncol(adjacency)) {
    stop("`adjacency` must be a square sparse Matrix object.", call. = FALSE)
  }
  adjacency <- methods::as(adjacency, "generalMatrix")
  entries <- Matrix::summary(adjacency)
  if (nrow(entries) > 0L &&
      (anyNA(entries$x) || any(!is.finite(entries$x)) || any(entries$x < 0))) {
    stop("`adjacency` must contain finite, non-negative weights.", call. = FALSE)
  }
  adjacency <- Matrix::drop0(adjacency)
  if (!isTRUE(Matrix::isSymmetric(adjacency))) {
    stop("`adjacency` must be symmetric.", call. = FALSE)
  }
  cluster <- .cc_complete_factor(cluster, "cluster")
  if (length(cluster) != nrow(adjacency)) {
    stop("`cluster` must have one value per adjacency row.", call. = FALSE)
  }
  entries <- Matrix::summary(adjacency)
  list(adjacency = adjacency, entries = entries, cluster = cluster)
}

.cc_nhood_from_entries <- function(entries, cluster, log2_fold_change, only_inter) {
  k <- nlevels(cluster)
  cluster_index <- as.integer(cluster)
  keep <- if (only_inter && nrow(entries) > 0L) {
    cluster_index[entries$i] != cluster_index[entries$j]
  } else {
    rep.int(TRUE, nrow(entries))
  }
  links <- matrix(0, nrow = k, ncol = k)
  if (any(keep)) {
    linear_index <- cluster_index[entries$i[keep]] +
      (cluster_index[entries$j[keep]] - 1L) * k
    aggregated <- rowsum(as.numeric(entries$x[keep]), linear_index, reorder = FALSE)
    links[as.integer(rownames(aggregated))] <- aggregated[, 1L]
  }
  counts <- tabulate(cluster_index, nbins = k)
  .cc_nhood_from_counts(
    links, counts, levels(cluster), log2_fold_change, only_inter
  )
}

.cc_nhood_from_counts <- function(links, counts, cluster_levels,
                                  log2_fold_change, only_inter) {
  source_totals <- rowSums(links)
  observed <- sweep(links, 1L, source_totals, FUN = "/")

  # v0.2.0 constructs a constant mean-degree vector. For asymmetric analysis
  # this simplifies to n_target / N whenever the graph has at least one edge.
  if (sum(links) == 0) {
    expected_proportion <- rep.int(NaN, length(counts))
  } else {
    expected_proportion <- counts / sum(counts)
  }
  expected <- matrix(
    expected_proportion,
    nrow = length(counts), ncol = length(counts), byrow = TRUE
  )
  dimnames(observed) <- dimnames(expected) <- list(cluster_levels, cluster_levels)
  enrichment <- if (log2_fold_change) log2(observed / expected) else observed - expected
  if (only_inter) {
    diag(observed) <- NA_real_
    diag(expected) <- NA_real_
    diag(enrichment) <- NA_real_
  }
  list(
    observed = observed,
    expected = expected,
    enrichment = enrichment,
    params = list(
      log2_fold_change = log2_fold_change,
      only_inter = only_inter,
      symmetric = FALSE,
      reference = "CellCharter 0.2.0"
    )
  )
}

.cc_nhood_library_stats <- function(entries, cluster, library, only_inter) {
  library <- droplevels(library)
  cluster_index <- as.integer(cluster)
  library_index <- as.integer(library)
  k <- nlevels(cluster)
  n_libraries <- nlevels(library)
  if (nrow(entries) > 0L &&
      any(library_index[entries$i] != library_index[entries$j])) {
    stop("`adjacency` contains edges between physical libraries.", call. = FALSE)
  }
  keep <- if (only_inter && nrow(entries) > 0L) {
    cluster_index[entries$i] != cluster_index[entries$j]
  } else {
    rep.int(TRUE, nrow(entries))
  }
  links <- array(0, dim = c(k, k, n_libraries))
  if (any(keep)) {
    linear_index <- cluster_index[entries$i[keep]] +
      (cluster_index[entries$j[keep]] - 1L) * k +
      (library_index[entries$i[keep]] - 1L) * k * k
    aggregated <- rowsum(as.numeric(entries$x[keep]), linear_index, reorder = FALSE)
    links[as.integer(rownames(aggregated))] <- aggregated[, 1L]
  }
  counts <- matrix(0L, nrow = k, ncol = n_libraries)
  count_index <- cluster_index + (library_index - 1L) * k
  counts_tab <- tabulate(count_index, nbins = k * n_libraries)
  counts[] <- counts_tab
  dimnames(links) <- list(levels(cluster), levels(cluster), levels(library))
  dimnames(counts) <- list(levels(cluster), levels(library))
  list(links = links, counts = counts, cluster_levels = levels(cluster))
}

.cc_nhood_from_library_stats <- function(stats, libraries,
                                          log2_fold_change, only_inter) {
  library_index <- match(libraries, dimnames(stats$links)[[3L]])
  if (length(library_index) == 0L || anyNA(library_index)) {
    stop("A permutation produced an empty or unknown library group.", call. = FALSE)
  }
  links <- apply(stats$links[, , library_index, drop = FALSE], c(1L, 2L), sum)
  counts <- rowSums(stats$counts[, library_index, drop = FALSE])
  .cc_nhood_from_counts(
    links, counts, stats$cluster_levels, log2_fold_change, only_inter
  )
}
