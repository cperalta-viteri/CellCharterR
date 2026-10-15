# Fixed-K clustering is deliberately separate from AutoK. This keeps a request
# for one K from expanding to the neighbouring models required by stability.

.cc_fixed_validate_args <- function(x, reserved, name) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x)))))) {
    stop("`", name, "` must be a named list.", call. = FALSE)
  }
  conflict <- intersect(names(x), reserved)
  if (length(conflict)) {
    stop("`", name, "` must not supply reserved argument(s): ",
         paste(conflict, collapse = ", "), ".", call. = FALSE)
  }
  x
}

.cc_fixed_atomic_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("fixed-k-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically publish fixed-K checkpoint `", path, "`.", call. = FALSE)
  }
  invisible(path)
}

#' Fit one fixed-K model and attach its domains
#'
#' Unlike [cellcharter_cluster()], this function performs exactly one GMM fit
#' at the requested `n_clusters`; it does not fit neighbouring K values for
#' AutoK stability. A completed checkpoint can be resumed after interruption.
#'
#' @section Experimental status:
#' Fixed-K fitting uses the experimental R Torch GMM path by default. For
#' multi-million-observation datasets, consider fitting with Python CellCharter
#' and storing its labels directly in `colData(sce)`.
#'
#' @param sce An aggregated `SingleCellExperiment` or `SpatialExperiment`.
#' @param representation Aggregated representation in `reducedDims(sce)`.
#' @param n_clusters Exact number of mixture components to fit.
#' @param domain_col Output column in `colData(sce)`.
#' @param checkpoint_dir Optional durable checkpoint directory.
#' @param seed Deterministic GMM initialization seed.
#' @param fit_args,predict_args Named additional arguments for `fit_fun` and
#'   `predict_fun`, respectively.
#' @param resume Reuse an exact compatible completed checkpoint when available.
#' @param fit_fun,predict_fun Fitting and prediction hooks. Their defaults are
#'   [torch_gmm_fit()] and [torch_gmm_predict()].
#' @return The updated experiment with factor labels in `domain_col`.
#' @examples
#' spe <- cellcharter_mouse_spleen_tiny(with_results = TRUE)
#' demo_fit <- function(x, n_clusters, seed, ...) list(k = n_clusters, seed = seed)
#' demo_predict <- function(model, newdata, ...) {
#'   (seq_len(nrow(newdata)) - 1L) %% model$k
#' }
#' spe <- cellcharter_cluster_fixed(
#'   spe, representation = "X_cellcharter_tiny", n_clusters = 4L,
#'   domain_col = "fixed_domain", fit_fun = demo_fit,
#'   predict_fun = demo_predict
#' )
#' table(spe$fixed_domain)
#' @export
cellcharter_cluster_fixed <- function(
    sce, representation = "X_cellcharter", n_clusters,
    domain_col = "cellcharter_domain", checkpoint_dir = NULL, seed = 1L,
    fit_args = list(), predict_args = list(), resume = TRUE,
    fit_fun = torch_gmm_fit, predict_fun = torch_gmm_predict) {
  x <- .cc_sce_reduced_dim(sce, representation)
  aggregation <- .cc_sce_require_stage(sce, paste0("aggregation:", representation))
  raw_n_clusters <- n_clusters
  raw_seed <- seed
  n_clusters <- as.integer(n_clusters)
  seed <- as.integer(seed)
  if (!is.numeric(raw_n_clusters) || length(raw_n_clusters) != 1L || is.na(raw_n_clusters) ||
      !is.finite(raw_n_clusters) || raw_n_clusters != n_clusters ||
      n_clusters < 1L ||
      n_clusters > nrow(x)) {
    stop("`n_clusters` must be one integer between 1 and the number of cells.",
         call. = FALSE)
  }
  if (!is.numeric(raw_seed) || length(raw_seed) != 1L || is.na(raw_seed) || !is.finite(raw_seed) ||
      raw_seed != seed || seed < 0L) {
    stop("`seed` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.character(domain_col) || length(domain_col) != 1L || is.na(domain_col) ||
      !nzchar(domain_col)) stop("`domain_col` must be one non-empty string.", call. = FALSE)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("`resume` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.function(fit_fun) || !is.function(predict_fun)) {
    stop("`fit_fun` and `predict_fun` must be functions.", call. = FALSE)
  }
  fit_args <- .cc_fixed_validate_args(
    fit_args, c("x", "n_clusters", "seed"), "fit_args"
  )
  predict_args <- .cc_fixed_validate_args(
    predict_args, c("model", "newdata"), "predict_args"
  )

  identity <- .cc_validate_sce_identity(sce, reduced_dims = representation)
  signature <- digest::digest(list(
    schema = 1L, cell_order_sha256 = identity$cell_order_sha256,
    representation = representation, dimensions = dim(x),
    aggregation = aggregation, n_clusters = n_clusters, seed = seed,
    fit_args = fit_args, predict_args = predict_args,
    implementation = cellcharter_implementation_digest(
      functions = list(fit = fit_fun, predict = predict_fun)
    )
  ), algo = "sha256")
  checkpoint_path <- NULL
  artifact <- NULL
  if (!is.null(checkpoint_dir)) {
    if (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L ||
        is.na(checkpoint_dir) || !nzchar(checkpoint_dir)) {
      stop("`checkpoint_dir` must be NULL or one non-empty path.", call. = FALSE)
    }
    checkpoint_path <- file.path(checkpoint_dir, paste0("fixed_k", n_clusters, ".rds"))
    if (resume && file.exists(checkpoint_path)) {
      candidate <- try(readRDS(checkpoint_path), silent = TRUE)
      if (!inherits(candidate, "try-error") && identical(candidate$signature, signature) &&
          identical(candidate$status, "complete")) artifact <- candidate
    }
  }

  reused <- !is.null(artifact)
  if (!reused) {
    model <- do.call(fit_fun, c(list(
      x = x, n_clusters = n_clusters, seed = seed
    ), fit_args))
    labels <- do.call(predict_fun, c(list(model = model, newdata = x), predict_args))
    labels <- as.integer(labels)
    if (length(labels) != ncol(sce) || anyNA(labels) ||
        any(labels < 0L | labels >= n_clusters)) {
      stop("Fixed-K prediction must return one zero-based label in [0, K-1] per cell.",
           call. = FALSE)
    }
    artifact <- list(
      status = "complete", signature = signature, model = model, labels = labels,
      n_clusters = n_clusters, seed = seed, representation = representation,
      cell_order_sha256 = identity$cell_order_sha256,
      created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    if (!is.null(checkpoint_path)) .cc_fixed_atomic_save(artifact, checkpoint_path)
  }

  labels <- artifact$labels
  if (length(labels) != ncol(sce) || anyNA(labels) ||
      any(labels < 0L | labels >= n_clusters)) {
    stop("Fixed-K checkpoint labels do not align one-to-one with experiment cells.",
         call. = FALSE)
  }
  sce[[domain_col]] <- factor(labels, levels = seq.int(0L, n_clusters - 1L))
  metadata <- .cc_sce_metadata(sce)
  metadata$fixed_k[[domain_col]] <- list(
    n_clusters = n_clusters, seed = seed, representation = representation,
    signature = signature, resumed = reused,
    checkpoint = if (is.null(checkpoint_path)) NULL else normalizePath(checkpoint_path),
    model = if (is.null(checkpoint_path)) artifact$model else NULL,
    labels = if (is.null(checkpoint_path)) labels else NULL
  )
  S4Vectors::metadata(sce)$CellCharterR <- metadata
  .cc_sce_record_stage(
    sce, paste0("clustering:", domain_col),
    list(representation = representation, n_clusters = n_clusters,
         domain_col = domain_col, seed = seed, fixed_k = TRUE,
         resumed = reused, checkpoint = metadata$fixed_k[[domain_col]]$checkpoint)
  )
}
