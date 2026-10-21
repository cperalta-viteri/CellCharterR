#' Fit a full-covariance GMM with the compiled CPU backend
#'
#' This accelerated CPU backend uses double precision, a globally centered
#' full-covariance EM update, and a private deterministic C++ random-number
#' generator. During each fit it caches one observation-by-component posterior
#' matrix and accumulates covariance sufficient statistics in row blocks with
#' the linked BLAS. Gaussian scoring precomputes the precision-Cholesky factor
#' once per component and EM iteration. Its K-means++ initialization is not
#' RNG-identical to the pinned
#' sklearn initializer, so fitted solutions require empirical parity checks and
#' are not described as exact seeded CellCharter reproductions. Its result
#' inherits from both `cpp_gmm` and `torch_gmm`: the dedicated `cpp_gmm` class
#' keeps generic prediction on the double-precision, Torch-free scoring path,
#' while the parent class preserves existing checkpoint compatibility. The
#' input is currently an in-memory numeric matrix. Long fits poll for user
#' interrupts and join all native worker threads before returning to R.
#' @param x Observation-by-feature numeric matrix-like object. The compiled
#'   backend materialises it once as a dense double matrix and makes one shared
#'   native row-major copy plus a centered copy shared by concurrent fits.
#'   Each active fit additionally stores one observation-by-component
#'   double-precision posterior matrix.
#' @param num_threads Positive number of native row-shard threads used by this
#'   fit. Values above 32 are capped with a warning. Unlike the Torch backend,
#'   `NULL` is not accepted.
#' @inheritParams torch_gmm_fit
#' @return A serializable `cpp_gmm` model inheriting from `torch_gmm`.
#' @useDynLib CellCharterR, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
cpp_gmm_fit <- function(x, n_clusters, max_iter = 100L, tolerance = 1e-3,
                        covariance_regularization = 1e-6, seed = 1L,
                        num_threads = 1L) {
  n_clusters <- .rc_positive_integer(n_clusters, "n_clusters")
  max_iter <- .rc_positive_integer(max_iter, "max_iter")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  tolerance <- .rc_nonnegative_number(tolerance, "tolerance", strictly_positive = TRUE)
  covariance_regularization <- .rc_nonnegative_number(
    covariance_regularization, "covariance_regularization", strictly_positive = TRUE)
  num_threads <- .rc_positive_integer(num_threads, "num_threads")
  if (num_threads > 32L) { warning("Capping `num_threads` at 32.", call. = FALSE); num_threads <- 32L }
  .cpp_gmm_require_single_thread_blas(1L, num_threads)
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || !nrow(x) || !ncol(x))
    stop("`x` must be a non-empty numeric matrix.", call. = FALSE)
  if (n_clusters > nrow(x)) stop("`n_clusters` cannot exceed nrow(x).", call. = FALSE)
  # Use the one-job form of the monitored native engine. This preserves the
  # exact same fit kernel and output while allowing R's main thread to poll for
  # Ctrl+C and cooperatively stop/join the row-shard workers.
  model <- .Call(
    "C_cpp_gmm_fit_many",
    x,
    n_clusters,
    as.double(seed),
    max_iter,
    tolerance,
    covariance_regularization,
    1L,
    num_threads,
    FALSE,
    PACKAGE = "CellCharterR"
  )[[1L]]
  model <- .cpp_gmm_compatibility(model)
  model$runtime <- list(device = "cpu", backend = "compiled-cpp-double",
                        seed = seed, cpu_threads = num_threads)
  class(model) <- c("cpp_gmm", "torch_gmm")
  model
}

#' Fit independent compiled GMM jobs concurrently
#'
#' Worker threads share one immutable copy of `x`; fitting kernels do not call
#' the R API or use R's global RNG. R objects are assembled after all workers
#' join. This is intended for the independent fits in an AutoK campaign.
#' @param n_clusters Non-empty integer vector of component counts, one per fit;
#'   every value must be positive and no greater than the number of observations.
#' @param seeds Integer seeds, recycled only when length one. Otherwise the
#'   vector must have one seed per requested fit.
#' @param workers Positive maximum number of independent fits run concurrently.
#'   Values above 32 are capped. This controls between-fit parallelism. When
#'   the effective native topology uses more than one thread, start R with
#'   `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`,
#'   `BLIS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, and `GOTO_NUM_THREADS` all
#'   set to `1`; the function fails early otherwise to prevent nested BLAS
#'   oversubscription.
#' @param threads_per_fit Positive native row-shard thread count within each
#'   independent fit. The product `workers * threads_per_fit` is capped at 32,
#'   with `threads_per_fit` reduced when necessary.
#' @param progress Whether to report native input preparation, completed fits,
#'   and a heartbeat during long periods without a completed fit.
#' @inheritParams cpp_gmm_fit
#' @return A list in input order containing one serializable `cpp_gmm` model per
#'   `(n_clusters, seed)` pair. Models inherit from `torch_gmm` for checkpoint
#'   compatibility but prediction remains on the Torch-free compiled path.
#' @export
cpp_gmm_fit_many <- function(x, n_clusters, seeds = seq_along(n_clusters), workers = 1L,
                             threads_per_fit = 1L,
                             max_iter = 100L, tolerance = 1e-3,
                             covariance_regularization = 1e-6,
                             progress = interactive()) {
  if (!length(n_clusters)) stop("`n_clusters` must not be empty.", call. = FALSE)
  workers <- .rc_positive_integer(workers, "workers")
  if (workers > 32L) { warning("Capping `workers` at 32.", call. = FALSE); workers <- 32L }
  threads_per_fit <- .rc_positive_integer(threads_per_fit, "threads_per_fit")
  if (as.double(workers) * as.double(threads_per_fit) > 32) {
    threads_per_fit <- max(1L, 32L %/% workers)
    warning(sprintf("Capping `threads_per_fit` at %d so the native topology does not exceed 32 threads.", threads_per_fit), call. = FALSE)
  }
  # The native engine owns fit- and row-level parallelism.  Its covariance
  # kernels call the linked BLAS, so a second BLAS thread pool would multiply
  # the requested topology (for example 8 x 4 x 32 OpenBLAS threads).  These
  # controls must be established before R starts; setting them here can be too
  # late once a BLAS implementation has initialized its worker pool.
  .cpp_gmm_require_single_thread_blas(workers, threads_per_fit)
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || !nrow(x) || !ncol(x) || any(!is.finite(x)))
    stop("`x` must be a finite, non-empty numeric matrix.", call. = FALSE)
  ks <- vapply(n_clusters, .rc_positive_integer, integer(1), name = "n_clusters")
  if (any(ks > nrow(x))) stop("Every value in `n_clusters` must be <= nrow(x).", call. = FALSE)
  if (length(seeds) == 1L) seeds <- rep(seeds, length(ks))
  if (length(seeds) != length(ks)) stop("`seeds` must have length one or length(n_clusters).", call. = FALSE)
  seeds <- vapply(seeds, .rc_positive_integer, integer(1), name = "seeds", allow_zero = TRUE)
  max_iter <- .rc_positive_integer(max_iter, "max_iter")
  tolerance <- .rc_nonnegative_number(tolerance, "tolerance", strictly_positive = TRUE)
  covariance_regularization <- .rc_nonnegative_number(
    covariance_regularization, "covariance_regularization", strictly_positive = TRUE)
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  models <- .Call("C_cpp_gmm_fit_many", x, ks, as.double(seeds), max_iter,
                  tolerance, covariance_regularization, workers, threads_per_fit,
                  progress,
                  PACKAGE = "CellCharterR")
  for (i in seq_along(models)) {
    models[[i]] <- .cpp_gmm_compatibility(models[[i]])
    models[[i]]$runtime <- list(device = "cpu", backend = "compiled-cpp-double",
                                seed = seeds[[i]], cpu_threads = threads_per_fit,
                                outer_workers = workers,
                                threads_per_fit = threads_per_fit)
    class(models[[i]]) <- c("cpp_gmm", "torch_gmm")
  }
  names(models) <- paste0("k", ks)
  models
}

.cpp_gmm_blas_environment <- function() {
  variables <- c(
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "GOTO_NUM_THREADS"
  )
  values <- Sys.getenv(variables, unset = "")
  stats::setNames(as.character(values), variables)
}

.cpp_gmm_blas_contract <- function() {
  versions <- extSoftVersion()
  version_value <- function(name) {
    if (name %in% names(versions)) unname(as.character(versions[[name]])) else NA_character_
  }
  list(
    library = version_value("BLAS"),
    lapack = version_value("LAPACK"),
    environment = as.list(.cpp_gmm_blas_environment())
  )
}

.cpp_gmm_require_single_thread_blas <- function(workers, threads_per_fit) {
  if (as.double(workers) * as.double(threads_per_fit) <= 1) {
    return(invisible(.cpp_gmm_blas_contract()))
  }
  values <- .cpp_gmm_blas_environment()
  unsafe <- names(values)[values != "1"]
  if (length(unsafe)) {
    assignments <- paste0(unsafe, "=1", collapse = " ")
    stop(
      "Parallel compiled GMM fitting requires a single-threaded BLAS to avoid ",
      "nested oversubscription. Start a fresh R process with all BLAS/OpenMP ",
      "controls set to 1 before loading CellCharterR. Missing or unsafe: ",
      paste(unsafe, collapse = ", "), ". For example: `", assignments,
      " Rscript your-script.R`.",
      call. = FALSE
    )
  }
  invisible(.cpp_gmm_blas_contract())
}

.cpp_gmm_compatibility <- function(model) {
  if (is.null(model$nll_history) || !length(model$nll_history) ||
      !is.finite(utils::tail(model$nll_history, 1L))) {
    stop("Compiled GMM model lacks a finite monitored NLL history.", call. = FALSE)
  }
  # CellCharter 0.2.0 ranks repeated fits using PyCave's `nll_`, which is the
  # last pre-M-step training metric. Keep the independently rescored final NLL
  # in `model$nll`, but expose the reference selection quantity explicitly.
  model$selection_nll <- as.numeric(utils::tail(model$nll_history, 1L))
  model$selection_nll_semantics <- "last monitored mean NLL before M-step"
  attempts <- max(0L, as.integer(round(log10(
    model$covariance_regularization / model$requested_covariance_regularization))))
  model$regularization_attempts <-
    model$requested_covariance_regularization * 10^(0:attempts)
  model$regularization_max_attempts <- 8L
  model$regularization_restart_errors <- character(attempts)
  model$covariance_storage <- paste(
    "Unregularized centered covariances; covariance_regularization",
    "is added exactly once when constructing the fitted distribution.")
  model$collapse_threshold <- 1e-6
  model$initialization <- list(
    method = "cpp-kmeans++-lloyd-hard", counts = model$initialization_counts,
    collapsed_components = model$initialization_collapsed_components,
    branch = paste(
      "full-data C++ KMeans; native threads parallelize independent fits,",
      "Lloyd assignments, and initialization covariance"
    ))
  model$reference_semantics <- paste(
    "CellCharter 0.2.0-style full-data centered full-covariance EM with",
    "precision-Cholesky scoring, cached posteriors, and blocked BLAS",
    "covariance sufficient statistics;",
    "compiled double-precision CPU backend with a distinct deterministic",
    "C++ K-means++/Lloyd initializer (not sklearn RNG-identical).")
  model
}

#' Predict from a compiled GMM without requiring Torch
#' @param model A model returned by [cpp_gmm_fit()].
#' @param newdata Numeric observation-by-feature matrix or `DelayedMatrix`.
#' @param batch_size Maximum rows read and scored at once. Label prediction
#'   retains only one `batch_size`-by-K score block in memory; posterior
#'   prediction explicitly returns the complete observation-by-K matrix.
#' @return Zero-based component labels.
#' @export
cpp_gmm_predict <- function(model, newdata, batch_size = 8192L) {
  scoring <- .cpp_gmm_scoring_input(model, newdata, batch_size)
  labels <- integer(scoring$reader$n)
  if (!length(labels)) return(labels)
  prepared <- .cpp_gmm_prepare_scoring(model)
  for (start in seq.int(1L, scoring$reader$n, by = scoring$batch_size)) {
    index <- start:min(scoring$reader$n, start + scoring$batch_size - 1L)
    # Retain only this batch's N_batch-by-K scores. Label prediction never
    # constructs or normalizes the complete N-by-K posterior matrix.
    block <- .cpp_gmm_get_block(scoring$reader, index)
    score <- .cpp_gmm_log_scores(block, prepared)
    labels[index] <- max.col(score, ties.method = "first") - 1L
  }
  labels
}

#' Predict labels from a compiled GMM
#'
#' Dispatches compiled models to the double-precision compiled-model scorer
#' without loading Torch or materialising a complete posterior matrix.
#' @param object A model returned by [cpp_gmm_fit()].
#' @param ... Additional arguments passed to [cpp_gmm_predict()].
#' @return Zero-based component labels.
#' @rdname cpp_gmm_predict
#' @export
predict.cpp_gmm <- function(object, newdata, ...) {
  cpp_gmm_predict(object, newdata, ...)
}

.cpp_gmm_scoring_input <- function(model, newdata, batch_size) {
  if (!inherits(model, "torch_gmm")) {
    stop("`model` must inherit from torch_gmm.", call. = FALSE)
  }
  if (!is.matrix(newdata) && !inherits(newdata, "DelayedMatrix")) {
    stop("`newdata` must be a numeric matrix or DelayedMatrix.", call. = FALSE)
  }
  dimensions <- dim(newdata)
  if (length(dimensions) != 2L || is.na(dimensions[[1L]]) ||
      is.na(dimensions[[2L]])) {
    stop("`newdata` must have two valid dimensions.", call. = FALSE)
  }
  if (is.matrix(newdata) &&
      (!is.numeric(newdata) || any(!is.finite(newdata)))) {
    stop("`newdata` must be a finite numeric matrix-like object.", call. = FALSE)
  }
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (dimensions[[2L]] != model$n_features) {
    stop("Feature dimensions do not match the fitted model.", call. = FALSE)
  }
  reader <- list(
    n = as.integer(dimensions[[1L]]),
    p = as.integer(dimensions[[2L]]),
    get = function(index) {
      if (is.matrix(newdata)) {
        newdata[index, , drop = FALSE]
      } else {
        as.matrix(newdata[index, , drop = FALSE])
      }
    }
  )
  list(reader = reader, batch_size = batch_size)
}

.cpp_gmm_get_block <- function(reader, index) {
  block <- reader$get(index)
  if (!is.matrix(block) || !is.numeric(block) || ncol(block) != reader$p ||
      nrow(block) != length(index) || any(!is.finite(block))) {
    stop("`newdata` must yield finite numeric matrix blocks.", call. = FALSE)
  }
  block
}

.cpp_gmm_prepare_scoring <- function(model) {
  factors <- lapply(model$covariances, function(covariance) {
    chol(covariance + diag(model$covariance_regularization, model$n_features))
  })
  constants <- vapply(seq_len(model$n_clusters), function(component) {
    log(model$weights[[component]]) - 0.5 * (
      model$n_features * log(2 * pi) + 2 * sum(log(diag(factors[[component]])))
    )
  }, numeric(1L))
  list(
    means = model$means,
    factors = factors,
    constants = constants,
    n_clusters = model$n_clusters
  )
}

.cpp_gmm_log_scores <- function(x, prepared) {
  score <- vapply(seq_len(prepared$n_clusters), function(component) {
    centered <- sweep(x, 2L, prepared$means[component, ], "-")
    solved <- forwardsolve(t(prepared$factors[[component]]), t(centered))
    prepared$constants[[component]] - 0.5 * colSums(solved * solved)
  }, numeric(nrow(x)))
  matrix(score, nrow = nrow(x), ncol = prepared$n_clusters)
}

#' @rdname cpp_gmm_predict
#' @return `cpp_gmm_predict_proba()` returns posterior probabilities.
#' @export
cpp_gmm_predict_proba <- function(model, newdata, batch_size = 8192L) {
  scoring <- .cpp_gmm_scoring_input(model, newdata, batch_size)
  out <- matrix(NA_real_, scoring$reader$n, model$n_clusters)
  if (!scoring$reader$n) return(out)
  prepared <- .cpp_gmm_prepare_scoring(model)
  for (start in seq.int(1L, scoring$reader$n, by = scoring$batch_size)) {
    index <- start:min(scoring$reader$n, start + scoring$batch_size - 1L)
    score <- .cpp_gmm_log_scores(
      .cpp_gmm_get_block(scoring$reader, index), prepared
    )
    maximum <- apply(score, 1L, max)
    probability <- exp(score - maximum)
    out[index, ] <- probability / rowSums(probability)
  }
  out
}

#' Mean log likelihood from a compiled GMM without requiring Torch
#' @inheritParams cpp_gmm_predict
#' @param batch_size Positive number of rows read and scored at once. Only one
#'   `batch_size`-by-K score block and a scalar running sum are retained; no
#'   complete observation-by-component matrix is allocated.
#' @return Mean per-observation log likelihood.
#' @export
cpp_gmm_log_likelihood <- function(model, newdata, batch_size = 8192L) {
  scoring <- .cpp_gmm_scoring_input(model, newdata, batch_size)
  if (!scoring$reader$n) {
    stop("`newdata` must contain at least one observation.", call. = FALSE)
  }
  prepared <- .cpp_gmm_prepare_scoring(model)
  total <- 0
  for (start in seq.int(1L, scoring$reader$n, by = scoring$batch_size)) {
    index <- start:min(scoring$reader$n, start + scoring$batch_size - 1L)
    score <- .cpp_gmm_log_scores(
      .cpp_gmm_get_block(scoring$reader, index), prepared
    )
    maximum <- apply(score, 1L, max)
    total <- total + sum(maximum + log(rowSums(exp(score - maximum))))
  }
  total / scoring$reader$n
}
