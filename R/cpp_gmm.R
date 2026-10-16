#' Fit a full-covariance GMM with the compiled CPU backend
#'
#' This experimental accelerated CPU backend uses double precision, centered
#' two-pass covariance updates, and a private deterministic C++ random-number
#' generator. Its K-means++ initialization is not RNG-identical to the pinned
#' sklearn initializer, so fitted solutions require empirical parity checks and
#' are not described as exact seeded CellCharter reproductions. Its result
#' has class `torch_gmm`, so the existing prediction and checkpoint interfaces
#' can consume it. The input is currently an in-memory numeric matrix.
#' @inheritParams torch_gmm_fit
#' @return A serializable `torch_gmm` model.
#' @useDynLib CellCharterR, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
cpp_gmm_fit <- function(x, n_clusters, max_iter = 100L, tolerance = 1e-3,
                        covariance_regularization = 1e-6, seed = 1L,
                        num_threads = 1L) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || !nrow(x) || !ncol(x))
    stop("`x` must be a non-empty numeric matrix.", call. = FALSE)
  n_clusters <- .rc_positive_integer(n_clusters, "n_clusters")
  max_iter <- .rc_positive_integer(max_iter, "max_iter")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  tolerance <- .rc_nonnegative_number(tolerance, "tolerance", strictly_positive = TRUE)
  covariance_regularization <- .rc_nonnegative_number(
    covariance_regularization, "covariance_regularization", strictly_positive = TRUE)
  num_threads <- .rc_positive_integer(num_threads, "num_threads")
  if (num_threads > 32L) { warning("Capping `num_threads` at 32.", call. = FALSE); num_threads <- 32L }
  if (n_clusters > nrow(x)) stop("`n_clusters` cannot exceed nrow(x).", call. = FALSE)
  model <- .Call("C_cpp_gmm_fit", x, n_clusters, as.double(seed), max_iter,
                 tolerance, covariance_regularization, num_threads, PACKAGE = "CellCharterR")
  model <- .cpp_gmm_compatibility(model)
  model$runtime <- list(device = "cpu", backend = "compiled-cpp-double",
                        seed = seed, cpu_threads = num_threads)
  class(model) <- "torch_gmm"
  model
}

#' Fit independent compiled GMM jobs concurrently
#'
#' Worker threads share one immutable copy of `x`; fitting kernels do not call
#' the R API or use R's global RNG. R objects are assembled after all workers
#' join. This is intended for the independent fits in an AutoK campaign.
#' @param n_clusters Integer vector of component counts, one per fit.
#' @param seeds Integer seeds, recycled when length one.
#' @param workers Maximum concurrent native worker threads.
#' @param threads_per_fit Native row-shard threads within each independent fit.
#' @inheritParams cpp_gmm_fit
#' @return A list of serializable `torch_gmm` models.
#' @export
cpp_gmm_fit_many <- function(x, n_clusters, seeds = seq_along(n_clusters), workers = 1L,
                             threads_per_fit = 1L,
                             max_iter = 100L, tolerance = 1e-3,
                             covariance_regularization = 1e-6) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || !nrow(x) || !ncol(x) || any(!is.finite(x)))
    stop("`x` must be a finite, non-empty numeric matrix.", call. = FALSE)
  if (!length(n_clusters)) stop("`n_clusters` must not be empty.", call. = FALSE)
  ks <- vapply(n_clusters, .rc_positive_integer, integer(1), name = "n_clusters")
  if (any(ks > nrow(x))) stop("Every value in `n_clusters` must be <= nrow(x).", call. = FALSE)
  if (length(seeds) == 1L) seeds <- rep(seeds, length(ks))
  if (length(seeds) != length(ks)) stop("`seeds` must have length one or length(n_clusters).", call. = FALSE)
  seeds <- vapply(seeds, .rc_positive_integer, integer(1), name = "seeds", allow_zero = TRUE)
  workers <- .rc_positive_integer(workers, "workers")
  if (workers > 32L) { warning("Capping `workers` at 32.", call. = FALSE); workers <- 32L }
  threads_per_fit <- .rc_positive_integer(threads_per_fit, "threads_per_fit")
  if (workers * threads_per_fit > 32L) {
    threads_per_fit <- max(1L, 32L %/% workers)
    warning(sprintf("Capping `threads_per_fit` at %d so the native topology does not exceed 32 threads.", threads_per_fit), call. = FALSE)
  }
  max_iter <- .rc_positive_integer(max_iter, "max_iter")
  tolerance <- .rc_nonnegative_number(tolerance, "tolerance", strictly_positive = TRUE)
  covariance_regularization <- .rc_nonnegative_number(
    covariance_regularization, "covariance_regularization", strictly_positive = TRUE)
  models <- .Call("C_cpp_gmm_fit_many", x, ks, as.double(seeds), max_iter,
                  tolerance, covariance_regularization, workers, threads_per_fit,
                  PACKAGE = "CellCharterR")
  for (i in seq_along(models)) {
    models[[i]] <- .cpp_gmm_compatibility(models[[i]])
    models[[i]]$runtime <- list(device = "cpu", backend = "compiled-cpp-double",
                                seed = seeds[[i]], cpu_threads = threads_per_fit,
                                outer_workers = workers,
                                threads_per_fit = threads_per_fit)
    class(models[[i]]) <- "torch_gmm"
  }
  names(models) <- paste0("k", ks)
  models
}

.cpp_gmm_compatibility <- function(model) {
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
    branch = "full-data C++ KMeans; native threads parallelize independent fits")
  model$reference_semantics <- paste(
    "CellCharter 0.2.0-style full-data two-pass centered full-covariance EM;",
    "compiled double-precision CPU backend with a distinct deterministic",
    "C++ K-means++/Lloyd initializer (not sklearn RNG-identical).")
  model
}

#' Predict from a compiled GMM without requiring Torch
#' @param model A model returned by [cpp_gmm_fit()].
#' @param newdata Numeric observation-by-feature matrix.
#' @param batch_size Maximum rows scored at once.
#' @return Zero-based component labels.
#' @export
cpp_gmm_predict <- function(model, newdata, batch_size = 8192L) {
  probability <- cpp_gmm_predict_proba(model, newdata, batch_size)
  max.col(probability, ties.method = "first") - 1L
}

#' @rdname cpp_gmm_predict
#' @return `cpp_gmm_predict_proba()` returns posterior probabilities.
#' @export
cpp_gmm_predict_proba <- function(model, newdata, batch_size = 8192L) {
  if (!inherits(model, "torch_gmm")) stop("`model` must inherit from torch_gmm.", call. = FALSE)
  x <- as.matrix(newdata)
  if (!is.numeric(x) || length(dim(x)) != 2L || any(!is.finite(x))) {
    stop("`newdata` must be a finite numeric matrix-like object.", call. = FALSE)
  }
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (ncol(x) != model$n_features) stop("Feature dimensions do not match the fitted model.", call. = FALSE)
  out <- matrix(NA_real_, nrow(x), model$n_clusters)
  if (!nrow(x)) return(out)
  for (start in seq.int(1L, nrow(x), by = batch_size)) {
    index <- start:min(nrow(x), start + batch_size - 1L)
    score <- vapply(seq_len(model$n_clusters), function(component) {
      covariance <- model$covariances[[component]] +
        diag(model$covariance_regularization, model$n_features)
      factor <- chol(covariance)
      centered <- sweep(x[index, , drop = FALSE], 2L, model$means[component, ], "-")
      solved <- forwardsolve(t(factor), t(centered))
      log(model$weights[[component]]) - 0.5 * (
        model$n_features * log(2 * pi) + 2 * sum(log(diag(factor))) +
          colSums(solved * solved))
    }, numeric(length(index)))
    score <- matrix(score, nrow = length(index), ncol = model$n_clusters)
    maximum <- apply(score, 1L, max)
    probability <- exp(score - maximum)
    out[index, ] <- probability / rowSums(probability)
  }
  out
}

#' Mean log likelihood from a compiled GMM without requiring Torch
#' @inheritParams cpp_gmm_predict
#' @return Mean per-observation log likelihood.
#' @export
cpp_gmm_log_likelihood <- function(model, newdata, batch_size = 8192L) {
  if (!inherits(model, "torch_gmm")) stop("`model` must inherit from torch_gmm.", call. = FALSE)
  x <- as.matrix(newdata)
  if (!is.numeric(x) || length(dim(x)) != 2L || any(!is.finite(x))) {
    stop("`newdata` must be a finite numeric matrix-like object.", call. = FALSE)
  }
  if (!nrow(x)) stop("`newdata` must contain at least one observation.", call. = FALSE)
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (ncol(x) != model$n_features) stop("Feature dimensions do not match the fitted model.", call. = FALSE)
  total <- 0
  for (start in seq.int(1L, nrow(x), by = batch_size)) {
    index <- start:min(nrow(x), start + batch_size - 1L)
    score <- vapply(seq_len(model$n_clusters), function(component) {
      factor <- chol(model$covariances[[component]] +
                       diag(model$covariance_regularization, model$n_features))
      centered <- sweep(x[index, , drop = FALSE], 2L, model$means[component, ], "-")
      solved <- forwardsolve(t(factor), t(centered))
      log(model$weights[[component]]) - 0.5 * (model$n_features * log(2*pi) +
        2 * sum(log(diag(factor))) + colSums(solved * solved))
    }, numeric(length(index)))
    score <- matrix(score, nrow = length(index), ncol = model$n_clusters)
    maximum <- apply(score, 1L, max)
    total <- total + sum(maximum + log(rowSums(exp(score - maximum))))
  }
  total / nrow(x)
}
