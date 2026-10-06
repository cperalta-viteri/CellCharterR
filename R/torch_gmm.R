#' Select the Torch device used by rCellCharter
#'
#' CUDA is used when it is visible to R Torch; otherwise the exact same batched
#' implementation runs on CPU.  This is deliberately separate from HDF5 I/O:
#' each chunk is read from disk on CPU, copied once to the device, and released.
#' @param device `"auto"`, `"cpu"`, or `"cuda"`.
#' @export
torch_gmm_device <- function(device = c("auto", "cuda", "cpu")) {
  device <- match.arg(device)
  cuda <- isTRUE(torch::cuda_is_available())
  if (device == "cuda" && !cuda) stop("CUDA was requested but is not available to R Torch.")
  if (device == "auto") device <- if (cuda) "cuda" else "cpu"
  torch::torch_device(device)
}

.rc_reader <- function(x) {
  if (!is.matrix(x) && !inherits(x, "DelayedMatrix"))
    stop("`x` must be a matrix or DelayedMatrix with observations in rows.")
  n <- nrow(x); p <- ncol(x)
  list(n = n, p = p, get = function(i) as.matrix(x[i, , drop = FALSE]))
}

#' Sample rows from a disk-backed feature matrix
#'
#' The sampled matrix is intentionally the only sizeable object retained on the
#' R heap during initialisation. It is used for k-means-like mean initialisation.
#' @param x observation-by-feature matrix or DelayedMatrix.
#' @param n number of rows to sample.
#' @param seed random seed.
#' @export
torch_gmm_sample_rows <- function(x, n = 10000L, seed = 1L) {
  r <- .rc_reader(x)
  set.seed(seed)
  r$get(sample.int(r$n, min(as.integer(n), r$n), replace = FALSE))
}

.rc_chunks <- function(n, batch_size) split(seq_len(n), ceiling(seq_len(n) / batch_size))
.rc_tensor <- function(x, device) torch::torch_tensor(x, dtype = torch::torch_float32(), device = device)
.rc_eye <- function(p, device) torch::torch_eye(p, dtype = torch::torch_float32(), device = device)

.rc_initial_means <- function(reader, k, init_rows, seed, device) {
  sample <- torch_gmm_sample_rows(init_rows, n = max(4096L, 32L * k), seed = seed)
  if (nrow(sample) < k) stop("More GMM components than observations.")
  # Random distinct starts plus a few Lloyd iterations avoid materialising x.
  set.seed(seed)
  means <- .rc_tensor(sample[sample.int(nrow(sample), k), , drop = FALSE], device)
  xs <- .rc_tensor(sample, device)
  for (iter in seq_len(8L)) {
    d2 <- torch::torch_stack(lapply(seq_len(k), function(j) {
      torch::torch_sum((xs - means[j, ])$pow(2), dim = 2)
    }), dim = 2)
    lab <- torch::torch_argmin(d2, dim = 2)
    updated <- lapply(seq_len(k), function(j) {
      member <- torch::torch_nonzero(lab == (j - 1L))$squeeze()
      # `as_array()` is a Torch R function (not a tensor method).  This small
      # CPU transfer only checks for an empty k-means component at start-up.
      if (length(torch::as_array(member$to(device = "cpu"))) == 0L) means[j, ] else xs[member, ]$mean(dim = 1)
    })
    means <- torch::torch_stack(updated, dim = 1)
  }
  means
}

.rc_logprob <- function(x, means, chol, weights) {
  k <- dim(means)[1]
  p <- dim(means)[2]
  norm <- p * log(2 * pi)
  pieces <- lapply(seq_len(k), function(j) {
    centered <- x - means[j, ]
    # R Torch exports the modern linalg wrappers without a `torch_` prefix.
    solved <- torch::linalg_solve_triangular(chol[[j]], centered$t(), upper = FALSE)
    maha <- torch::torch_sum(solved$pow(2), dim = 1)
    logdet <- 2 * torch::torch_log(torch::torch_diag(chol[[j]]))$sum()$item()
    -0.5 * (norm + logdet + maha) + torch::torch_log(weights[j])
  })
  torch::torch_stack(pieces, dim = 2)
}

.rc_cholesky <- function(covariances, reg, eye) {
  lapply(covariances, function(s) {
    # A failed factorisation is repaired locally rather than rerunning all EM.
    current <- reg
    for (attempt in seq_len(8L)) {
      out <- try(torch::linalg_cholesky(s + current * eye), silent = TRUE)
      if (!inherits(out, "try-error")) return(out)
      current <- current * 10
    }
    stop("Covariance remained non-positive-definite after regularisation.")
  })
}

.rc_estep_stats <- function(reader, chunks, means, chol, weights, device) {
  k <- dim(means)[1]; p <- dim(means)[2]
  nk <- torch::torch_zeros(k, dtype = torch::torch_float32(), device = device)
  sx <- torch::torch_zeros(c(k, p), dtype = torch::torch_float32(), device = device)
  sxx <- lapply(seq_len(k), function(.) torch::torch_zeros(c(p, p), dtype = torch::torch_float32(), device = device))
  nll <- 0
  for (idx in chunks) {
    xb <- .rc_tensor(reader$get(idx), device)
    lp <- .rc_logprob(xb, means, chol, weights)
    lse <- torch::torch_logsumexp(lp, dim = 2)
    resp <- torch::torch_exp(lp - lse$unsqueeze(2))
    nk <- nk + resp$sum(dim = 1)
    sx <- sx + resp$t()$matmul(xb)
    for (j in seq_len(k)) {
      w <- resp[, j]$unsqueeze(2)
      sxx[[j]] <- sxx[[j]] + (xb * w)$t()$matmul(xb)
    }
    nll <- nll - lse$sum()$item()
  }
  list(nk = nk, sx = sx, sxx = sxx, nll = nll)
}

#' Fit a batched full-covariance GMM with R Torch
#'
#' This is EM, not gradient-based neural-network training. Thus Torch performs
#' the dense linear algebra and optional GPU execution directly; Luz is not put
#' in the critical path, where its optimiser/training-loop abstraction would add
#' overhead without contributing to closed-form E/M updates.
#'
#' @param x observation-by-feature matrix or DelayedMatrix. Rows are read in
#'   chunks and never materialised jointly.
#' @param n_clusters number of mixture components.
#' @param batch_size maximum number of observations transferred per E-step.
#' @param max_iter maximum EM iterations.
#' @param tolerance stop when per-observation NLL change is below this value.
#' @param covariance_regularization diagonal ridge added before Cholesky.
#' @param device Torch device selection.
#' @param seed initialisation seed.
#' @return A `torch_gmm` model.  Only mixture parameters are retained in memory.
#' @export
torch_gmm_fit <- function(x, n_clusters, batch_size = 8192L, max_iter = 100L,
                          tolerance = 1e-3, covariance_regularization = 1e-6,
                          device = "auto", seed = 1L, num_threads = NULL) {
  reader <- .rc_reader(x); k <- as.integer(n_clusters)
  # CPU Torch otherwise may inherit a conservative one-thread configuration.
  # Do not combine a large thread count with outer AutoK parallel workers.
  if (!is.null(num_threads)) torch::torch_set_num_threads(as.integer(num_threads))
  if (k < 2L || k > reader$n) stop("`n_clusters` must be between 2 and nrow(x).")
  dev <- torch_gmm_device(device); chunks <- .rc_chunks(reader$n, as.integer(batch_size))
  means <- .rc_initial_means(reader, k, x, seed, dev)
  eye <- .rc_eye(reader$p, dev)
  # Identity starts are robust and discarded after the first streamed E-step.
  covariances <- replicate(k, eye, simplify = FALSE)
  chol <- .rc_cholesky(covariances, covariance_regularization, eye)
  weights <- torch::torch_ones(k, dtype = torch::torch_float32(), device = dev) / k
  old_nll <- Inf; converged <- FALSE
  for (iter in seq_len(as.integer(max_iter))) {
    st <- .rc_estep_stats(reader, chunks, means, chol, weights, dev)
    weights <- st$nk / reader$n
    means <- st$sx / st$nk$unsqueeze(2)
    covariances <- lapply(seq_len(k), function(j) {
      st$sxx[[j]] / st$nk[j] - means[j, ]$unsqueeze(2)$matmul(means[j, ]$unsqueeze(1))
    })
    chol <- .rc_cholesky(covariances, covariance_regularization, eye)
    mean_nll <- st$nll / reader$n
    if (is.finite(old_nll) && abs(old_nll - mean_nll) < tolerance) { converged <- TRUE; break }
    old_nll <- mean_nll
  }
  structure(list(means = means$to(device = "cpu"), covariances = lapply(covariances, function(z) z$to(device = "cpu")),
                 weights = weights$to(device = "cpu"), n_features = reader$p, n_clusters = k,
                 nll = mean_nll, iterations = iter, converged = converged,
                 covariance_regularization = covariance_regularization), class = "torch_gmm")
}

.rc_model_device <- function(model, device) {
  dev <- torch_gmm_device(device); eye <- .rc_eye(model$n_features, dev)
  means <- model$means$to(device = dev); weights <- model$weights$to(device = dev)
  covs <- lapply(model$covariances, function(z) z$to(device = dev))
  list(device = dev, means = means, weights = weights,
       chol = .rc_cholesky(covs, model$covariance_regularization, eye))
}

#' Posterior GMM probabilities in bounded-size batches
#' @param model a fitted `torch_gmm`.
#' @param newdata observation-by-feature matrix or DelayedMatrix.
#' @param batch_size disk-to-device chunk size.
#' @param device Torch device selection.
#' @export
torch_gmm_predict_proba <- function(model, newdata, batch_size = 8192L, device = "auto") {
  stopifnot(inherits(model, "torch_gmm")); reader <- .rc_reader(newdata)
  if (reader$p != model$n_features) stop("Feature dimensions do not match the fitted model.")
  pars <- .rc_model_device(model, device); out <- vector("list", ceiling(reader$n / batch_size)); b <- 0L
  for (idx in .rc_chunks(reader$n, as.integer(batch_size))) {
    b <- b + 1L; xb <- .rc_tensor(reader$get(idx), pars$device)
    lp <- .rc_logprob(xb, pars$means, pars$chol, pars$weights)
    out[[b]] <- torch::as_array(torch::torch_exp(lp - torch::torch_logsumexp(lp, dim = 2)$unsqueeze(2))$to(device = "cpu"))
  }
  do.call(rbind, out)
}

#' Predict CellCharter domain labels from a Torch GMM
#' @inheritParams torch_gmm_predict_proba
#' @export
torch_gmm_predict <- function(model, newdata, batch_size = 8192L, device = "auto") {
  max.col(torch_gmm_predict_proba(model, newdata, batch_size, device), ties.method = "first") - 1L
}

#' @export
predict.torch_gmm <- function(object, newdata, ...) torch_gmm_predict(object, newdata, ...)

#' Streamed mean negative log likelihood
#' @inheritParams torch_gmm_predict_proba
#' @export
torch_gmm_log_likelihood <- function(model, newdata, batch_size = 8192L, device = "auto") {
  reader <- .rc_reader(newdata); pars <- .rc_model_device(model, device); value <- 0
  for (idx in .rc_chunks(reader$n, as.integer(batch_size))) {
    lp <- .rc_logprob(.rc_tensor(reader$get(idx), pars$device), pars$means, pars$chol, pars$weights)
    value <- value + torch::torch_logsumexp(lp, dim = 2)$sum()$item()
  }
  value / reader$n
}
