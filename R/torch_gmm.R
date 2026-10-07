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
  if (length(n) != 1L || length(p) != 1L || n < 1L || p < 1L)
    stop("`x` must have at least one observation and one feature.")
  list(n = n, p = p, get = function(i) as.matrix(x[i, , drop = FALSE]))
}

.rc_positive_integer <- function(x, name, allow_null = FALSE, allow_zero = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  minimum <- if (allow_zero) 0 else 1
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < minimum || x != as.integer(x))
    stop(sprintf("`%s` must be one %sinteger.", name,
                 if (allow_zero) "non-negative " else "positive "), call. = FALSE)
  as.integer(x)
}

.rc_nonnegative_number <- function(x, name, strictly_positive = FALSE) {
  lower_ok <- if (strictly_positive) x > 0 else x >= 0
  if (length(x) != 1L || is.na(x) || !is.finite(x) || !lower_ok)
    stop(sprintf("`%s` must be one finite %s number.", name,
                 if (strictly_positive) "positive" else "non-negative"), call. = FALSE)
  as.numeric(x)
}

.rc_chunk_starts <- function(n, batch_size) seq.int(1L, n, by = batch_size)
.rc_chunk_index <- function(start, n, batch_size) seq.int(start, min(n, start + batch_size - 1L))

.rc_get <- function(reader, idx) {
  value <- reader$get(idx)
  if (!is.numeric(value) || ncol(value) != reader$p)
    stop("Every block of `x` must be a numeric matrix with a stable feature dimension.", call. = FALSE)
  if (any(!is.finite(value)))
    stop("`x` contains non-finite values; impute or remove them before fitting the GMM.", call. = FALSE)
  value
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
  n <- .rc_positive_integer(n, "n")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  base::set.seed(seed)
  .rc_get(r, sample.int(r$n, min(n, r$n), replace = FALSE))
}

.rc_tensor <- function(x, device) torch::torch_tensor(x, dtype = torch::torch_float32(), device = device)
.rc_eye <- function(p, device) torch::torch_eye(p, dtype = torch::torch_float32(), device = device)

.rc_initial_means <- function(reader, k, x, seed, device) {
  sample <- torch_gmm_sample_rows(x, n = max(4096L, 32L * k), seed = seed)
  if (any(!is.finite(sample)))
    stop("`x` contains non-finite values; impute or remove them before fitting the GMM.", call. = FALSE)
  if (nrow(sample) < k) stop("More GMM components than observations.")
  # Random distinct starts plus a few Lloyd iterations avoid materialising x.
  base::set.seed(seed)
  means <- .rc_tensor(sample[sample.int(nrow(sample), k), , drop = FALSE], device)
  xs <- .rc_tensor(sample, device)
  for (iter in seq_len(8L)) {
    d2 <- torch::torch_stack(lapply(seq_len(k), function(j) {
      torch::torch_sum((xs - means[j, ])$pow(2), dim = 2)
    }), dim = 2)
    lab <- torch::torch_argmin(d2, dim = 2)
    updated <- lapply(seq_len(k), function(j) {
      # R torch exposes index-valued reductions as one-based indices. Keep the
      # comparison one-based here, then convert final public labels to the
      # Python/CellCharter zero-based convention in `torch_gmm_predict()`.
      member <- lab == j
      member_count <- member$sum()$item()
      # A weighted sum preserves the feature dimension even when exactly one
      # sampled observation belongs to a component; squeezed index tensors do
      # not have that guarantee in R torch.
      if (member_count == 0) {
        means[j, ]
      } else {
        (xs * member$unsqueeze(2))$sum(dim = 1) / member_count
      }
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

.rc_estep_stats <- function(reader, batch_size, means, chol, weights, device) {
  k <- dim(means)[1]; p <- dim(means)[2]
  nk <- torch::torch_zeros(k, dtype = torch::torch_float32(), device = device)
  sx <- torch::torch_zeros(c(k, p), dtype = torch::torch_float32(), device = device)
  sxx <- lapply(seq_len(k), function(.) torch::torch_zeros(c(p, p), dtype = torch::torch_float32(), device = device))
  nll <- 0
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    xb <- .rc_tensor(.rc_get(reader, idx), device)
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

.rc_gpu_name <- function(cuda) {
  if (!cuda || !nzchar(Sys.which("nvidia-smi"))) return(NA_character_)
  value <- try(system2("nvidia-smi", c("--query-gpu=name", "--format=csv,noheader"),
                       stdout = TRUE, stderr = FALSE), silent = TRUE)
  if (inherits(value, "try-error") || !length(value)) NA_character_ else value[[1L]]
}

#' Report the Torch execution environment
#'
#' @param device Selected Torch device.
#' @param batch_size Maximum observations per device batch.
#' @param seed Initialization seed.
#' @param num_threads Requested CPU thread count, or `NULL` to retain Torch's setting.
#' @return A named list suitable for checkpoint provenance.
#' @export
torch_gmm_runtime_info <- function(device = "auto", batch_size = 8192L,
                                   seed = 1L, num_threads = NULL) {
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  num_threads <- .rc_positive_integer(num_threads, "num_threads", allow_null = TRUE)
  if (!is.null(num_threads)) torch::torch_set_num_threads(num_threads)
  cuda <- isTRUE(torch::cuda_is_available())
  selected <- torch_gmm_device(device)$type
  list(
    device = selected,
    torch_version = as.character(utils::packageVersion("torch")),
    cuda_available = cuda,
    gpu_name = .rc_gpu_name(cuda),
    batch_size = batch_size,
    seed = seed,
    cpu_threads = torch::torch_get_num_threads()
  )
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
#' @param num_threads Requested CPU thread count, or `NULL` to retain Torch's
#'   current setting.
#' @return A `torch_gmm` model.  Only mixture parameters are retained in memory.
#' @export
torch_gmm_fit <- function(x, n_clusters, batch_size = 8192L, max_iter = 100L,
                          tolerance = 1e-3, covariance_regularization = 1e-6,
                          device = "auto", seed = 1L, num_threads = NULL) {
  reader <- .rc_reader(x)
  k <- .rc_positive_integer(n_clusters, "n_clusters")
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  max_iter <- .rc_positive_integer(max_iter, "max_iter")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  tolerance <- .rc_nonnegative_number(tolerance, "tolerance", strictly_positive = TRUE)
  covariance_regularization <- .rc_nonnegative_number(
    covariance_regularization, "covariance_regularization", strictly_positive = TRUE
  )
  num_threads <- .rc_positive_integer(num_threads, "num_threads", allow_null = TRUE)
  # CPU Torch otherwise may inherit a conservative one-thread configuration.
  # Do not combine a large thread count with outer AutoK parallel workers.
  if (!is.null(num_threads)) torch::torch_set_num_threads(num_threads)
  # K = 1 is needed as the lower neighbour when AutoK evaluates a declared
  # range beginning at K = 2, matching CellCharter v0.2.0.
  if (k < 1L || k > reader$n) stop("`n_clusters` must be between 1 and nrow(x).")
  # Validate every disk-backed block before starting an expensive EM fit. This
  # is a streaming scan and never realizes the complete matrix.
  for (start in .rc_chunk_starts(reader$n, batch_size))
    .rc_get(reader, .rc_chunk_index(start, reader$n, batch_size))

  runtime <- torch_gmm_runtime_info(device, batch_size, seed, num_threads)
  message(sprintf(
    "rCellCharter GMM: device=%s torch=%s cuda=%s gpu=%s batch_size=%d seed=%d cpu_threads=%d",
    runtime$device, runtime$torch_version, runtime$cuda_available,
    ifelse(is.na(runtime$gpu_name), "NA", runtime$gpu_name), runtime$batch_size,
    runtime$seed, runtime$cpu_threads
  ))
  dev <- torch_gmm_device(runtime$device)
  means <- .rc_initial_means(reader, k, x, seed, dev)
  eye <- .rc_eye(reader$p, dev)
  # Identity starts are robust and discarded after the first streamed E-step.
  covariances <- replicate(k, eye, simplify = FALSE)
  chol <- .rc_cholesky(covariances, covariance_regularization, eye)
  weights <- torch::torch_ones(k, dtype = torch::torch_float32(), device = dev) / k
  old_nll <- Inf; converged <- FALSE
  collapsed_history <- vector("list", max_iter)
  for (iter in seq_len(max_iter)) {
    st <- .rc_estep_stats(reader, batch_size, means, chol, weights, dev)
    nk_cpu <- as.numeric(torch::as_array(st$nk$to(device = "cpu")))
    collapsed <- which(!is.finite(nk_cpu) | nk_cpu < 1e-6)
    collapsed_history[[iter]] <- collapsed
    safe_nk <- .rc_tensor(pmax(nk_cpu, 1), dev)
    means <- st$sx / safe_nk$unsqueeze(2)
    if (length(collapsed)) {
      for (j in collapsed) {
        row <- ((as.double(seed) + iter + j - 2) %% reader$n) + 1
        means[j, ] <- .rc_tensor(.rc_get(reader, as.integer(row)), dev)$squeeze(1)
      }
    }
    covariances <- lapply(seq_len(k), function(j) {
      if (j %in% collapsed) {
        # Re-seed collapsed components from a deterministic observation and an
        # identity covariance so repeated fits remain reproducible.
        eye
      } else {
        value <- st$sxx[[j]] / st$nk[j] -
          means[j, ]$unsqueeze(2)$matmul(means[j, ]$unsqueeze(1))
        (value + value$t()) / 2
      }
    })
    weight_values <- pmax(nk_cpu, 1e-6)
    weights <- .rc_tensor(weight_values / sum(weight_values), dev)
    chol <- .rc_cholesky(covariances, covariance_regularization, eye)
    mean_nll <- st$nll / reader$n
    if (is.finite(old_nll) && abs(old_nll - mean_nll) < tolerance) { converged <- TRUE; break }
    old_nll <- mean_nll
  }
  # The loop's E-step evaluates the parameters from before its M-step. Re-run a
  # streamed E-step so the persisted likelihood describes the final model.
  final_stats <- .rc_estep_stats(reader, batch_size, means, chol, weights, dev)
  mean_nll <- final_stats$nll / reader$n
  # Persist ordinary R arrays, not Torch external pointers. Torch tensors do
  # not survive saveRDS/readRDS across sessions, while AutoK relies on durable
  # per-fit model checkpoints and lazy resumption.
  structure(list(
                 means = torch::as_array(means$to(device = "cpu")),
                 covariances = lapply(covariances, function(z) {
                   torch::as_array(z$to(device = "cpu"))
                 }),
                 weights = as.numeric(torch::as_array(weights$to(device = "cpu"))),
                 n_features = reader$p, n_clusters = k,
                 nll = mean_nll, iterations = iter, converged = converged,
                 covariance_regularization = covariance_regularization,
                 collapsed_components = collapsed_history[seq_len(iter)],
                 runtime = runtime), class = "torch_gmm")
}

.rc_model_device <- function(model, device) {
  dev <- torch_gmm_device(device); eye <- .rc_eye(model$n_features, dev)
  means <- .rc_tensor(model$means, dev)
  weights <- .rc_tensor(model$weights, dev)
  covs <- lapply(model$covariances, .rc_tensor, device = dev)
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
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (reader$p != model$n_features) stop("Feature dimensions do not match the fitted model.")
  pars <- .rc_model_device(model, device); out <- vector("list", ceiling(reader$n / batch_size)); b <- 0L
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    b <- b + 1L; xb <- .rc_tensor(.rc_get(reader, idx), pars$device)
    lp <- .rc_logprob(xb, pars$means, pars$chol, pars$weights)
    out[[b]] <- torch::as_array(torch::torch_exp(lp - torch::torch_logsumexp(lp, dim = 2)$unsqueeze(2))$to(device = "cpu"))
  }
  do.call(rbind, out)
}

#' Predict CellCharter domain labels from a Torch GMM
#' @inheritParams torch_gmm_predict_proba
#' @return A zero-based integer label vector matching CellCharter's Python
#'   component convention.
#' @export
torch_gmm_predict <- function(model, newdata, batch_size = 8192L, device = "auto") {
  stopifnot(inherits(model, "torch_gmm"))
  reader <- .rc_reader(newdata)
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (reader$p != model$n_features) stop("Feature dimensions do not match the fitted model.")
  pars <- .rc_model_device(model, device)
  labels <- integer(reader$n)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    lp <- .rc_logprob(.rc_tensor(.rc_get(reader, idx), pars$device),
                      pars$means, pars$chol, pars$weights)
    # Only this batch's K scores exist at once; the complete N-by-K posterior is
    # never collected just to obtain labels.
    labels[idx] <- as.integer(torch::as_array(
      torch::torch_argmax(lp, dim = 2)$to(device = "cpu")
    )) - 1L
  }
  labels
}

#' @export
predict.torch_gmm <- function(object, newdata, ...) torch_gmm_predict(object, newdata, ...)

#' Streamed mean log likelihood
#' @inheritParams torch_gmm_predict_proba
#' @export
torch_gmm_log_likelihood <- function(model, newdata, batch_size = 8192L, device = "auto") {
  stopifnot(inherits(model, "torch_gmm"))
  reader <- .rc_reader(newdata)
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (reader$p != model$n_features) stop("Feature dimensions do not match the fitted model.")
  pars <- .rc_model_device(model, device); value <- 0
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    lp <- .rc_logprob(.rc_tensor(.rc_get(reader, idx), pars$device), pars$means, pars$chol, pars$weights)
    value <- value + torch::torch_logsumexp(lp, dim = 2)$sum()$item()
  }
  value / reader$n
}
