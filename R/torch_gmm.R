#' Select the Torch device used by CellCharterR
#'
#' CUDA is used when it is visible to R Torch; otherwise the exact same batched
#' implementation runs on CPU.  This is deliberately separate from HDF5 I/O:
#' each chunk is read from disk on CPU, copied once to the device, and released.
#' @param device `"auto"`, `"cpu"`, or `"cuda"`.
#' @return A `torch_device` selecting the resolved CPU or CUDA device.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' cpu_device <- torch_gmm_device("cpu")
#' cpu_device$type
#' @export
torch_gmm_device <- function(device = c("auto", "cuda", "cpu")) {
  .cc_require_torch()
  device <- match.arg(device)
  cuda <- isTRUE(torch::cuda_is_available())
  if (device == "cuda" && !cuda) stop("CUDA was requested but is not available to R Torch.")
  if (device == "auto") device <- if (cuda) "cuda" else "cpu"
  torch::torch_device(device)
}

.cc_require_torch <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(
      paste0(
        "The experimental Torch GMM backend requires the optional `torch` ",
        "package. Install it with `install.packages(\"torch\")`, or use ",
        "the CPU C++ backend with `backend = \"cpp\"`."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_reader <- function(x) {
  if (inherits(x, "cellcharter_torch_cache")) {
    valid <- is.list(x) && identical(x$format_version, 2L) &&
      isTRUE(x$validated) && inherits(x$tensor, "torch_tensor") &&
      identical(as.integer(dim(x$tensor)), c(as.integer(x$n), as.integer(x$p))) &&
      identical(x$tensor$device$type, "cpu") &&
      identical(as.character(x$tensor$dtype), "Float")
    if (!isTRUE(valid)) {
      stop("Malformed CellCharterR representation cache; rebuild it from the source representation.",
           call. = FALSE)
    }
    return(list(
      n = x$n, p = x$p, trusted = TRUE, tensor = x$tensor,
      get = function(i) torch::as_array(x$tensor[i, , drop = FALSE])
    ))
  }
  if (!is.matrix(x) && !inherits(x, "DelayedMatrix"))
    stop("`x` must be a matrix or DelayedMatrix with observations in rows.")
  n <- nrow(x); p <- ncol(x)
  if (length(n) != 1L || length(p) != 1L || n < 1L || p < 1L)
    stop("`x` must have at least one observation and one feature.")
  list(n = n, p = p, trusted = FALSE,
       get = function(i) as.matrix(x[i, , drop = FALSE]))
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
  if (!isTRUE(reader$trusted) && any(!is.finite(value)))
    stop("`x` contains non-finite values; impute or remove them before fitting the GMM.", call. = FALSE)
  # R matrices are usually double precision, whereas the Torch backend is
  # deliberately float32 for CellCharter/PyCave parity. Detect overflow before
  # conversion so it is reported as an input problem rather than as repeated
  # Cholesky failures later in EM.
  float32_max <- 3.4028234663852886e38
  if (!isTRUE(reader$trusted) && any(abs(value) > float32_max)) {
    observed <- max(abs(value))
    stop(sprintf(
      paste0(
        "`x` contains magnitude %.6g, outside the finite float32 range ",
        "(maximum %.6g) required by the Torch GMM; centre or rescale the ",
        "representation before fitting."
      ),
      observed, float32_max
    ), call. = FALSE)
  }
  value
}

.rc_ram_cache <- function(x, batch_size = 8192L) {
  source <- .rc_reader(x)
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  cached <- torch::torch_empty(
    c(source$n, source$p), dtype = torch::torch_float32(), device = "cpu"
  )
  for (start in .rc_chunk_starts(source$n, batch_size)) {
    idx <- .rc_chunk_index(start, source$n, batch_size)
    block <- .rc_get(source, idx)
    block_tensor <- torch::torch_tensor(block, dtype = torch::torch_float32())
    if (!isTRUE(torch::torch_isfinite(block_tensor)$all()$item())) {
      stop("`x` contains values outside the finite float32 range required by the GMM.",
           call. = FALSE)
    }
    # Integer-vector indexing creates an advanced-indexing copy in R torch.
    # narrow() is a writable contiguous view and also preserves a 1 x p tail.
    cached$narrow(1L, start, length(idx))$copy_(block_tensor)
  }
  structure(
    list(format_version = 2L, validated = TRUE,
         n = source$n, p = source$p, tensor = cached),
    class = "cellcharter_torch_cache"
  )
}

.rc_device_block <- function(reader, idx, device) {
  if (!is.null(reader$tensor)) {
    return(reader$tensor[idx, , drop = FALSE]$to(device = device))
  }
  .rc_tensor(.rc_get(reader, idx), device)
}

#' Sample rows from a disk-backed feature matrix
#'
#' The sampled matrix is intentionally the only sizeable object retained on the
#' R heap during initialisation. It is used for k-means-like mean initialisation.
#' @param x observation-by-feature matrix or DelayedMatrix.
#' @param n number of rows to sample.
#' @param seed random seed.
#' @return A numeric matrix containing at most `n` sampled rows.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' x <- matrix(seq_len(20), nrow = 10)
#' identical(
#'   torch_gmm_sample_rows(x, n = 3, seed = 7),
#'   torch_gmm_sample_rows(x, n = 3, seed = 7)
#' )
#' @export
torch_gmm_sample_rows <- function(x, n = 10000L, seed = 1L) {
  r <- .rc_reader(x)
  n <- .rc_positive_integer(n, "n")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  withr::with_seed(
    seed,
    .rc_get(r, sample.int(r$n, min(n, r$n), replace = FALSE))
  )
}

.rc_tensor <- function(x, device) torch::torch_tensor(x, dtype = torch::torch_float32(), device = device)
.rc_eye <- function(p, device) torch::torch_eye(p, dtype = torch::torch_float32(), device = device)

.rc_kmeans_distances <- function(x, centers) {
  value <- vapply(seq_len(nrow(centers)), function(j) {
    rowSums(sweep(x, 2L, centers[j, ], FUN = "-")^2)
  }, numeric(nrow(x)))
  matrix(value, nrow = nrow(x), ncol = nrow(centers))
}

.rc_kmeans_tolerance <- function(reader, batch_size, relative_tolerance = 1e-4) {
  seen <- 0
  running_mean <- numeric(reader$p)
  running_m2 <- numeric(reader$p)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    block <- .rc_get(reader, .rc_chunk_index(start, reader$n, batch_size))
    block_n <- nrow(block)
    block_mean <- colMeans(block)
    block_m2 <- colSums(sweep(block, 2L, block_mean, FUN = "-")^2)
    if (seen == 0) {
      running_mean <- block_mean
      running_m2 <- block_m2
      seen <- block_n
    } else {
      total <- seen + block_n
      delta <- block_mean - running_mean
      running_m2 <- running_m2 + block_m2 + delta^2 * seen * block_n / total
      running_mean <- running_mean + delta * block_n / total
      seen <- total
    }
  }
  mean(running_m2 / seen) * relative_tolerance
}

.rc_kmeans_plus_plus <- function(reader, k, batch_size, seed) {
  withr::with_seed(seed, .rc_kmeans_plus_plus_impl(reader, k, batch_size, seed))
}

.rc_kmeans_plus_plus_impl <- function(reader, k, batch_size, seed) {
  centers <- matrix(NA_real_, nrow = k, ncol = reader$p)
  centers[1L, ] <- .rc_get(reader, sample.int(reader$n, 1L))[1L, ]
  if (k == 1L) return(centers)

  # This is the bounded-memory form of sklearn's greedy k-means++ seeding:
  # D2 candidates and their potentials are recomputed by streaming passes.
  n_local_trials <- 2L + floor(log(k))
  for (component in 2:k) {
    active <- centers[seq_len(component - 1L), , drop = FALSE]
    total <- 0
    for (start in .rc_chunk_starts(reader$n, batch_size)) {
      block <- .rc_get(reader, .rc_chunk_index(start, reader$n, batch_size))
      total <- total + sum(apply(.rc_kmeans_distances(block, active), 1L, min))
    }
    if (!is.finite(total) || total <= 0) {
      row <- ((as.double(seed) + component - 2) %% reader$n) + 1
      centers[component, ] <- .rc_get(reader, as.integer(row))[1L, ]
      next
    }

    thresholds <- sort(stats::runif(n_local_trials, min = 0, max = total))
    candidates <- matrix(NA_real_, nrow = n_local_trials, ncol = reader$p)
    found <- rep(FALSE, n_local_trials)
    cumulative <- 0
    for (start in .rc_chunk_starts(reader$n, batch_size)) {
      block <- .rc_get(reader, .rc_chunk_index(start, reader$n, batch_size))
      nearest <- apply(.rc_kmeans_distances(block, active), 1L, min)
      within <- cumulative + cumsum(nearest)
      for (candidate in which(!found)) {
        location <- which(within >= thresholds[[candidate]])[1L]
        if (!is.na(location)) {
          candidates[candidate, ] <- block[location, ]
          found[[candidate]] <- TRUE
        }
      }
      cumulative <- cumulative + sum(nearest)
      if (all(found)) break
    }
    if (any(!found)) candidates[!found, ] <- .rc_get(reader, reader$n)[1L, ]

    potentials <- numeric(n_local_trials)
    for (start in .rc_chunk_starts(reader$n, batch_size)) {
      block <- .rc_get(reader, .rc_chunk_index(start, reader$n, batch_size))
      nearest <- apply(.rc_kmeans_distances(block, active), 1L, min)
      for (candidate in seq_len(n_local_trials)) {
        candidate_d2 <- rowSums(sweep(
          block, 2L, candidates[candidate, ], FUN = "-"
        )^2)
        potentials[[candidate]] <- potentials[[candidate]] +
          sum(pmin(nearest, candidate_d2))
      }
    }
    centers[component, ] <- candidates[which.min(potentials), ]
  }
  centers
}

.rc_stream_kmeans <- function(reader, k, batch_size, seed,
                              max_iter = 300L, tolerance = 1e-4) {
  centers <- .rc_kmeans_plus_plus(reader, k, batch_size, seed)
  absolute_tolerance <- .rc_kmeans_tolerance(reader, batch_size, tolerance)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    sums <- matrix(0, nrow = k, ncol = reader$p)
    counts <- numeric(k)
    farthest_d2 <- numeric()
    farthest_rows <- matrix(numeric(), nrow = 0L, ncol = reader$p)
    for (start in .rc_chunk_starts(reader$n, batch_size)) {
      block <- .rc_get(reader, .rc_chunk_index(start, reader$n, batch_size))
      distances <- .rc_kmeans_distances(block, centers)
      labels <- max.col(-distances, ties.method = "first")
      for (component in seq_len(k)) {
        member <- labels == component
        if (any(member)) {
          counts[[component]] <- counts[[component]] + sum(member)
          sums[component, ] <- sums[component, ] + colSums(block[member, , drop = FALSE])
        }
      }
      nearest <- distances[cbind(seq_len(nrow(block)), labels)]
      keep <- order(nearest, decreasing = TRUE)[seq_len(min(k, length(nearest)))]
      farthest_d2 <- c(farthest_d2, nearest[keep])
      farthest_rows <- rbind(farthest_rows, block[keep, , drop = FALSE])
      retain <- order(farthest_d2, decreasing = TRUE)[seq_len(min(k, length(farthest_d2)))]
      farthest_d2 <- farthest_d2[retain]
      farthest_rows <- farthest_rows[retain, , drop = FALSE]
    }
    updated <- centers
    nonempty <- which(counts > 0)
    updated[nonempty, ] <- sums[nonempty, , drop = FALSE] / counts[nonempty]
    empty <- which(counts == 0)
    if (length(empty)) {
      for (position in seq_along(empty)) {
        source <- min(position, nrow(farthest_rows))
        updated[empty[[position]], ] <- farthest_rows[source, ]
      }
    }
    shift <- sum((updated - centers)^2)
    centers <- updated
    if (shift <= absolute_tolerance) {
      converged <- TRUE
      break
    }
  }
  list(
    means = centers, iterations = iter, converged = converged,
    method = "R-native streamed greedy k-means++/Lloyd",
    reference_method = "sklearn.cluster.KMeans(n_init=1)",
    rng_note = paste(
      "Base R and sklearn use different RNG consumption and numerical kernels;",
      "component identities and exact centroids are not promised to be bit-identical."
    )
  )
}

.rc_logprob <- function(x, means, chol, weights) {
  k <- dim(means)[1]
  p <- dim(means)[2]
  norm <- p * log(2 * pi)
  # Solve all components in one batched kernel. Keeping log-determinants on the
  # device is important on CUDA: `$item()` here used to introduce one device
  # synchronization per component and per input batch.
  chol_tensor <- torch::torch_stack(chol, dim = 1)
  centered <- x$unsqueeze(2) - means$unsqueeze(1) # observations x K x features
  rhs <- centered$permute(c(2, 3, 1))             # K x features x observations
  solved <- torch::linalg_solve_triangular(chol_tensor, rhs, upper = FALSE)
  maha <- solved$pow(2)$sum(dim = 2)$transpose(1, 2)
  logdet <- 2 * torch::torch_log(
    torch::torch_diagonal(chol_tensor, dim1 = 2, dim2 = 3)
  )$sum(dim = 2)
  -0.5 * (norm + maha + logdet$unsqueeze(1)) +
    torch::torch_log(weights)$unsqueeze(1)
}

.rc_cholesky_failure <- function(message) {
  structure(list(message = message, call = NULL),
            class = c("rc_gmm_cholesky_error", "error", "condition"))
}

.rc_cholesky <- function(covariances, reg, eye) {
  lapply(seq_along(covariances), function(component) {
    out <- try(
      torch::linalg_cholesky(covariances[[component]] + reg * eye),
      silent = TRUE
    )
    if (inherits(out, "try-error")) {
      stop(.rc_cholesky_failure(sprintf(
        "Cholesky failed for component %d at covariance regularization %.17g.",
        component, reg
      )))
    }
    out
  })
}

.rc_estep_stats <- function(reader, batch_size, means, chol, weights, device) {
  k <- dim(means)[1]; p <- dim(means)[2]
  nk <- torch::torch_zeros(k, dtype = torch::torch_float32(), device = device)
  sx <- torch::torch_zeros(c(k, p), dtype = torch::torch_float32(), device = device)
  # Keep the historical semantics: each batch is reduced in float32, while
  # batch totals accumulate in float64. Keeping the accumulator on-device
  # removes per-batch synchronization without changing that precision split.
  nll <- torch::torch_zeros(1L, dtype = torch::torch_float64(), device = device)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    xb <- .rc_device_block(reader, idx, device)
    lp <- .rc_logprob(xb, means, chol, weights)
    lse <- torch::torch_logsumexp(lp, dim = 2)
    resp <- torch::torch_exp(lp - lse$unsqueeze(2))
    nk <- nk + resp$sum(dim = 1)
    sx <- sx + resp$t()$matmul(xb)
    nll$sub_(lse$sum()$to(dtype = torch::torch_float64()))
  }
  # Synchronize CUDA only once per complete data pass, not once per batch.
  list(nk = nk, sx = sx, nll = nll$item())
}

.rc_mstep_means_weights <- function(nk, sx, nk_cpu, collapsed,
                                    collapse_threshold, device) {
  k <- length(nk_cpu)
  means <- torch::torch_zeros_like(sx)
  active <- setdiff(seq_len(k), collapsed)
  for (j in active) means[j, ] <- sx[j, ] / nk[j]

  # Active components use their exact effective membership everywhere. Only a
  # collapsed component receives the explicit floor needed to keep it alive.
  weight_mass <- nk_cpu
  weight_mass[collapsed] <- collapse_threshold
  weights <- .rc_tensor(weight_mass / sum(weight_mass), device)
  list(means = means, weights = weights)
}

.rc_hard_initialize <- function(reader, batch_size, initialized_means,
                                device, eye, collapse_threshold = 1e-6) {
  k <- nrow(initialized_means); p <- reader$p
  centroids <- .rc_tensor(initialized_means, device)
  nk <- torch::torch_zeros(k, dtype = torch::torch_float32(), device = device)
  sx <- torch::torch_zeros(c(k, p), dtype = torch::torch_float32(), device = device)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    block <- .rc_device_block(
      reader, .rc_chunk_index(start, reader$n, batch_size), device
    )
    distances <- torch::torch_stack(lapply(seq_len(k), function(component) {
      torch::torch_sum((block - centroids[component, ])$pow(2), dim = 2)
    }), dim = 2)
    labels <- torch::torch_argmin(distances, dim = 2)
    for (component in seq_len(k)) {
      member <- (labels == component)$to(dtype = torch::torch_float32())
      nk[component] <- nk[component] + member$sum()
      sx[component, ] <- sx[component, ] +
        (block * member$unsqueeze(2))$sum(dim = 1)
    }
  }
  nk_cpu <- as.numeric(torch::as_array(nk$to(device = "cpu")))
  collapsed <- which(!is.finite(nk_cpu) | nk_cpu < collapse_threshold)
  update <- .rc_mstep_means_weights(
    nk, sx, nk_cpu, collapsed, collapse_threshold, device
  )
  means <- update$means
  if (length(collapsed)) {
    for (component in collapsed) {
      row <- ((component - 1L) %% reader$n) + 1L
      means[component, ] <- .rc_tensor(.rc_get(reader, row), device)$squeeze(1)
    }
  }

  covariance_sum <- lapply(seq_len(k), function(.) {
    torch::torch_zeros(c(p, p), dtype = torch::torch_float32(), device = device)
  })
  active <- setdiff(seq_len(k), collapsed)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    block <- .rc_device_block(
      reader, .rc_chunk_index(start, reader$n, batch_size), device
    )
    distances <- torch::torch_stack(lapply(seq_len(k), function(component) {
      torch::torch_sum((block - centroids[component, ])$pow(2), dim = 2)
    }), dim = 2)
    labels <- torch::torch_argmin(distances, dim = 2)
    for (component in active) {
      member <- (labels == component)$to(dtype = torch::torch_float32())
      centered <- block - means[component, ]
      covariance_sum[[component]] <- covariance_sum[[component]] +
        (centered * member$unsqueeze(2))$t()$matmul(centered)
    }
  }
  covariances <- lapply(seq_len(k), function(component) {
    if (component %in% collapsed) return(eye)
    value <- covariance_sum[[component]] / nk[component]
    (value + value$t()) / 2
  })
  list(
    means = means, weights = update$weights, covariances = covariances,
    collapsed = collapsed, counts = nk_cpu
  )
}

.rc_centered_covariances <- function(reader, batch_size, old_means, old_chol,
                                     old_weights, new_means, nk, collapsed,
                                     device, eye) {
  k <- dim(old_means)[1]; p <- dim(old_means)[2]
  covariance_sum <- lapply(seq_len(k), function(.) {
    torch::torch_zeros(c(p, p), dtype = torch::torch_float32(), device = device)
  })
  active <- setdiff(seq_len(k), collapsed)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    xb <- .rc_device_block(reader, idx, device)
    # Replay the pre-update model so both M-step passes use identical
    # responsibilities, matching PyCave's batched full-covariance update.
    lp <- .rc_logprob(xb, old_means, old_chol, old_weights)
    resp <- torch::torch_exp(lp - torch::torch_logsumexp(lp, dim = 2)$unsqueeze(2))
    for (j in active) {
      centered <- xb - new_means[j, ]
      covariance_sum[[j]] <- covariance_sum[[j]] +
        (centered * resp[, j]$unsqueeze(2))$t()$matmul(centered)
    }
  }
  lapply(seq_len(k), function(j) {
    if (j %in% collapsed) return(eye)
    value <- covariance_sum[[j]] / nk[j]
    (value + value$t()) / 2
  })
}

.rc_select_gpu_name <- function(value, visible = "") {
  if (!length(value)) return(NA_character_)
  fields <- strsplit(as.character(value), ",", fixed = TRUE)
  fields <- lapply(fields, trimws)
  valid <- lengths(fields) >= 3L
  if (!any(valid)) return(NA_character_)
  fields <- fields[valid]
  index <- vapply(fields, `[[`, character(1L), 1L)
  uuid <- vapply(fields, `[[`, character(1L), 2L)
  name <- vapply(fields, function(x) paste(x[-c(1L, 2L)], collapse = ","),
                 character(1L))
  selector <- strsplit(as.character(visible %||% ""), ",", fixed = TRUE)[[1L]][1L]
  selector <- trimws(selector)
  # `strsplit("")[[1L]][1L]` is NA on supported R releases.  An unset or
  # empty CUDA_VISIBLE_DEVICES means the first device returned by nvidia-smi,
  # not a literal device named "NA".
  if (!length(selector) || is.na(selector) || !nzchar(selector)) {
    return(name[[1L]])
  }
  selected <- if (grepl("^[0-9]+$", selector)) {
    which(index == selector)
  } else {
    which(uuid == selector | startsWith(uuid, selector) |
            startsWith(selector, uuid))
  }
  if (length(selected)) return(name[[selected[[1L]]]])
  # Never report physical GPU 0 when CUDA exposes a UUID/MIG device that the
  # physical-device query cannot resolve.
  paste0("visible CUDA device ", selector)
}

.rc_gpu_name <- function(cuda) {
  if (!cuda || !nzchar(Sys.which("nvidia-smi"))) return(NA_character_)
  value <- try(system2(
    "nvidia-smi",
    c("--query-gpu=index,uuid,name", "--format=csv,noheader"),
    stdout = TRUE, stderr = FALSE
  ), silent = TRUE)
  if (inherits(value, "try-error")) return(NA_character_)
  .rc_select_gpu_name(value, Sys.getenv("CUDA_VISIBLE_DEVICES", ""))
}

#' Report the Torch execution environment
#'
#' @param device Selected Torch device.
#' @param batch_size Maximum observations per device batch.
#' @param seed Initialization seed.
#' @param num_threads Requested CPU thread count, or `NULL` to retain Torch's setting.
#' @return A named list suitable for checkpoint provenance.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' info <- torch_gmm_runtime_info(
#'   device = "cpu", batch_size = 8, seed = 7, num_threads = 1
#' )
#' info[c("device", "batch_size", "seed")]
#' @export
torch_gmm_runtime_info <- function(device = "auto", batch_size = 8192L,
                                   seed = 1L, num_threads = NULL) {
  .cc_require_torch()
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  num_threads <- .rc_positive_integer(num_threads, "num_threads", allow_null = TRUE)
  if (!is.null(num_threads)) {
    old_num_threads <- torch::torch_get_num_threads()
    on.exit(torch::torch_set_num_threads(old_num_threads), add = TRUE)
    torch::torch_set_num_threads(num_threads)
  }
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
#' The parity target is the versioned CellCharter 0.2.0 CODEX notebook. Its
#' `Cluster` call omits Python's `batch_size`, so CellCharter selects sklearn
#' `KMeans(n_init = 1)`. This implementation uses bounded-memory, R-native
#' greedy k-means++ and streamed Lloyd passes. The algorithm is equivalent in
#' role, but sklearn and base R consume random numbers and perform reductions
#' differently, so exact initial centroids and component identities are not
#' promised. No Python runtime is required.
#'
#' PyCave 3.2.1 tests `number_of_loader_batches == 1` when selecting its
#' so-called batch-training path. The notebook therefore takes the two-pass
#' initialization/EM path with early-stopping patience two. `batch_size` here
#' only bounds storage and device transfers; streamed passes reproduce that
#' full-loader behavior. PyCave's opposite routing for an explicitly supplied
#' batch size that creates multiple loader batches is not the parity target.
#'
#' @section Experimental status:
#' R Torch GMM fitting is experimental. It is suitable for tests, teaching
#' data, and controlled validation, but its runtime on multi-million-row
#' datasets may be substantially longer than Python CellCharter. The graph,
#' aggregation, downstream statistics, and plotting APIs are not affected by
#' this status.
#'
#' @param x observation-by-feature matrix or DelayedMatrix. Rows are read in
#'   chunks and never materialised jointly.
#' @param n_clusters number of mixture components.
#' @param batch_size Maximum observations transferred per E-step, or `"auto"`
#'   for conservative CUDA-memory-aware selection.
#' @param max_iter maximum EM iterations.
#' @param tolerance stop when per-observation NLL change is below this value.
#' @param covariance_regularization diagonal ridge added before Cholesky.
#' @param device Torch device selection.
#' @param seed initialisation seed.
#' @param num_threads Requested CPU thread count, or `NULL` to retain Torch's
#'   current setting.
#' @param initialization_method `"legacy"` for the original base-R streamed
#'   initializer or `"torch"` for the accelerated full-data Torch initializer.
#' @param validate_input Whether to scan every input block for finite values
#'   before fitting. AutoK may perform this scan once for the whole campaign.
#' @param .oom_retries Internal maximum number of complete-fit CUDA OOM retries.
#'   AutoK sets this to zero because it owns campaign-level batch-size retries.
#' @return A `torch_gmm` model.  Only mixture parameters are retained in memory.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' # A tiny, deterministic, CPU-only fit.
#' x <- rbind(
#'   c(-2.0, -2.0), c(-2.2, -1.8), c(-1.8, -2.2),
#'   c(-2.1, -2.3), c(-1.9, -1.7), c(-2.3, -2.1),
#'   c( 2.0,  2.0), c( 2.2,  1.8), c( 1.8,  2.2),
#'   c( 2.1,  2.3), c( 1.9,  1.7), c( 2.3,  2.1)
#' )
#' fit <- torch_gmm_fit(
#'   x, n_clusters = 2, batch_size = 6, max_iter = 3,
#'   device = "cpu", seed = 7, num_threads = 1
#' )
#' fit$n_clusters
#' @export
torch_gmm_fit <- function(x, n_clusters, batch_size = 8192L, max_iter = 100L,
                          tolerance = 1e-3, covariance_regularization = 1e-6,
                          device = "auto", seed = 1L, num_threads = NULL,
                          initialization_method = c("legacy", "torch"),
                          validate_input = TRUE, .oom_retries = 10L) {
  .cc_require_torch()
  reader <- .rc_reader(x)
  k <- .rc_positive_integer(n_clusters, "n_clusters")
  if (k > reader$n) {
    stop("`n_clusters` must be between 1 and nrow(x).", call. = FALSE)
  }
  if (is.character(batch_size) && length(batch_size) == 1L &&
      identical(batch_size, "auto")) {
    tuning <- torch_gmm_autotune_batch_size(
      reader$p, k, reader$n, device = device
    )
    batch_size <- tuning$batch_size
    device <- tuning$device
  } else {
    batch_size <- .rc_positive_integer(batch_size, "batch_size")
    device <- torch_gmm_device(device)$type
  }

  # CUDA memory reported by the driver is necessarily only a snapshot. A
  # concurrent process or allocator fragmentation can invalidate even a
  # conservative estimate. Restart the complete deterministic fit with a
  # smaller transfer batch; never continue from a partially updated model.
  .rc_torch_oom_retry(
    call = function(current_batch) .rc_torch_gmm_fit_once(
      x = x, n_clusters = k, batch_size = current_batch,
      max_iter = max_iter, tolerance = tolerance,
      covariance_regularization = covariance_regularization,
      device = device, seed = seed, num_threads = num_threads,
      initialization_method = initialization_method,
      validate_input = validate_input
    ),
    batch_size = batch_size, device = device, max_retries = .oom_retries
  )
}

.rc_torch_cuda_cleanup <- function() {
  invisible(gc(verbose = FALSE))
  if (!requireNamespace("torch", quietly = TRUE)) return(invisible(NULL))
  ns <- asNamespace("torch")
  if (exists("cuda_empty_cache", envir = ns, inherits = FALSE)) {
    try(get("cuda_empty_cache", envir = ns)(), silent = TRUE)
  }
  invisible(NULL)
}

.rc_torch_is_cuda_oom <- function(condition) {
  grepl(
    "CUDA out of memory|CUDACachingAllocator|cuda.*alloc",
    conditionMessage(condition), ignore.case = TRUE
  )
}

.rc_torch_oom_retry <- function(call, batch_size, device, max_retries = 10L) {
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  max_retries <- .rc_positive_integer(
    max_retries, "max_retries", allow_zero = TRUE
  )
  retries <- 0L
  repeat {
    outcome <- tryCatch(
      list(value = call(batch_size), error = NULL),
      error = function(error) list(value = NULL, error = error)
    )
    if (is.null(outcome$error)) return(outcome$value)
    retry <- identical(device, "cuda") &&
      .rc_torch_is_cuda_oom(outcome$error) &&
      retries < max_retries && batch_size > 1L
    if (!retry) stop(outcome$error)
    next_batch <- as.integer(max(1L, floor(batch_size / 2)))
    retries <- retries + 1L
    message(sprintf(
      paste0(
        "CellCharterR GMM: CUDA OOM at batch_size=%d; discarding the ",
        "partial fit, clearing the cache, and retrying from initialization ",
        "with batch_size=%d (retry %d/%d)."
      ),
      batch_size, next_batch, retries, max_retries
    ))
    batch_size <- next_batch
    .rc_torch_cuda_cleanup()
  }
}

.rc_torch_gmm_fit_once <- function(
    x, n_clusters, batch_size = 8192L, max_iter = 100L,
    tolerance = 1e-3, covariance_regularization = 1e-6,
    device = "auto", seed = 1L, num_threads = NULL,
    initialization_method = c("legacy", "torch"), validate_input = TRUE) {
  reader <- .rc_reader(x)
  k <- .rc_positive_integer(n_clusters, "n_clusters")
  if (is.character(batch_size) && length(batch_size) == 1L &&
      identical(batch_size, "auto")) {
    tuning <- torch_gmm_autotune_batch_size(
      reader$p, k, reader$n, device = device
    )
    batch_size <- tuning$batch_size
    device <- tuning$device
  } else {
    batch_size <- .rc_positive_integer(batch_size, "batch_size")
  }
  max_iter <- .rc_positive_integer(max_iter, "max_iter")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  tolerance <- .rc_nonnegative_number(tolerance, "tolerance", strictly_positive = TRUE)
  covariance_regularization <- .rc_nonnegative_number(
    covariance_regularization, "covariance_regularization", strictly_positive = TRUE
  )
  num_threads <- .rc_positive_integer(num_threads, "num_threads", allow_null = TRUE)
  initialization_method <- match.arg(initialization_method)
  if (!is.logical(validate_input) || length(validate_input) != 1L ||
      is.na(validate_input)) {
    stop("`validate_input` must be TRUE or FALSE.", call. = FALSE)
  }
  # CPU Torch otherwise may inherit a conservative one-thread configuration.
  # Do not combine a large thread count with outer AutoK parallel workers.
  if (!is.null(num_threads)) {
    old_num_threads <- torch::torch_get_num_threads()
    on.exit(torch::torch_set_num_threads(old_num_threads), add = TRUE)
    torch::torch_set_num_threads(num_threads)
  }
  # K = 1 is needed as the lower neighbour when AutoK evaluates a declared
  # range beginning at K = 2, matching CellCharter v0.2.0.
  if (k < 1L || k > reader$n) stop("`n_clusters` must be between 1 and nrow(x).")
  # Validate every disk-backed block before starting an expensive EM fit. This
  # is a streaming scan and never realizes the complete matrix.
  if (validate_input) {
    for (start in .rc_chunk_starts(reader$n, batch_size))
      .rc_get(reader, .rc_chunk_index(start, reader$n, batch_size))
  }

  runtime <- torch_gmm_runtime_info(device, batch_size, seed, num_threads)
  message(sprintf(
    "CellCharterR GMM: device=%s torch=%s cuda=%s gpu=%s batch_size=%d seed=%d cpu_threads=%d",
    runtime$device, runtime$torch_version, runtime$cuda_available,
    ifelse(is.na(runtime$gpu_name), "NA", runtime$gpu_name), runtime$batch_size,
    runtime$seed, runtime$cpu_threads
  ))
  dev <- torch_gmm_device(runtime$device)
  initialization <- if (identical(initialization_method, "torch")) {
    .rc_torch_stream_kmeans(reader, k, batch_size, seed, dev)
  } else {
    .rc_stream_kmeans(reader, k, batch_size, seed)
  }
  eye <- .rc_eye(reader$p, dev)
  collapse_threshold <- 1e-6
  maximum_regularization_attempts <- 8L
  attempted_regularization <- numeric()
  restart_errors <- character()
  fit <- NULL
  for (attempt in seq_len(maximum_regularization_attempts)) {
    effective_regularization <- covariance_regularization * 10^(attempt - 1L)
    if (!is.finite(effective_regularization)) break
    attempted_regularization <- c(attempted_regularization, effective_regularization)
    candidate <- tryCatch({
      initialized <- .rc_hard_initialize(
        reader, batch_size, initialization$means, dev, eye, collapse_threshold
      )
      means <- initialized$means
      weights <- initialized$weights
      covariances <- initialized$covariances
      chol <- .rc_cholesky(covariances, effective_regularization, eye)
      converged <- FALSE
      collapsed_history <- vector("list", max_iter)
      nll_history <- numeric(max_iter)
      best_nll <- Inf
      early_stopping_wait <- 0L
      for (iter in seq_len(max_iter)) {
        # PyCave's first epoch computes weights/means; its second replays the
        # same responsibilities for centered covariances and logs the old-model
        # NLL. These two streamed passes are one complete EM iteration here.
        st <- .rc_estep_stats(reader, batch_size, means, chol, weights, dev)
        nk_cpu <- as.numeric(torch::as_array(st$nk$to(device = "cpu")))
        collapsed <- which(!is.finite(nk_cpu) | nk_cpu < collapse_threshold)
        collapsed_history[[iter]] <- collapsed
        old_means <- means
        old_chol <- chol
        old_weights <- weights
        update <- .rc_mstep_means_weights(
          st$nk, st$sx, nk_cpu, collapsed, collapse_threshold, dev
        )
        means <- update$means
        if (length(collapsed)) {
          for (component in collapsed) {
            row <- ((as.double(seed) + iter + component - 2) %% reader$n) + 1
            means[component, ] <-
              .rc_tensor(.rc_get(reader, as.integer(row)), dev)$squeeze(1)
          }
        }
        covariances <- .rc_centered_covariances(
          reader, batch_size, old_means, old_chol, old_weights, means, st$nk,
          collapsed, dev, eye
        )
        weights <- update$weights
        chol <- .rc_cholesky(covariances, effective_regularization, eye)

        mean_nll <- st$nll / reader$n
        nll_history[[iter]] <- mean_nll
        if (mean_nll < best_nll - tolerance) {
          best_nll <- mean_nll
          early_stopping_wait <- 0L
        } else {
          early_stopping_wait <- early_stopping_wait + 1L
          if (early_stopping_wait >= 2L) {
            converged <- TRUE
            break
          }
        }
      }
      # The monitored NLL precedes each M-step. Persist an explicit evaluation
      # of the final fitted parameters instead of the last callback value.
      final_stats <- .rc_estep_stats(reader, batch_size, means, chol, weights, dev)
      list(
        means = means, weights = weights, covariances = covariances,
        nll = final_stats$nll / reader$n, iterations = iter,
        converged = converged,
        collapsed_components = collapsed_history[seq_len(iter)],
        initialization_collapsed_components = initialized$collapsed,
        initialization_counts = initialized$counts,
        nll_history = nll_history[seq_len(iter)],
        early_stopping = list(
          monitor = "mean_nll_before_m_step", min_delta = tolerance,
          patience = 2L, wait_count = early_stopping_wait,
          best_nll = best_nll
        )
      )
    }, rc_gmm_cholesky_error = function(error) error)
    if (!inherits(candidate, "rc_gmm_cholesky_error")) {
      fit <- candidate
      break
    }
    restart_errors <- c(restart_errors, conditionMessage(candidate))
    if (attempt < maximum_regularization_attempts) {
      message(sprintf(
        "Restarting the complete GMM fit after Cholesky failure; regularization %.17g -> %.17g.",
        effective_regularization, effective_regularization * 10
      ))
    }
  }
  if (is.null(fit)) {
    stop(sprintf(
      "GMM Cholesky failed after %d complete-fit attempts (last regularization %.17g).",
      length(attempted_regularization), utils::tail(attempted_regularization, 1L)
    ), call. = FALSE)
  }
  initialization$centers <- initialization$means
  initialization$means <- NULL
  initialization$branch <-
    "full-data KMeans; R batch_size is internal chunking only"
  # Persist ordinary R arrays, not Torch external pointers. Torch tensors do
  # not survive saveRDS/readRDS across sessions, while AutoK relies on durable
  # per-fit model checkpoints and lazy resumption.
  structure(list(
                 means = torch::as_array(fit$means$to(device = "cpu")),
                 covariances = lapply(fit$covariances, function(z) {
                   torch::as_array(z$to(device = "cpu"))
                 }),
                 weights = as.numeric(torch::as_array(fit$weights$to(device = "cpu"))),
                 n_features = reader$p, n_clusters = k,
                 nll = fit$nll, iterations = fit$iterations,
                 converged = fit$converged,
                 covariance_regularization = utils::tail(attempted_regularization, 1L),
                 requested_covariance_regularization = covariance_regularization,
                 regularization_attempts = attempted_regularization,
                 regularization_max_attempts = maximum_regularization_attempts,
                 regularization_restart_errors = restart_errors,
                 covariance_storage = paste(
                   "Unregularized centered covariances; covariance_regularization",
                   "is added exactly once when constructing the fitted distribution."
                 ),
                 collapsed_components = fit$collapsed_components,
                 collapse_threshold = collapse_threshold,
                 initialization_collapsed_components =
                   fit$initialization_collapsed_components,
                 initialization_counts = fit$initialization_counts,
                 nll_history = fit$nll_history,
                 early_stopping = fit$early_stopping,
                 initialization = initialization,
                 reference_semantics = paste(
                   "CellCharter 0.2.0 CODEX omitted batch_size; PyCave 3.2.1",
                   "one-loader-batch == 1 two-pass path (patience=2)."
                 ),
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
#' @return A numeric observation-by-component posterior-probability matrix.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' model <- structure(list(
#'   means = rbind(c(-1, -1), c(1, 1)),
#'   covariances = list(diag(2), diag(2)),
#'   weights = c(0.5, 0.5), n_features = 2L, n_clusters = 2L,
#'   covariance_regularization = 1e-6
#' ), class = "torch_gmm")
#' probability <- torch_gmm_predict_proba(
#'   model, rbind(c(-1, -1), c(1, 1)), batch_size = 1, device = "cpu"
#' )
#' rowSums(probability)
#' @export
torch_gmm_predict_proba <- function(model, newdata, batch_size = 8192L, device = "auto") {
  .cc_require_torch()
  stopifnot(inherits(model, "torch_gmm")); reader <- .rc_reader(newdata)
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (reader$p != model$n_features) stop("Feature dimensions do not match the fitted model.")
  pars <- .rc_model_device(model, device); out <- vector("list", ceiling(reader$n / batch_size)); b <- 0L
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    b <- b + 1L; xb <- .rc_device_block(reader, idx, pars$device)
    lp <- .rc_logprob(xb, pars$means, pars$chol, pars$weights)
    out[[b]] <- torch::as_array(torch::torch_exp(lp - torch::torch_logsumexp(lp, dim = 2)$unsqueeze(2))$to(device = "cpu"))
  }
  do.call(rbind, out)
}

#' Predict CellCharter domain labels from a Torch GMM
#' @inheritParams torch_gmm_predict_proba
#' @return A zero-based integer label vector matching CellCharter's Python
#'   component convention.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' model <- structure(list(
#'   means = rbind(c(-1, -1), c(1, 1)),
#'   covariances = list(diag(2), diag(2)),
#'   weights = c(0.5, 0.5), n_features = 2L, n_clusters = 2L,
#'   covariance_regularization = 1e-6
#' ), class = "torch_gmm")
#' torch_gmm_predict(
#'   model, rbind(c(-1, -1), c(1, 1)), batch_size = 1, device = "cpu"
#' )
#' @export
torch_gmm_predict <- function(model, newdata, batch_size = 8192L, device = "auto") {
  .cc_require_torch()
  stopifnot(inherits(model, "torch_gmm"))
  reader <- .rc_reader(newdata)
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (reader$p != model$n_features) stop("Feature dimensions do not match the fitted model.")
  pars <- .rc_model_device(model, device)
  labels <- integer(reader$n)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    lp <- .rc_logprob(.rc_device_block(reader, idx, pars$device),
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
#' @return Mean per-observation log likelihood under `model`.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' model <- structure(list(
#'   means = rbind(c(-1, -1), c(1, 1)),
#'   covariances = list(diag(2), diag(2)),
#'   weights = c(0.5, 0.5), n_features = 2L, n_clusters = 2L,
#'   covariance_regularization = 1e-6
#' ), class = "torch_gmm")
#' torch_gmm_log_likelihood(
#'   model, rbind(c(-1, -1), c(1, 1)), batch_size = 1, device = "cpu"
#' )
#' @export
torch_gmm_log_likelihood <- function(model, newdata, batch_size = 8192L, device = "auto") {
  .cc_require_torch()
  stopifnot(inherits(model, "torch_gmm"))
  reader <- .rc_reader(newdata)
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  if (reader$p != model$n_features) stop("Feature dimensions do not match the fitted model.")
  pars <- .rc_model_device(model, device); value <- 0
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    lp <- .rc_logprob(.rc_device_block(reader, idx, pars$device), pars$means, pars$chol, pars$weights)
    value <- value + torch::torch_logsumexp(lp, dim = 2)$sum()$item()
  }
  value / reader$n
}
