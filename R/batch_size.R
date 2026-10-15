.cc_cuda_free_memory <- function() {
  executable <- Sys.which("nvidia-smi")
  if (!nzchar(executable)) return(NA_real_)

  # Query rather than allocate: a failed probe must never risk a CUDA OOM.
  output <- try(system2(
    executable,
    c("--query-gpu=memory.free", "--format=csv,noheader,nounits"),
    stdout = TRUE, stderr = FALSE
  ), silent = TRUE)
  if (inherits(output, "try-error") || !length(output)) return(NA_real_)

  visible <- strsplit(Sys.getenv("CUDA_VISIBLE_DEVICES", ""), ",", fixed = TRUE)[[1L]][1L]
  index <- suppressWarnings(as.integer(visible)) + 1L
  if (!is.finite(index) || index < 1L || index > length(output)) index <- 1L
  mib <- suppressWarnings(as.numeric(trimws(output[[index]])))
  if (length(mib) != 1L || !is.finite(mib) || mib <= 0) return(NA_real_)
  mib * 1024^2
}

.cc_gmm_batch_estimate <- function(n_features, n_clusters, free_memory_bytes,
                                   safety_fraction, limit) {
  # Float32 peak: centered/weighted/solve workspaces plus logits and input.
  bytes_per_observation <- 4 * (6 * n_clusters * n_features +
                                  3 * n_clusters + 2 * n_features)
  fixed_bytes <- 4 * (4 * n_clusters * n_features^2 +
                        8 * n_clusters * n_features)
  usable <- max(0, free_memory_bytes * safety_fraction - fixed_bytes)
  estimated <- floor(usable / bytes_per_observation)
  list(
    batch_size = as.integer(max(1, min(limit, estimated))),
    usable_memory_bytes = as.numeric(usable),
    estimated_bytes_per_observation = as.numeric(bytes_per_observation)
  )
}

#' Select a conservative batch size for the Torch GMM
#'
#' Estimates the largest bounded-memory batch that fits within a conservative
#' fraction of currently free CUDA memory. The function only queries device
#' status; it never allocates a trial tensor and therefore cannot cause an
#' out-of-memory error while tuning. Explicit `free_memory_bytes` is useful for
#' schedulers, reproducible tests, and systems where `nvidia-smi` is unavailable.
#'
#' The estimate covers the largest observation-by-component-by-feature
#' temporaries used by the full-covariance E/M steps and reserves memory for
#' model covariance matrices. It is intentionally conservative, not a promise
#' that every external CUDA workload will remain unchanged.
#'
#' @param n_features Positive number of representation features.
#' @param n_clusters Positive number of mixture components.
#' @param n_obs Optional positive number of observations; limits the result.
#' @param device `"auto"`, `"cuda"`, or `"cpu"`.
#' @param cap Maximum returned batch size.
#' @param safety_fraction Fraction of reported free CUDA memory available to
#'   the estimate. Must be greater than zero and no greater than one.
#' @param free_memory_bytes Optional finite positive byte count overriding the
#'   device query.
#' @param cpu_batch_size Conservative CPU result when the resolved device is
#'   CPU.
#' @return A list containing `batch_size`, resolved `device`, memory estimate,
#'   query source, and tuning parameters.
#' @examplesIf requireNamespace("torch", quietly = TRUE)
#' torch_gmm_autotune_batch_size(30, 12, device = "cpu")$batch_size
#' if (torch::cuda_is_available()) {
#'   torch_gmm_autotune_batch_size(
#'     30, 12, device = "cuda", free_memory_bytes = 8 * 1024^3
#'   )$batch_size
#' }
#' @export
torch_gmm_autotune_batch_size <- function(
    n_features, n_clusters, n_obs = NULL,
    device = c("auto", "cuda", "cpu"), cap = 32768L,
    safety_fraction = 0.35, free_memory_bytes = NULL,
    cpu_batch_size = 8192L) {
  positive_integer <- function(x, name, allow_null = FALSE) {
    if (allow_null && is.null(x)) return(NULL)
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || x != as.integer(x))
      stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
    as.integer(x)
  }
  n_features <- positive_integer(n_features, "n_features")
  n_clusters <- positive_integer(n_clusters, "n_clusters")
  n_obs <- positive_integer(n_obs, "n_obs", allow_null = TRUE)
  cap <- positive_integer(cap, "cap")
  cpu_batch_size <- positive_integer(cpu_batch_size, "cpu_batch_size")
  if (length(safety_fraction) != 1L || !is.finite(safety_fraction) ||
      safety_fraction <= 0 || safety_fraction > 1)
    stop("`safety_fraction` must be one finite number in (0, 1].", call. = FALSE)
  if (!is.null(free_memory_bytes) &&
      (length(free_memory_bytes) != 1L || !is.finite(free_memory_bytes) ||
       free_memory_bytes <= 0))
    stop("`free_memory_bytes` must be NULL or one finite positive number.", call. = FALSE)

  requested_device <- match.arg(device)
  # An explicit CPU estimate is pure arithmetic and must not require Lantern or
  # a working CUDA runtime merely to choose a conservative host batch.
  resolved_device <- if (identical(requested_device, "cpu")) {
    "cpu"
  } else {
    torch_gmm_device(requested_device)$type
  }
  limit <- min(cap, if (is.null(n_obs)) cap else n_obs)

  if (identical(resolved_device, "cpu")) {
    return(list(
      batch_size = as.integer(min(limit, cpu_batch_size)), device = "cpu",
      free_memory_bytes = NA_real_, usable_memory_bytes = NA_real_,
      estimated_bytes_per_observation = NA_real_, memory_source = "cpu_default",
      cap = cap, safety_fraction = safety_fraction
    ))
  }

  source <- if (is.null(free_memory_bytes)) "nvidia-smi" else "argument"
  if (is.null(free_memory_bytes)) free_memory_bytes <- .cc_cuda_free_memory()
  if (!is.finite(free_memory_bytes)) {
    # A deterministic conservative fallback is preferable to guessing total VRAM.
    return(list(
      batch_size = as.integer(min(limit, 8192L)), device = "cuda",
      free_memory_bytes = NA_real_, usable_memory_bytes = NA_real_,
      estimated_bytes_per_observation = NA_real_, memory_source = "fallback",
      cap = cap, safety_fraction = safety_fraction
    ))
  }

  estimate <- .cc_gmm_batch_estimate(
    n_features, n_clusters, free_memory_bytes, safety_fraction, limit
  )

  list(
    batch_size = estimate$batch_size, device = "cuda",
    free_memory_bytes = as.numeric(free_memory_bytes),
    usable_memory_bytes = estimate$usable_memory_bytes,
    estimated_bytes_per_observation = estimate$estimated_bytes_per_observation,
    memory_source = source, cap = cap, safety_fraction = safety_fraction
  )
}
