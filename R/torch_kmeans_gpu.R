# Torch-native, bounded-memory k-means used by the opt-in accelerated GMM
# initialization path. The legacy R initializer remains the default so that
# existing AutoK checkpoints retain their original numerical contract.

.rc_torch_squared_distances <- function(x, centers) {
  distances <- x$pow(2)$sum(dim = 2, keepdim = TRUE) +
    centers$pow(2)$sum(dim = 2)$unsqueeze(1) -
    2 * torch::torch_matmul(x, centers$t())
  torch::torch_clamp(distances, min = 0)
}

.rc_torch_nearest <- function(x, centers) {
  result <- torch::torch_min(.rc_torch_squared_distances(x, centers), dim = 2)
  list(distance = result[[1L]], label = result[[2L]])
}

.rc_torch_kmeans_plus_plus <- function(reader, k, batch_size, seed, device) {
  # Candidate thresholds deliberately come from R's locally scoped RNG. This
  # avoids dependence on the global Torch RNG state and is reproducible for a
  # fixed device/runtime.
  first <- withr::with_seed(seed, sample.int(reader$n, 1L))
  centers <- .rc_tensor(.rc_get(reader, first), device)
  if (k == 1L) return(centers)

  n_local_trials <- 2L + floor(log(k))
  withr::with_seed(seed, {
    # Consume the same first-row draw on every invocation before drawing the
    # weighted candidates.
    sample.int(reader$n, 1L)
    for (component in 2:k) {
      total <- 0
      for (start in .rc_chunk_starts(reader$n, batch_size)) {
        idx <- .rc_chunk_index(start, reader$n, batch_size)
        xb <- .rc_tensor(.rc_get(reader, idx), device)
        total <- total + .rc_torch_nearest(xb, centers)$distance$sum()$item()
      }

      if (!is.finite(total) || total <= 0) {
        row <- as.integer(((as.double(seed) + component - 2) %% reader$n) + 1)
        centers <- torch::torch_cat(
          list(centers, .rc_tensor(.rc_get(reader, row), device)), dim = 1
        )
        next
      }

      thresholds <- sort(stats::runif(n_local_trials, min = 0, max = total))
      candidates <- vector("list", n_local_trials)
      cumulative <- 0
      for (start in .rc_chunk_starts(reader$n, batch_size)) {
        idx <- .rc_chunk_index(start, reader$n, batch_size)
        block <- .rc_get(reader, idx)
        xb <- .rc_tensor(block, device)
        nearest <- as.numeric(torch::as_array(
          .rc_torch_nearest(xb, centers)$distance$to(device = "cpu")
        ))
        within <- cumulative + cumsum(nearest)
        for (candidate in which(vapply(candidates, is.null, logical(1L)))) {
          location <- which(within >= thresholds[[candidate]])[1L]
          if (!is.na(location)) candidates[[candidate]] <- block[location, , drop = FALSE]
        }
        cumulative <- cumulative + sum(nearest)
        if (all(!vapply(candidates, is.null, logical(1L)))) break
      }
      missing <- which(vapply(candidates, is.null, logical(1L)))
      if (length(missing)) {
        fallback <- .rc_get(reader, reader$n)
        for (candidate in missing) candidates[[candidate]] <- fallback
      }
      candidate_tensor <- .rc_tensor(do.call(rbind, candidates), device)

      potentials <- torch::torch_zeros(
        n_local_trials, dtype = torch::torch_float32(), device = device
      )
      for (start in .rc_chunk_starts(reader$n, batch_size)) {
        idx <- .rc_chunk_index(start, reader$n, batch_size)
        xb <- .rc_tensor(.rc_get(reader, idx), device)
        nearest <- .rc_torch_nearest(xb, centers)$distance$unsqueeze(2)
        trial_distances <- .rc_torch_squared_distances(xb, candidate_tensor)
        potentials$add_(torch::torch_minimum(trial_distances, nearest)$sum(dim = 1))
      }
      selected <- torch::torch_argmin(potentials)$item()
      centers <- torch::torch_cat(
        list(centers, candidate_tensor[selected, ]$unsqueeze(1)), dim = 1
      )
    }
  })
  centers
}

.rc_torch_kmeans_tolerance <- function(reader, batch_size, device,
                                        relative_tolerance) {
  total <- torch::torch_zeros(reader$p, dtype = torch::torch_float32(), device = device)
  total_sq <- torch::torch_zeros(reader$p, dtype = torch::torch_float32(), device = device)
  for (start in .rc_chunk_starts(reader$n, batch_size)) {
    idx <- .rc_chunk_index(start, reader$n, batch_size)
    xb <- .rc_tensor(.rc_get(reader, idx), device)
    total$add_(xb$sum(dim = 1))
    total_sq$add_(xb$pow(2)$sum(dim = 1))
  }
  variance <- total_sq / reader$n - (total / reader$n)$pow(2)
  max(0, variance$mean()$item()) * relative_tolerance
}

.rc_torch_stream_kmeans <- function(reader, k, batch_size, seed, device,
                                    max_iter = 300L, tolerance = 1e-4) {
  k <- .rc_positive_integer(k, "k")
  batch_size <- .rc_positive_integer(batch_size, "batch_size")
  seed <- .rc_positive_integer(seed, "seed", allow_zero = TRUE)
  max_iter <- .rc_positive_integer(max_iter, "max_iter")
  tolerance <- .rc_nonnegative_number(tolerance, "tolerance")
  if (k > reader$n) {
    stop("`k` cannot exceed the number of observations.", call. = FALSE)
  }
  if (!inherits(device, "torch_device")) {
    stop("`device` must be a resolved `torch_device`.", call. = FALSE)
  }

  centers <- .rc_torch_kmeans_plus_plus(reader, k, batch_size, seed, device)
  absolute_tolerance <- .rc_torch_kmeans_tolerance(
    reader, batch_size, device, tolerance
  )
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    counts <- torch::torch_zeros(k, dtype = torch::torch_float32(), device = device)
    sums <- torch::torch_zeros(
      c(k, reader$p), dtype = torch::torch_float32(), device = device
    )
    farthest_d2 <- numeric()
    farthest_rows <- matrix(numeric(), nrow = 0L, ncol = reader$p)

    for (start in .rc_chunk_starts(reader$n, batch_size)) {
      idx <- .rc_chunk_index(start, reader$n, batch_size)
      block <- .rc_get(reader, idx)
      xb <- .rc_tensor(block, device)
      nearest <- .rc_torch_nearest(xb, centers)
      membership <- torch::nnf_one_hot(nearest$label, num_classes = k)$to(
        dtype = torch::torch_float32()
      )
      counts$add_(membership$sum(dim = 1))
      sums$add_(torch::torch_matmul(membership$t(), xb))

      # Retain only K rows for deterministic empty-cluster recovery.
      nearest_cpu <- as.numeric(torch::as_array(
        nearest$distance$to(device = "cpu")
      ))
      keep <- order(nearest_cpu, decreasing = TRUE)[seq_len(min(k, length(nearest_cpu)))]
      farthest_d2 <- c(farthest_d2, nearest_cpu[keep])
      farthest_rows <- rbind(farthest_rows, block[keep, , drop = FALSE])
      retain <- order(farthest_d2, decreasing = TRUE)[seq_len(min(k, length(farthest_d2)))]
      farthest_d2 <- farthest_d2[retain]
      farthest_rows <- farthest_rows[retain, , drop = FALSE]
    }

    counts_cpu <- as.numeric(torch::as_array(counts$to(device = "cpu")))
    updated <- centers$clone()
    nonempty <- which(counts_cpu > 0)
    if (length(nonempty)) updated[nonempty, ] <- sums[nonempty, ] / counts[nonempty]$unsqueeze(2)
    empty <- which(counts_cpu == 0)
    if (length(empty)) {
      replacements <- farthest_rows[seq_len(length(empty)), , drop = FALSE]
      updated[empty, ] <- .rc_tensor(replacements, device)
    }
    shift <- (updated - centers)$pow(2)$sum()$item()
    centers <- updated
    if (shift <= absolute_tolerance) {
      converged <- TRUE
      break
    }
  }

  list(
    means = torch::as_array(centers$to(device = "cpu")),
    iterations = iter,
    converged = converged,
    method = "Torch streamed full-data greedy k-means++/Lloyd",
    reference_method = "sklearn.cluster.KMeans(n_init=1)",
    rng_note = paste(
      "Candidate draws use a locally scoped R seed; fixed-device repeats are",
      "deterministic subject to the selected Torch backend's deterministic guarantees."
    )
  )
}
