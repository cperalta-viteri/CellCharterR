#' Fowlkes--Mallows agreement between two clusterings
#'
#' This is the pair-counting definition used by
#' `sklearn.metrics.fowlkes_mallows_score`, the default similarity function in
#' CellCharter 0.2.0. Cluster names do not need to agree between the inputs.
#'
#' @param labels_a,labels_b Cluster labels of equal length.
#' @return A scalar in `[0, 1]`.
#' @examples
#' labels_a <- c(1L, 1L, 2L, 2L)
#' labels_b <- c("a", "a", "b", "b")
#' cellcharter_fowlkes_mallows(labels_a, labels_b)
#' @export
cellcharter_fowlkes_mallows <- function(labels_a, labels_b) {
  if (length(labels_a) != length(labels_b)) {
    stop("`labels_a` and `labels_b` must have the same length.")
  }
  if (!length(labels_a) || anyNA(labels_a) || anyNA(labels_b)) {
    stop("Cluster labels must be non-empty and cannot contain missing values.")
  }

  contingency <- table(labels_a, labels_b)
  choose2 <- function(x) x * (x - 1) / 2
  true_positive <- sum(choose2(contingency))
  predicted_positive <- sum(choose2(rowSums(contingency)))
  condition_positive <- sum(choose2(colSums(contingency)))

  denominator <- sqrt(predicted_positive * condition_positive)
  if (denominator == 0) 0 else as.numeric(true_positive / denominator)
}

.autok_validate_integer <- function(x, name, lower = 1L) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x != as.integer(x) || x < lower) {
    stop(sprintf("`%s` must be one integer >= %d.", name, lower))
  }
  as.integer(x)
}

.autok_k_grid <- function(n_clusters) {
  # A length-two numeric vector is the R equivalent of CellCharter's Python
  # tuple: it declares selectable bounds and is expanded by one K at each end.
  if (is.numeric(n_clusters) && !is.list(n_clusters) && length(n_clusters) == 2L) {
    bounds <- as.integer(n_clusters)
    if (anyNA(bounds) || any(bounds != n_clusters) || bounds[[1L]] > bounds[[2L]] || bounds[[1L]] < 1L) {
      stop("A two-element `n_clusters` must contain increasing positive integer bounds.")
    }
    fit_k <- seq.int(max(1L, bounds[[1L]] - 1L), bounds[[2L]] + 1L)
    requested_k <- seq.int(bounds[[1L]], bounds[[2L]])
  } else {
    # A list mirrors CellCharter's explicit-list form. It must already contain
    # both neighboring boundary values because only interior K are selectable.
    fit_k <- as.integer(unlist(n_clusters, use.names = FALSE))
    if (length(fit_k) < 3L || anyNA(fit_k) || any(fit_k < 1L) ||
        anyDuplicated(fit_k) || !identical(fit_k, seq.int(min(fit_k), max(fit_k)))) {
      stop("An explicit `n_clusters` list must contain at least three consecutive positive integers.")
    }
    requested_k <- fit_k[-c(1L, length(fit_k))]
  }
  list(fit_k = fit_k, requested_k = requested_k)
}

.autok_checkpoint_path <- function(checkpoint_dir, k, run) {
  file.path(checkpoint_dir, "fits", sprintf("k_%03d", k), sprintf("run_%03d.rds", run))
}

.autok_atomic_save_rds <- function(object, path) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "-"), tmpdir = directory, fileext = ".tmp"
  )
  on.exit(unlink(temporary), add = TRUE)
  base::saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically move checkpoint into place: ", path)
  }
  invisible(path)
}

.autok_read_checkpoint <- function(path, expected, n_obs, compatible_data_ids = character()) {
  artifact <- tryCatch(
    readRDS(path),
    error = function(e) stop("Cannot read AutoK checkpoint `", path, "`: ", conditionMessage(e))
  )
  required <- c(
    "format_version", "k", "run", "seed", "n_obs", "n_features",
    "data_id", "nll", "model", "labels"
  )
  missing_fields <- setdiff(required, names(artifact))
  if (length(missing_fields)) {
    stop("AutoK checkpoint `", path, "` is missing: ", paste(missing_fields, collapse = ", "))
  }
  for (field in c("k", "run", "seed", "n_obs", "n_features")) {
    if (!identical(artifact[[field]], expected[[field]])) {
      stop("AutoK checkpoint `", path, "` does not match current `", field, "`.")
    }
  }
  accepted_ids <- unique(c(expected$data_id, compatible_data_ids))
  if (!artifact$data_id %in% accepted_ids) {
    stop("AutoK checkpoint `", path, "` does not match current `data_id`.")
  }
  if (!artifact$format_version %in% c(1L, 2L) ||
      length(artifact$labels) != n_obs || anyNA(artifact$labels)) {
    stop("AutoK checkpoint `", path, "` has invalid version or labels.")
  }
  if (length(artifact$nll) != 1L || !is.numeric(artifact$nll) || !is.finite(artifact$nll)) {
    stop("AutoK checkpoint `", path, "` has an invalid negative log-likelihood.")
  }
  artifact
}

.autok_runtime_args <- function(args) {
  args[setdiff(names(args), c("batch_size", "device", "num_threads", "validate_input"))]
}

.autok_resume_contract <- function(input_data_id, fit_args, predict_args,
                                    default_torch, implementation_digest) {
  list(
    schema = "CellCharterR_autok_resume_identity",
    schema_version = 2L,
    input_data_id = input_data_id,
    algorithm = if (default_torch) {
      # Increment only for a change that makes fitted models scientifically
      # incompatible. Package releases and performance-only changes do not.
      "torch-full-covariance-em-v1"
    } else {
      implementation_digest
    },
    fit_args = if (default_torch) .autok_runtime_args(fit_args) else fit_args,
    predict_args = if (default_torch) .autok_runtime_args(predict_args) else predict_args
  )
}

.autok_legacy_data_id <- function(checkpoint_dir, current) {
  path <- file.path(checkpoint_dir, "autok_state.rds")
  if (!file.exists(path)) {
    # Version-1 fit files did not store their input checksum separately and
    # state was written only after a complete K round. Permit a one-time,
    # clearly warned migration of an interrupted partial round after checking
    # that every artifact has one consistent legacy identity. K/run/seed and
    # dimensions are still checked individually when each file is opened.
    files <- list.files(file.path(checkpoint_dir, "fits"), pattern = "^run_[0-9]+[.]rds$",
                        recursive = TRUE, full.names = TRUE)
    ids <- unique(unlist(lapply(files, function(file) {
      value <- tryCatch(readRDS(file), error = function(e) NULL)
      if (is.null(value) || !identical(value$format_version, 1L)) character()
      else as.character(value$data_id)
    }), use.names = FALSE))
    if (length(ids) == 1L) {
      warning(
        "Migrating a partial version-1 AutoK run without `autok_state.rds`. ",
        "The legacy format cannot re-verify the input checksum; only continue ",
        "if this checkpoint directory belongs to the same representation and parameters.",
        call. = FALSE
      )
      return(ids)
    }
    return(character())
  }
  old <- tryCatch(readRDS(path), error = function(e) NULL)
  required <- c("data_id", "input_data_id", "fit_args", "predict_args", "fit_k",
                "candidate_k", "seed", "n_obs", "n_features")
  if (is.null(old) || !all(required %in% names(old))) return(character())
  equal <- identical(old$input_data_id, current$input_data_id) &&
    identical(as.integer(old$fit_k), as.integer(current$fit_k)) &&
    identical(as.integer(old$candidate_k), as.integer(current$candidate_k)) &&
    identical(as.integer(old$seed), as.integer(current$seed)) &&
    identical(as.integer(old$n_obs), as.integer(current$n_obs)) &&
    identical(as.integer(old$n_features), as.integer(current$n_features)) &&
    identical(.cc_contract_canonical(.autok_runtime_args(old$fit_args)),
              .cc_contract_canonical(.autok_runtime_args(current$fit_args))) &&
    identical(.cc_contract_canonical(.autok_runtime_args(old$predict_args)),
              .cc_contract_canonical(.autok_runtime_args(current$predict_args)))
  if (isTRUE(equal)) unique(c(as.character(old$data_id), old$compatible_data_ids))
  else character()
}

.autok_cuda_cleanup <- function() {
  invisible(gc(verbose = FALSE))
  if (!requireNamespace("torch", quietly = TRUE)) return(invisible(NULL))
  ns <- asNamespace("torch")
  if (exists("cuda_empty_cache", envir = ns, inherits = FALSE)) {
    try(get("cuda_empty_cache", envir = ns)(), silent = TRUE)
  }
  invisible(NULL)
}

.autok_is_cuda_oom <- function(condition) {
  text <- conditionMessage(condition)
  grepl("CUDA out of memory|CUDACachingAllocator|cuda.*alloc", text,
        ignore.case = TRUE)
}

.autok_oom_retry <- function(call, args, stage, minimum_batch_size = 256L) {
  repeat {
    outcome <- tryCatch(list(value = call(args), error = NULL),
                        error = function(e) list(value = NULL, error = e))
    if (is.null(outcome$error)) return(list(value = outcome$value, args = args))
    batch <- args$batch_size
    can_retry <- .autok_is_cuda_oom(outcome$error) && is.numeric(batch) &&
      length(batch) == 1L && is.finite(batch) && batch > minimum_batch_size
    if (!can_retry) stop(outcome$error)
    next_batch <- as.integer(max(minimum_batch_size, floor(batch / 2)))
    message(sprintf(
      "CellCharterR AutoK: CUDA OOM during %s at batch_size=%d; clearing cache and retrying with %d.",
      stage, as.integer(batch), next_batch
    ))
    args$batch_size <- next_batch
    .autok_cuda_cleanup()
  }
}

.autok_format_duration <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("unknown")
  seconds <- as.integer(round(seconds))
  days <- seconds %/% 86400L
  hours <- (seconds %% 86400L) %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  secs <- seconds %% 60L
  if (days > 0L) sprintf("%dd %02dh %02dm", days, hours, minutes)
  else if (hours > 0L) sprintf("%dh %02dm", hours, minutes)
  else if (minutes > 0L) sprintf("%dm %02ds", minutes, secs)
  else sprintf("%ds", secs)
}

.autok_read_labels <- function(path, n_obs) {
  artifact <- tryCatch(
    readRDS(path),
    error = function(e) stop(
      "Cannot read AutoK checkpoint `", path, "`: ", conditionMessage(e)
    )
  )
  value <- artifact$labels
  if (length(value) != n_obs || anyNA(value)) {
    stop("AutoK checkpoint `", path, "` has invalid labels.")
  }
  value
}

.autok_mirror_stability <- function(edge_stability, fit_k) {
  candidate_k <- fit_k[-c(1L, length(fit_k))]
  if (!nrow(edge_stability)) {
    return(matrix(numeric(), nrow = length(candidate_k), ncol = 0L,
                  dimnames = list(as.character(candidate_k), NULL)))
  }

  mirrored <- lapply(candidate_k, function(k) {
    # CellCharter 0.2.0 concatenates the K -> K+1 edge before K-1 -> K.
    right <- edge_stability$score[edge_stability$lower_k == k]
    left <- edge_stability$score[edge_stability$upper_k == k]
    c(right, left)
  })
  lengths <- vapply(mirrored, length, integer(1))
  if (length(unique(lengths)) != 1L) {
    stop("Internal error: neighboring-K stability edges are incomplete.")
  }
  result <- do.call(rbind, mirrored)
  rownames(result) <- as.character(candidate_k)
  result
}

.autok_mape <- function(previous, current) {
  if (length(previous) != length(current) || any(!is.finite(previous)) || any(!is.finite(current))) {
    stop("Cannot compute AutoK convergence from invalid stability means.")
  }
  # Match sklearn's float64 MAPE convention, including its finite treatment of
  # a zero previous value.
  mean(abs(current - previous) / pmax(abs(previous), .Machine$double.eps))
}

.autok_default_fit <- function(x, n_clusters, seed, ...) {
  torch_gmm_fit(x, n_clusters = n_clusters, seed = seed, ...)
}

.autok_default_predict <- function(model, newdata, ...) {
  torch_gmm_predict(model, newdata, ...)
}

.autok_default_nll <- function(model) {
  value <- model$nll
  if (is.null(value)) value <- model$nll_
  value
}

.autok_validate_hook_args <- function(args, reserved, name) {
  if (!is.list(args)) stop("`", name, "` must be a list.")
  duplicates <- intersect(names(args), reserved)
  if (length(duplicates)) {
    stop("`", name, "` cannot override: ", paste(duplicates, collapse = ", "))
  }
  args
}

#' Select a stable number of CellCharter domains
#'
#' Reproduces the neighboring-K stability algorithm in CellCharter 0.2.0. For
#' every run, labels at K are compared with all earlier-run labels at K+1. The
#' two edges incident to each interior K are then mirrored and combined. Fits
#' are checkpointed individually before stability is updated, so an interrupted
#' analysis resumes without repeating completed fits.
#'
#' @param x Observation-by-feature matrix accepted by `fit_fun`.
#' @param n_clusters Two integer bounds for selectable K. The fitted grid is
#'   expanded by one at each boundary, as in a CellCharter Python tuple. An R
#'   list instead denotes an already-expanded consecutive grid.
#' @param checkpoint_dir Directory for atomic per-fit and state checkpoints.
#' @param max_runs Maximum repeated fits for each K.
#' @param convergence_tol Stop when the MAPE between consecutive stability
#'   curves is below this non-negative threshold. The first possible stop is
#'   after run three, matching CellCharter 0.2.0.
#' @param seed Base seed. Run `i` uses `seed + i - 1` for every K.
#' @param fit_fun Function called as `fit_fun(x, n_clusters, seed, ...)`.
#' @param predict_fun Function called as `predict_fun(model, newdata, ...)`.
#' @param nll_fun Function returning the model's scalar negative log-likelihood.
#' @param similarity_fun Label-agreement function. Defaults to the
#'   Fowlkes--Mallows score.
#' @param fit_args,predict_args Additional named arguments for the hooks.
#' @param data_id Optional immutable representation checksum or provenance ID.
#'   It is combined with scientifically material fit/predict arguments and an
#'   explicit algorithm-compatibility version. Runtime-only controls (device,
#'   batch size, and thread count) and package upgrades do not invalidate a
#'   compatible completed fit; the exact implementation digest is retained as
#'   provenance.
#' @param resume Reuse matching completed fit checkpoints.
#' @return A `cellcharter_autok` object containing stability, convergence
#'   history, and paths to the lowest-NLL model checkpoint for every K. Label
#'   vectors remain in the per-fit RDS files rather than being duplicated in RAM.
#' @examples
#' # Small deterministic hooks keep this example fast and independent of a GPU.
#' toy <- rbind(c(-2, -2), c(-2, -1), c(-1, -2),
#'              c(2, 2), c(2, 1), c(1, 2))
#' toy_fit <- function(x, n_clusters, seed, device) {
#'   fit <- withr::with_seed(
#'     seed, stats::kmeans(x, centers = n_clusters, nstart = 1)
#'   )
#'   list(labels = fit$cluster, nll = fit$tot.withinss, device = device)
#' }
#' toy_predict <- function(model, newdata) model$labels
#' checkpoint_dir <- tempfile("cellcharter-autok-")
#' autok <- cellcharter_autok_fit(
#'   toy, n_clusters = c(2L, 2L), checkpoint_dir = checkpoint_dir,
#'   max_runs = 2L, seed = 7L, fit_fun = toy_fit,
#'   predict_fun = toy_predict, nll_fun = function(model) model$nll,
#'   fit_args = list(device = "cpu")
#' )
#' autok$best_k
#' unlink(checkpoint_dir, recursive = TRUE)
#' @export
cellcharter_autok_fit <- function(
    x, n_clusters, checkpoint_dir, max_runs = 10L, convergence_tol = 0.01,
    seed = 1L, fit_fun = NULL, predict_fun = NULL, nll_fun = NULL,
    similarity_fun = cellcharter_fowlkes_mallows,
    fit_args = list(), predict_args = list(), data_id = NULL, resume = TRUE) {
  grid <- .autok_k_grid(n_clusters)
  fit_k <- grid$fit_k
  candidate_k <- fit_k[-c(1L, length(fit_k))]
  max_runs <- .autok_validate_integer(max_runs, "max_runs")
  seed <- .autok_validate_integer(seed, "seed", lower = 0L)
  if (length(convergence_tol) != 1L || !is.finite(convergence_tol) || convergence_tol < 0) {
    stop("`convergence_tol` must be one finite non-negative number.")
  }
  if (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L || !nzchar(checkpoint_dir)) {
    stop("`checkpoint_dir` must be one non-empty path.")
  }
  if (!is.null(data_id) && (length(data_id) != 1L || is.na(data_id))) {
    stop("`data_id` must be NULL or one non-missing value.")
  }
  input_data_id <- if (is.null(data_id)) NA_character_ else as.character(data_id)
  if (!is.function(similarity_fun)) stop("`similarity_fun` must be a function.")
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("`resume` must be TRUE or FALSE.")
  }
  fit_args <- .autok_validate_hook_args(fit_args, c("x", "n_clusters", "seed"), "fit_args")
  predict_args <- .autok_validate_hook_args(predict_args, c("model", "newdata"), "predict_args")
  fit_fun <- fit_fun %||% .autok_default_fit
  predict_fun <- predict_fun %||% .autok_default_predict
  nll_fun <- nll_fun %||% .autok_default_nll
  if (!is.function(fit_fun) || !is.function(predict_fun) || !is.function(nll_fun)) {
    stop("`fit_fun`, `predict_fun`, and `nll_fun` must be functions.")
  }
  # The default Torch fitter otherwise performs the same complete finite-value
  # scan for every K/run. Validate once for this immutable AutoK invocation and
  # pass an explicit trusted flag to each fit. Custom fit hooks retain their own
  # validation policy.
  if (identical(fit_fun, .autok_default_fit)) {
    validation_batch_size <- fit_args$batch_size %||% 8192L
    if (is.character(validation_batch_size) &&
        length(validation_batch_size) == 1L &&
        identical(validation_batch_size, "auto")) {
      tuning <- torch_gmm_autotune_batch_size(
        n_features = ncol(x), n_clusters = max(fit_k), n_obs = nrow(x),
        device = fit_args$device %||% "auto"
      )
      validation_batch_size <- tuning$batch_size
      fit_args$batch_size <- tuning$batch_size
      fit_args$device <- tuning$device
    }
    validation_batch_size <- .rc_positive_integer(
      validation_batch_size, "fit_args$batch_size"
    )
    validation_reader <- .rc_reader(x)
    for (start in .rc_chunk_starts(validation_reader$n, validation_batch_size)) {
      .rc_get(
        validation_reader,
        .rc_chunk_index(start, validation_reader$n, validation_batch_size)
      )
    }
    fit_args$validate_input <- FALSE
  }
  implementation_functions <- list(
    fit = fit_fun, predict = predict_fun, nll = nll_fun
  )
  # The public defaults are deliberately small wrappers.  Hash their
  # transitive Torch implementation as well, otherwise a numerical change in
  # an internal EM helper could incorrectly reuse old per-fit checkpoints even
  # though the wrapper body itself had not changed.
  if (identical(fit_fun, .autok_default_fit) ||
      identical(predict_fun, .autok_default_predict)) {
    implementation_functions <- c(
      implementation_functions,
      list(
        torch_gmm_fit = torch_gmm_fit,
        torch_gmm_predict = torch_gmm_predict,
        torch_gmm_predict_proba = torch_gmm_predict_proba,
        reader = .rc_reader,
        tensor = .rc_tensor,
        kmeans_distances = .rc_kmeans_distances,
        kmeans_tolerance = .rc_kmeans_tolerance,
        kmeans_plus_plus = .rc_kmeans_plus_plus,
        stream_kmeans = .rc_stream_kmeans,
        torch_squared_distances = .rc_torch_squared_distances,
        torch_nearest = .rc_torch_nearest,
        torch_kmeans_plus_plus = .rc_torch_kmeans_plus_plus,
        torch_kmeans_tolerance = .rc_torch_kmeans_tolerance,
        torch_stream_kmeans = .rc_torch_stream_kmeans,
        log_probability = .rc_logprob,
        cholesky = .rc_cholesky,
        estep = .rc_estep_stats,
        mstep_means_weights = .rc_mstep_means_weights,
        hard_initialization = .rc_hard_initialize,
        centered_covariances = .rc_centered_covariances,
        model_device = .rc_model_device
      )
    )
  }
  implementation_digest <- cellcharter_implementation_digest(
    functions = implementation_functions
  )
  default_torch <- identical(fit_fun, .autok_default_fit) &&
    identical(predict_fun, .autok_default_predict) && identical(nll_fun, .autok_default_nll)
  resume_contract <- .autok_resume_contract(
    input_data_id, fit_args, predict_args, default_torch, implementation_digest
  )
  data_id <- digest::digest(
    .cc_contract_canonical(resume_contract), algo = "sha256", serialize = TRUE
  )

  n_obs <- nrow(x)
  n_features <- ncol(x)
  if (length(n_obs) != 1L || !is.finite(n_obs) || n_obs < 1L ||
      length(n_features) != 1L || !is.finite(n_features) || n_features < 1L) {
    stop("`x` must have at least one row and one column.")
  }
  checkpoint_dir <- normalizePath(checkpoint_dir, winslash = "/", mustWork = FALSE)
  dir.create(file.path(checkpoint_dir, "fits"), recursive = TRUE, showWarnings = FALSE)
  legacy_state_exists <- file.exists(file.path(checkpoint_dir, "autok_state.rds"))
  legacy_data_ids <- if (resume && default_torch) .autok_legacy_data_id(
    checkpoint_dir,
    list(input_data_id = input_data_id, fit_args = fit_args, predict_args = predict_args,
         fit_k = fit_k, candidate_k = candidate_k, seed = seed,
         n_obs = n_obs, n_features = n_features)
  ) else character()
  if (length(legacy_data_ids)) {
    message(if (legacy_state_exists) {
      "Validated legacy AutoK state; completed fits will be resumed across the package upgrade."
    } else {
      "Recognized a partial legacy AutoK run; completed fits will be resumed across the package upgrade."
    })
  }

  edge_stability <- data.frame(
    run = integer(), previous_run = integer(), lower_k = integer(),
    upper_k = integer(), score = numeric()
  )
  fit_index <- data.frame(
    k = integer(), run = integer(), seed = integer(), nll = numeric(),
    path = character(), resumed = logical(), stringsAsFactors = FALSE
  )
  convergence_history <- data.frame(run = integer(), stability_change = numeric())
  previous_stability <- NULL
  converged <- FALSE
  runs_completed <- 0L
  campaign_total <- length(fit_k) * max_runs
  campaign_completed <- 0L
  new_fit_seconds <- numeric()

  for (run in seq_len(max_runs)) {
    run_seed <- as.integer(seed + run - 1L)
    if (is.na(run_seed)) stop("Run seed exceeds R's integer range.")
    current_round_labels <- if (run > 1L) {
      stats::setNames(vector("list", length(fit_k) - 1L), fit_k[-length(fit_k)])
    } else {
      NULL
    }

    for (k in fit_k) {
      path <- .autok_checkpoint_path(checkpoint_dir, k, run)
      expected <- list(
        k = as.integer(k), run = as.integer(run), seed = run_seed,
        n_obs = as.integer(n_obs), n_features = as.integer(n_features), data_id = data_id
      )
      reused <- resume && file.exists(path)
      if (reused) {
        artifact <- .autok_read_checkpoint(
          path, expected, n_obs, compatible_data_ids = legacy_data_ids
        )
      } else {
        fit_started <- proc.time()[["elapsed"]]
        fit_attempt <- .autok_oom_retry(
          function(args) withr::with_seed(
            run_seed,
            do.call(
              fit_fun,
              c(list(x = x, n_clusters = as.integer(k), seed = run_seed), args)
            )
          ),
          fit_args, stage = sprintf("K=%d run=%d fitting", k, run)
        )
        model <- fit_attempt$value
        # Retain a successful reduced batch for later fits in this campaign so
        # the same known-too-large allocation is not retried for every K.
        fit_args$batch_size <- fit_attempt$args$batch_size %||% fit_args$batch_size
        predict_attempt <- .autok_oom_retry(
          function(args) do.call(
            predict_fun,
            c(list(model = model, newdata = x), args)
          ),
          predict_args, stage = sprintf("K=%d run=%d prediction", k, run)
        )
        labels <- predict_attempt$value
        predict_args$batch_size <- predict_attempt$args$batch_size %||% predict_args$batch_size
        if (length(labels) != n_obs || anyNA(labels)) {
          stop("`predict_fun` must return one non-missing label per observation.")
        }
        nll <- nll_fun(model)
        if (length(nll) != 1L || !is.numeric(nll) || !is.finite(nll)) {
          stop("`nll_fun` must return one finite numeric negative log-likelihood.")
        }
        artifact <- list(
          format_version = 2L, k = as.integer(k), run = as.integer(run),
          seed = run_seed, n_obs = as.integer(n_obs), n_features = as.integer(n_features),
          data_id = data_id, nll = as.numeric(nll), model = model, labels = labels,
          effective_fit_batch_size = fit_attempt$args$batch_size %||% NA_integer_,
          effective_predict_batch_size = predict_attempt$args$batch_size %||% NA_integer_
        )
        .autok_atomic_save_rds(artifact, path)
        fit_elapsed <- proc.time()[["elapsed"]] - fit_started
        new_fit_seconds <- c(new_fit_seconds, fit_elapsed)
      }
      campaign_completed <- campaign_completed + 1L
      fit_index <- rbind(
        fit_index,
        data.frame(k = k, run = run, seed = run_seed, nll = artifact$nll,
                   path = path, resumed = reused, stringsAsFactors = FALSE)
      )
      if (run > 1L && k != fit_k[[length(fit_k)]]) {
        current_round_labels[[as.character(k)]] <- artifact$labels
      }
      rm(artifact)
      if (!reused) {
        .autok_cuda_cleanup()
        mean_fit_seconds <- mean(new_fit_seconds)
        eta_seconds <- mean_fit_seconds * (campaign_total - campaign_completed)
        message(sprintf(
          "CellCharterR AutoK: %d/%d fits complete; last=%s; ETA through max_runs=%s.",
          campaign_completed, campaign_total, .autok_format_duration(fit_elapsed),
          .autok_format_duration(eta_seconds)
        ))
      }
    }

    if (run > 1L) {
      scores <- vector("list", (run - 1L) * (length(fit_k) - 1L))
      score_index <- 0L
      # The current round's lower-K labels are reused for every prior run.
      # Read each vector once, then release the cache after this round. This
      # preserves CellCharter's comparison order while avoiding repeated RDS
      # decompression and bounds memory to one round of label vectors.
      current_labels <- current_round_labels
      # Preserve v0.2.0's comparison order: prior run outside, K edge inside.
      for (previous_run in seq_len(run - 1L)) {
        for (edge in seq_len(length(fit_k) - 1L)) {
          lower_k <- fit_k[[edge]]
          upper_k <- fit_k[[edge + 1L]]
          previous_path <- .autok_checkpoint_path(checkpoint_dir, upper_k, previous_run)
          previous_labels <- .autok_read_labels(previous_path, n_obs)
          score <- similarity_fun(current_labels[[as.character(lower_k)]], previous_labels)
          if (length(score) != 1L || !is.numeric(score) || !is.finite(score)) {
            stop("`similarity_fun` must return one finite numeric value.")
          }
          score_index <- score_index + 1L
          scores[[score_index]] <- data.frame(
            run = run, previous_run = previous_run, lower_k = lower_k,
            upper_k = upper_k, score = as.numeric(score)
          )
          rm(previous_labels)
        }
      }
      rm(current_labels)
      edge_stability <- rbind(edge_stability, do.call(rbind, scores))
      stability <- .autok_mirror_stability(edge_stability, fit_k)
      stability_mean <- rowMeans(stability)

      if (!is.null(previous_stability)) {
        stability_change <- .autok_mape(rowMeans(previous_stability), stability_mean)
        convergence_history <- rbind(
          convergence_history,
          data.frame(run = run, stability_change = stability_change)
        )
        if (stability_change < convergence_tol) converged <- TRUE
      }
      previous_stability <- stability
    }
    rm(current_round_labels)

    runs_completed <- run
    partial <- list(
      format_version = 2L, fit_k = fit_k, candidate_k = candidate_k,
      max_runs = max_runs, convergence_tol = convergence_tol, seed = seed,
      fit_args = fit_args, predict_args = predict_args,
      resume_contract = resume_contract,
      compatible_data_ids = legacy_data_ids,
      input_data_id = input_data_id, implementation_digest = implementation_digest,
      data_id = data_id, n_obs = n_obs, n_features = n_features,
      runs_completed = runs_completed, converged = converged,
      edge_stability = edge_stability, convergence_history = convergence_history
    )
    .autok_atomic_save_rds(partial, file.path(checkpoint_dir, "autok_state.rds"))
    if (converged) break
  }

  stability <- if (runs_completed > 1L) {
    .autok_mirror_stability(edge_stability, fit_k)
  } else {
    NULL
  }
  stability_mean <- if (is.null(stability)) NULL else rowMeans(stability)
  best_k <- if (is.null(stability_mean)) NA_integer_ else candidate_k[[which.max(stability_mean)]]

  best_rows <- lapply(fit_k, function(k) {
    entries <- fit_index[fit_index$k == k, , drop = FALSE]
    entries[which.min(entries$nll), , drop = FALSE]
  })
  best_index <- do.call(rbind, best_rows)
  rownames(best_index) <- as.character(best_index$k)

  result <- structure(
    list(
      reference_version = "CellCharter 0.2.0",
      reference_commit = "2a8fdde8213e32e5a00f1f9571bffd3dbb825a79",
      fit_k = fit_k, candidate_k = candidate_k,
      requested_k = grid$requested_k,
      max_runs = max_runs, convergence_tol = convergence_tol, seed = seed,
      fit_args = fit_args, predict_args = predict_args,
      resume_contract = resume_contract,
      compatible_data_ids = legacy_data_ids,
      input_data_id = input_data_id, implementation_digest = implementation_digest,
      data_id = data_id, n_obs = n_obs, n_features = n_features,
      runs_completed = runs_completed, converged = converged,
      stability = stability, stability_mean = stability_mean,
      edge_stability = edge_stability,
      convergence_history = convergence_history,
      best_k = as.integer(best_k), best_index = best_index,
      fit_index = fit_index, checkpoint_dir = checkpoint_dir
    ),
    class = "cellcharter_autok"
  )
  .autok_atomic_save_rds(result, file.path(checkpoint_dir, "autok_result.rds"))
  result
}

# Base R has no null-coalescing operator; keeping this helper local avoids an
# additional package dependency in the checkpointing path.
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Return the most stable K from an AutoK fit
#'
#' @param object A fitted `cellcharter_autok` object.
#' @return The first K attaining the maximum mean mirrored stability.
#' @examples
#' autok <- structure(
#'   list(candidate_k = 2:4, stability_mean = c(0.7, 0.9, 0.8)),
#'   class = "cellcharter_autok"
#' )
#' cellcharter_autok_best_k(autok)
#' @export
cellcharter_autok_best_k <- function(object) {
  if (!inherits(object, "cellcharter_autok")) stop("`object` must be a `cellcharter_autok` fit.")
  if (is.null(object$stability_mean)) {
    stop("Stability and best K require `max_runs` greater than one.")
  }
  as.integer(object$candidate_k[[which.max(object$stability_mean)]])
}

#' Load the lowest-NLL model fitted for K
#'
#' @param object A fitted `cellcharter_autok` object.
#' @param k Cluster count. Defaults to the most stable K.
#' @return The checkpointed fitted model.
#' @examples
#' checkpoint <- tempfile(fileext = ".rds")
#' saveRDS(list(model = list(n_clusters = 2L, device = "cpu")), checkpoint)
#' autok <- structure(
#'   list(
#'     candidate_k = 2L, stability_mean = 1, best_k = 2L,
#'     best_index = data.frame(k = 2L, path = checkpoint)
#'   ),
#'   class = "cellcharter_autok"
#' )
#' cellcharter_autok_model(autok)
#' unlink(checkpoint)
#' @export
cellcharter_autok_model <- function(object, k = NULL) {
  if (!inherits(object, "cellcharter_autok")) stop("`object` must be a `cellcharter_autok` fit.")
  if (is.null(k)) k <- cellcharter_autok_best_k(object)
  k <- .autok_validate_integer(k, "k")
  row <- object$best_index[object$best_index$k == k, , drop = FALSE]
  if (nrow(row) != 1L) stop("No fitted model is available for K = ", k, ".")
  readRDS(row$path[[1L]])$model
}

#' Prepare AutoK stability diagnostics
#'
#' Converts the repeated neighboring-K Fowlkes--Mallows scores in an AutoK fit
#' into tidy point and summary tables. A requested K outside the selectable
#' AutoK range (notably K = 1) is retained with an `NA` mean, because stability
#' cannot be estimated without fitted solutions on both neighboring K values.
#'
#' @param object A fitted `cellcharter_autok` object.
#' @param k_range Integer K values to display. Defaults to the inclusive range
#'   from one through the largest selectable K.
#' @param highlight_k K values to mark for interpretation, such as a published
#'   sensitivity analysis. Highlighting does not alter selection.
#' @return A list with `summary` (one row per displayed K) and `scores`
#'   (the saved repeated stability values in long form).
#' @examples
#' autok <- structure(
#'   list(
#'     candidate_k = 2:4,
#'     stability = rbind(c(0.7, 0.8), c(0.9, 0.95), c(0.6, 0.7)),
#'     stability_mean = c(0.75, 0.925, 0.65), best_k = 3L
#'   ),
#'   class = "cellcharter_autok"
#' )
#' cellcharter_autok_stability_data(autok, k_range = 1:4)$summary
#' @export
cellcharter_autok_stability_data <- function(
    object, k_range = NULL, highlight_k = integer()) {
  if (!inherits(object, "cellcharter_autok")) {
    stop("`object` must be a `cellcharter_autok` fit.")
  }
  if (is.null(object$stability) || is.null(object$stability_mean)) {
    stop("Stability diagnostics require at least two completed AutoK runs.")
  }
  if (is.null(k_range)) k_range <- seq.int(1L, max(object$candidate_k))
  k_range <- as.integer(k_range)
  if (!length(k_range) || anyNA(k_range) || any(k_range < 1L) ||
      anyDuplicated(k_range) || !identical(k_range, sort(k_range))) {
    stop("`k_range` must contain unique, increasing positive integers.")
  }
  highlight_k <- as.integer(highlight_k)
  if (anyNA(highlight_k) || any(highlight_k < 1L)) {
    stop("`highlight_k` must contain positive integers.")
  }

  candidate_k <- as.integer(object$candidate_k)
  stability <- object$stability
  if (nrow(stability) != length(candidate_k) ||
      length(object$stability_mean) != length(candidate_k)) {
    stop("The AutoK stability arrays do not match `candidate_k`.")
  }
  scores <- data.frame(
    k = rep(candidate_k, each = ncol(stability)),
    comparison = rep(seq_len(ncol(stability)), times = nrow(stability)),
    fowlkes_mallows = as.numeric(t(stability))
  )
  scores <- scores[scores$k %in% k_range, , drop = FALSE]

  mean_lookup <- stats::setNames(as.numeric(object$stability_mean), candidate_k)
  means <- unname(mean_lookup[as.character(k_range)])
  estimable <- !is.na(means)
  local_peak <- rep(FALSE, length(k_range))
  estimable_index <- which(estimable)
  if (length(estimable_index) >= 3L) {
    # Only interior points with two estimable neighbours are local peaks.
    for (i in estimable_index[-c(1L, length(estimable_index))]) {
      local_peak[[i]] <- means[[i]] > means[[i - 1L]] &&
        means[[i]] > means[[i + 1L]]
    }
  }
  summary <- data.frame(
    k = k_range,
    mean_fowlkes_mallows = means,
    estimable = estimable,
    local_peak = local_peak,
    selected = k_range == as.integer(object$best_k),
    highlighted = k_range %in% highlight_k
  )
  list(summary = summary, scores = scores)
}

#' Plot repeated AutoK stability across cluster counts
#'
#' Displays every saved neighboring-K Fowlkes--Mallows comparison, their mean
#' curve, local maxima, the independently selected K, and optional highlighted
#' K values. It never recomputes or substitutes a proxy for AutoK stability.
#'
#' @inheritParams cellcharter_autok_stability_data
#' @param title Plot title.
#' @return A `ggplot` object. The exact diagnostic tables are attached as the
#'   `cellcharter_stability_data` attribute for durable checkpointing.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' autok <- structure(
#'   list(
#'     candidate_k = 2:4,
#'     stability = rbind(c(0.7, 0.8), c(0.9, 0.95), c(0.6, 0.7)),
#'     stability_mean = c(0.75, 0.925, 0.65), best_k = 3L
#'   ),
#'   class = "cellcharter_autok"
#' )
#' cellcharter_plot_autok_stability(autok, k_range = 1:4, highlight_k = 4L)
#' @export
cellcharter_plot_autok_stability <- function(
    object, k_range = NULL, highlight_k = integer(),
    title = "CellCharter AutoK stability") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required to plot AutoK stability.")
  }
  diagnostic <- cellcharter_autok_stability_data(object, k_range, highlight_k)
  summary <- diagnostic$summary
  estimable <- summary[summary$estimable, , drop = FALSE]
  peaks <- summary[summary$local_peak, , drop = FALSE]

  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = diagnostic$scores,
      ggplot2::aes(x = k, y = fowlkes_mallows),
      colour = "grey65", alpha = 0.25, size = 1,
      position = ggplot2::position_jitter(width = 0.10, height = 0)
    ) +
    ggplot2::geom_line(
      data = estimable,
      ggplot2::aes(x = k, y = mean_fowlkes_mallows),
      colour = "#24557A", linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = estimable,
      ggplot2::aes(x = k, y = mean_fowlkes_mallows),
      colour = "#24557A", size = 2
    ) +
    ggplot2::geom_point(
      data = peaks,
      ggplot2::aes(x = k, y = mean_fowlkes_mallows),
      shape = 21, fill = "#E69F00", colour = "black", size = 3
    ) +
    ggplot2::scale_x_continuous(breaks = summary$k) +
    ggplot2::labs(
      x = "Number of clusters (K)",
      y = "Fowlkes--Mallows stability",
      title = title,
      subtitle = "Points: repeated neighboring-K comparisons; line: mean stability"
    ) +
    ggplot2::theme_bw(base_size = 16)

  marked <- summary[summary$estimable & (summary$selected | summary$highlighted), , drop = FALSE]
  if (nrow(marked)) {
    marked$label <- ifelse(
      marked$selected,
      paste0("selected K=", marked$k),
      paste0("K=", marked$k)
    )
    plot <- plot +
      ggplot2::geom_vline(
        data = marked, ggplot2::aes(xintercept = k, colour = selected),
        linetype = "dashed", linewidth = 0.55, show.legend = FALSE
      ) +
      ggplot2::geom_label(
        data = marked,
        ggplot2::aes(x = k, y = mean_fowlkes_mallows, label = label),
        nudge_y = 0.025, size = 5.6, show.legend = FALSE
      )
  }
  non_estimable <- summary[!summary$estimable, , drop = FALSE]
  if (nrow(non_estimable)) {
    plot <- plot + ggplot2::annotate(
      "text", x = non_estimable$k,
      y = min(diagnostic$scores$fowlkes_mallows, na.rm = TRUE),
      label = "not estimable", hjust = 0, angle = 90,
      size = 5.3, colour = "grey40"
    )
  }
  attr(plot, "cellcharter_stability_data") <- diagnostic
  plot
}

utils::globalVariables(c(
  "k", "fowlkes_mallows", "mean_fowlkes_mallows", "selected", "label"
))
