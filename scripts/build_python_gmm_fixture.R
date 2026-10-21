#!/usr/bin/env Rscript

# Build a small, permanent prediction-parity fixture from the checksum-pinned
# official PyCave K=11 model and its exact published X_cellcharter input.
source("R/spe_workflow.R")
source("R/torch_gmm.R")
source("R/autok.R")

args <- commandArgs(trailingOnly = TRUE)
archive <- if (length(args) >= 1L) args[[1L]] else
  "results/codex_mouse_spleen/input/BALBc-MRL_trVAE_CellCharter_autok.zip"
h5ad <- if (length(args) >= 2L) args[[2L]] else
  "results/codex_mouse_spleen/input/codex_mouse_spleen_clustered.h5ad"
output <- if (length(args) >= 3L) args[[3L]] else
  "inst/extdata/python_gmm_k11_fixture.rds"
if (!all(file.exists(c(archive, h5ad)))) {
  stop("Official AutoK ZIP and clustered H5AD are required to build the fixture.")
}

json_path <- tempfile(fileext = ".json")
on.exit(unlink(json_path), add = TRUE)
status <- system2(
  "python3",
  c("scripts/export_pycave_gmm_parameters.py", shQuote(archive), shQuote(json_path), "--k", "11")
)
if (!identical(status, 0L)) stop("Could not export the official PyCave parameters.")
state <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
means <- do.call(rbind, lapply(state$means, unlist, use.names = FALSE))
precision_cholesky <- lapply(
  state$precisions_cholesky,
  function(component) do.call(rbind, lapply(component, unlist, use.names = FALSE))
)
weights <- unlist(state$component_probs, use.names = FALSE)

published <- DelayedArray::t(HDF5Array::HDF5Array(h5ad, "obsm/X_cellcharter"))
base::set.seed(20230817L)
rows <- sort(sample.int(nrow(published), 512L))
x <- as.matrix(published[rows, , drop = FALSE])
obs <- .cc_h5ad_dataframe(h5ad, "obs")
published_labels <- as.integer(obs$cluster_cellcharter[rows]) - 1L

# Reproduce PyCave 3.2.1's jit_log_normal formula independently in base R:
# each row-vector residual is right-multiplied by the stored upper-triangular
# precision Cholesky factor before its squared norm is taken.
log_scores <- vapply(seq_len(nrow(means)), function(component) {
  centered <- sweep(x, 2L, means[component, ], FUN = "-")
  mahalanobis <- rowSums((centered %*% precision_cholesky[[component]])^2)
  log(weights[[component]]) + sum(log(diag(precision_cholesky[[component]]))) -
    0.5 * (ncol(x) * log(2 * pi) + mahalanobis)
}, numeric(nrow(x)))
expected <- max.col(log_scores, ties.method = "first") - 1L

candidate_covariances <- list(
  transpose_times = lapply(precision_cholesky, function(value) solve(crossprod(value))),
  times_transpose = lapply(precision_cholesky, function(value) solve(tcrossprod(value)))
)
candidate_predictions <- lapply(candidate_covariances, function(covariances) {
  model <- structure(list(
    means = means, covariances = covariances, weights = weights,
    n_features = ncol(means), n_clusters = nrow(means),
    covariance_regularization = 1e-12
  ), class = "torch_gmm")
  torch_gmm_predict(model, x, batch_size = 128L, device = "cpu")
})
fmi <- vapply(
  candidate_predictions,
  function(value) cellcharter_fowlkes_mallows(value, expected),
  numeric(1)
)
agreement <- fmi > 1 - 1e-7
if (sum(agreement) != 1L) {
  accuracy <- vapply(candidate_predictions, function(value) mean(value == expected), numeric(1))
  stop(
    "Could not uniquely recover the PyCave covariance convention; FMI/accuracies: ",
    paste(names(fmi), signif(fmi, 6), signif(accuracy, 6), collapse = ", ")
  )
}
convention <- names(agreement)[agreement]
fixture <- list(
  x = x, expected_labels = expected,
  model = structure(list(
    means = means, covariances = candidate_covariances[[convention]], weights = weights,
    n_features = ncol(means), n_clusters = nrow(means),
    covariance_regularization = 1e-12
  ), class = "torch_gmm"),
  rows = rows,
    covariance_recovery = convention,
  published_label_fmi = cellcharter_fowlkes_mallows(expected, published_labels),
  provenance = list(
    cellcharter_version = "0.2.0",
    cellcharter_commit = "2a8fdde8213e32e5a00f1f9571bffd3dbb825a79",
    backend = "PyCave 3.2.1 / PyTorch 1.12.1",
    pycave_commit = "6df56fa15797940a2063f3e70b3a7641f592a9f6",
    archive_sha256 = cellcharter_file_sha256(archive),
    h5ad_sha256 = cellcharter_file_sha256(h5ad),
    fixture_seed = 20230817L,
    agreement_metric = "Fowlkes-Mallows (label-permutation invariant)",
    expected_agreement = 1
  )
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, output, compress = "xz", version = 2L)
message("Wrote pinned Python GMM prediction fixture: ", output)
