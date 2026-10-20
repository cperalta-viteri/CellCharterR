#!/usr/bin/env Rscript

# A deliberately tiny live compiled-CPU fit used by Bioconductor-devel CI.
# Torch is optional and must not be initialized by the supported CPU workflow,
# so the enclosing workflow can also reject network access during validation.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "cpu-smoke.rds"

if (!identical(Sys.getenv("CELLCHARTERR_CPU_ONLY"), "true")) {
    stop("CELLCHARTERR_CPU_ONLY=true is required")
}

suppressPackageStartupMessages(library(CellCharterR))
set.seed(481L)
x <- rbind(
    matrix(stats::rnorm(48L * 3L, -2, 0.25), ncol = 3L),
    matrix(stats::rnorm(48L * 3L, 2, 0.25), ncol = 3L)
)

fit <- cpp_gmm_fit(
    x,
    n_clusters = 2L,
    seed = 481L,
    max_iter = 25L,
    num_threads = 1L
)
labels <- cpp_gmm_predict(fit, x, batch_size = 32L)

result <- list(
    passed = is.finite(fit$nll) &&
        length(labels) == nrow(x) && !anyNA(labels),
    device = fit$runtime$device,
    backend = fit$runtime$backend,
    n = nrow(x),
    p = ncol(x),
    k = fit$n_clusters,
    negative_log_likelihood = fit$nll,
    labels_complete = !anyNA(labels),
    session_info = utils::sessionInfo()
)
saveRDS(result, output)
if (!isTRUE(result$passed) || !identical(result$device, "cpu") ||
    !identical(result$backend, "compiled-cpp-double")) {
    stop("The live compiled-CPU smoke test failed")
}
