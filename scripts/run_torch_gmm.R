#!/usr/bin/env Rscript

# Command-line wrapper deliberately keeps model fitting resumable. The calling
# report writes an RDS checkpoint after each selected K/run, so a scheduler
# restart does not repeat completed models.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: run_torch_gmm.R <features.rds> <n_clusters> <model.rds>")
}
source("R/torch_gmm.R")
features <- readRDS(args[[1]])
model <- torch_gmm_fit(features, as.integer(args[[2]]), device = "auto")
saveRDS(model, args[[3]])
