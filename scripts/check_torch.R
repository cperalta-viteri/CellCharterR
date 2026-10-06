#!/usr/bin/env Rscript

# Small, fast pre-flight check. Run this before a full tutorial to avoid wasting
# a long job on a missing libtorch installation or an unavailable GPU.
library(torch)
cat("R Torch version:", as.character(packageVersion("torch")), "\n")
cat("CUDA visible to R Torch:", cuda_is_available(), "\n")
if (cuda_is_available()) {
  cat("Device:", torch_device("cuda")$type, "\n")
}
