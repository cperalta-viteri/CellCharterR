#!/usr/bin/env Rscript

# Download immutable public inputs used by the CODEX reproduction. Existing
# files are reused only after their published MD5 checksum has been verified.
input_dir <- Sys.getenv(
  "CODEX_INPUT_DIR",
  file.path("results", "codex_mouse_spleen", "input")
)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

artifacts <- data.frame(
  filename = c(
    "codex_mouse_spleen.h5ad",
    "codex_mouse_spleen_clustered.h5ad",
    "BALBc-MRL_trVAE_CellCharter_autok.zip"
  ),
  file_id = c("38538101", "46832722", "46833700"),
  md5 = c(
    "7d5ad7727ac214bcd61a55b2572ceb36",
    "b3d50ea1181b59b7b8fc0522633badb2",
    "3833ec2ca271e2a57de2fbe2322b64b1"
  ),
  sha256 = c(
    "b386cb72b20d483d3c84675bffcdd22cae1e4db6e4333271bc2c3bdf4ecef946",
    "a23434e23cbfb421b2f51593eb71fe5d642c545aca705e4607a4b2a081047bd6",
    "c7ea0c53d8c103c9752e343142e4ee632090faf661fb1c482760bef1aefb6e01"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(artifacts))) {
  destination <- file.path(input_dir, artifacts$filename[[i]])
  if (file.exists(destination) &&
      unname(tools::md5sum(destination)) == artifacts$md5[[i]] &&
      digest::digest(file = destination, algo = "sha256", serialize = FALSE) ==
        artifacts$sha256[[i]]) {
    message("Verified existing input: ", destination)
    next
  }

  temporary <- paste0(destination, ".partial")
  on.exit(unlink(temporary), add = TRUE)
  url <- paste0("https://api.figshare.com/v2/file/download/", artifacts$file_id[[i]])
  message("Downloading ", url, " -> ", destination)
  utils::download.file(url, temporary, mode = "wb", method = "libcurl", quiet = FALSE)
  observed <- unname(tools::md5sum(temporary))
  observed_sha256 <- digest::digest(
    file = temporary, algo = "sha256", serialize = FALSE
  )
  if (!identical(observed, artifacts$md5[[i]]) ||
      !identical(observed_sha256, artifacts$sha256[[i]])) {
    stop("Checksum mismatch for ", artifacts$filename[[i]],
         ": expected MD5/SHA-256 ", artifacts$md5[[i]], "/",
         artifacts$sha256[[i]], ", observed ", observed, "/", observed_sha256)
  }
  if (!file.rename(temporary, destination)) {
    stop("Could not move verified download into place: ", destination)
  }
}

message("All CODEX inputs are present and checksum-verified in ", input_dir)
