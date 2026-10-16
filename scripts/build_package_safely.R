#!/usr/bin/env Rscript

# Build CellCharterR from a small, explicit package-source staging tree.
#
# `R CMD build .` copies the repository to R's temporary directory before it
# applies .Rbuildignore. In this workspace that can duplicate the very large
# transitionR H5AD even though transitionR is ignored. This launcher copies
# only standard R-package inputs first, so R never sees private data/results.

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- arguments[startsWith(arguments, "--file=")]
if (length(script_argument) != 1L) {
    stop("Run this file with Rscript, for example `Rscript scripts/build_package_safely.R`.")
}
script_path <- sub("^--file=", "", script_argument)
package_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
if (!file.exists(file.path(package_root, "DESCRIPTION"))) {
    stop("Could not locate package DESCRIPTION relative to this script.")
}

top_level_files <- c(
    "DESCRIPTION", "NAMESPACE", "LICENSE", "LICENCE", "COPYING",
    "NEWS", "NEWS.md", "README", "README.md", ".Rbuildignore", ".Rinstignore",
    "configure", "configure.ac", "cleanup"
)
package_directories <- c(
    "R", "data", "demo", "exec", "inst", "man", "po", "src",
    "tests", "vignettes"
)
source_entries <- c(top_level_files, package_directories)
source_entries <- source_entries[file.exists(file.path(package_root, source_entries))]

stage_parent <- tempfile("cellcharterr-package-source-")
stage_package <- file.path(stage_parent, basename(package_root))
dir.create(stage_package, recursive = TRUE)
on.exit(unlink(stage_parent, recursive = TRUE, force = TRUE), add = TRUE)

copied <- file.copy(
    from = file.path(package_root, source_entries),
    to = stage_package,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
)
if (!all(copied)) {
    stop("Failed to construct the safe package-source staging tree.")
}

message("Building from data-free staging tree: ", stage_package)
old_directory <- setwd(stage_parent)
on.exit(setwd(old_directory), add = TRUE)
status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "build", shQuote(basename(stage_package)))
)
if (!identical(status, 0L)) stop("R CMD build failed with status ", status, ".")

description <- read.dcf(file.path(package_root, "DESCRIPTION"))
tarball_name <- sprintf(
    "%s_%s.tar.gz", description[1L, "Package"], description[1L, "Version"]
)
built_tarball <- file.path(stage_parent, tarball_name)
if (!file.exists(built_tarball)) stop("R CMD build did not create `", tarball_name, "`.")
destination <- file.path(package_root, tarball_name)
if (!file.copy(built_tarball, destination, overwrite = TRUE, copy.date = TRUE)) {
    stop("Could not publish the source tarball at `", destination, "`.")
}
message("Built: ", destination)
