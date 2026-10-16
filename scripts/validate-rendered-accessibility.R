#!/usr/bin/env Rscript

# Inventory rendered plot alternatives and create grayscale/protanopia/
# deuteranopia simulations. This validates documentation artifacts only; it
# does not read or recompute scientific data.

args <- commandArgs(trailingOnly = TRUE)
html_paths <- if (length(args)) args else c(
  "docs/CellCharterR-beginner-workflow.html",
  "docs/CellCharterR-mouse-spleen-reproducibility.html",
  "docs/CellCharterR-differential-neighborhood-networks.html",
  "docs/CellCharterR-h5ad-to-sce.html"
)
output_dir <- "results/accessibility_validation"

required_packages <- c("base64enc", "digest", "png", "xml2")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing accessibility package(s): ", paste(missing_packages, collapse = ", "))
}
if (any(!file.exists(html_paths))) {
  stop("Rendered HTML is absent: ", paste(html_paths[!file.exists(html_paths)], collapse = ", "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cvd_matrices <- list(
  protanopia = matrix(c(
    0.152286, 1.052583, -0.204868,
    0.114503, 0.786281, 0.099216,
    -0.003882, -0.048116, 1.051998
  ), nrow = 3L, byrow = TRUE),
  deuteranopia = matrix(c(
    0.367322, 0.860646, -0.227968,
    0.280085, 0.672501, 0.047413,
    -0.011820, 0.042940, 0.968881
  ), nrow = 3L, byrow = TRUE)
)

simulate_rgb <- function(image, transform) {
  dimensions <- dim(image)
  rgb <- matrix(image[, , seq_len(3L), drop = FALSE], ncol = 3L)
  simulated <- rgb %*% t(transform)
  simulated[simulated < 0] <- 0
  simulated[simulated > 1] <- 1
  output <- array(simulated, dim = c(dimensions[[1L]], dimensions[[2L]], 3L))
  if (dimensions[[3L]] == 4L) {
    output <- abind_alpha(output, image[, , 4L])
  }
  output
}

abind_alpha <- function(rgb, alpha) {
  output <- array(0, dim = c(dim(rgb)[1:2], 4L))
  output[, , seq_len(3L)] <- rgb
  output[, , 4L] <- alpha
  output
}

records <- list()
record_index <- 0L
for (html_path in html_paths) {
  document <- xml2::read_html(html_path)
  images <- xml2::xml_find_all(document, "//img[starts-with(@src, 'data:image/png;base64,')]")
  for (index in seq_along(images)) {
    alt <- xml2::xml_attr(images[[index]], "alt")
    if (is.na(alt) || !nzchar(trimws(alt))) {
      stop("Rendered PNG lacks alternative text in ", html_path, " at image ", index)
    }
    encoded <- sub(
      "^data:image/png;base64,", "", xml2::xml_attr(images[[index]], "src")
    )
    raw_png <- base64enc::base64decode(encoded)
    source_png <- tempfile(fileext = ".png")
    writeBin(raw_png, source_png)
    image <- png::readPNG(source_png)
    if (length(dim(image)) != 3L || dim(image)[[3L]] < 3L) {
      stop("Rendered image is not an RGB/RGBA PNG: ", html_path, " image ", index)
    }

    stem <- paste0(
      tools::file_path_sans_ext(basename(html_path)), "-", sprintf("%02d", index)
    )
    grayscale <- 0.2126 * image[, , 1L] +
      0.7152 * image[, , 2L] + 0.0722 * image[, , 3L]
    grayscale_rgb <- array(rep(grayscale, 3L), dim = c(dim(image)[1:2], 3L))
    if (dim(image)[[3L]] == 4L) {
      grayscale_rgb <- abind_alpha(grayscale_rgb, image[, , 4L])
    }
    png::writePNG(grayscale_rgb, file.path(output_dir, paste0(stem, "-grayscale.png")))
    for (mode in names(cvd_matrices)) {
      png::writePNG(
        simulate_rgb(image, cvd_matrices[[mode]]),
        file.path(output_dir, paste0(stem, "-", mode, ".png"))
      )
    }
    record_index <- record_index + 1L
    records[[record_index]] <- data.frame(
      html = basename(html_path), image_index = index, alt = alt,
      width_px = dim(image)[[2L]], height_px = dim(image)[[1L]],
      grayscale = paste0(stem, "-grayscale.png"),
      protanopia = paste0(stem, "-protanopia.png"),
      deuteranopia = paste0(stem, "-deuteranopia.png"),
      stringsAsFactors = FALSE
    )
    unlink(source_png)
  }
}
inventory <- do.call(rbind, records)
if (is.null(inventory) || !nrow(inventory)) {
  stop("No embedded PNG figures were found in the rendered reports.")
}
utils::write.csv(inventory, file.path(output_dir, "figure_inventory.csv"), row.names = FALSE)
result <- list(
  status = "PASS",
  validation_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  html_sha256 = stats::setNames(vapply(
    html_paths,
    function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ),
    character(1)
  ), basename(html_paths)),
  figure_count = nrow(inventory),
  all_alt_text_nonempty = TRUE,
  simulations = c("grayscale", "protanopia", "deuteranopia"),
  note = paste(
    "Every embedded PNG has non-empty alternative text and was rendered to",
    "grayscale, protanopia, and deuteranopia simulations for manual review.",
    "Some scientific panels encode values primarily by colour; inspect the",
    "underlying saved result objects when exact numeric values are required."
  )
)
saveRDS(result, file.path(output_dir, "accessibility_validation.rds"), version = 3L)
cat("Rendered accessibility validation: PASS\n")
cat("Figures inventoried:", nrow(inventory), "\n")
cat("Artifacts:", output_dir, "\n")
