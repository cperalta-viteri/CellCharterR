# rCellCharter

`rCellCharter` is a disk-backed R implementation of the **post-representation**
CellCharter workflow. It accepts PCA (or another explicitly supplied embedding),
aggregates features over a spatial graph, then finds spatial domains with a
full-covariance GMM. It intentionally does not reimplement scVI.

## Why the GMM uses Torch directly

The GMM is an expectation-maximisation (EM) algorithm. Its E-step evaluates
Gaussian probabilities and its M-step has closed-form weighted means and
covariances. These are large matrix operations, so this implementation uses
R Torch tensors and CUDA when available. `luz` is intentionally a suggested
dependency, rather than a requirement for the inner loop: Luz is designed for
gradient-based neural-network training and adding an optimiser/training loop
to closed-form EM would make it slower, not faster.

The GMM never materialises the whole feature matrix. Each iteration:

1. reads an HDF5/DelayedArray row chunk;
2. transfers the chunk to CPU/GPU Torch memory;
3. calculates responsibilities and sufficient statistics;
4. releases the chunk and advances to the next one;
5. updates only the small K-by-P means and K P-by-P covariance matrices.

## Minimal example

```r
library(HDF5Array)
library(rCellCharter)

# `X_cellcharter` has observations in rows and lives on disk.
gmm <- torch_gmm_fit(
  X_cellcharter,
  n_clusters = 8,
  batch_size = 8192,
  max_iter = 100,
  device = "auto", # CUDA if available, CPU otherwise
  seed = 42
)

# Prediction is also chunked; labels are the only full-length in-memory result.
domains <- torch_gmm_predict(gmm, X_cellcharter, batch_size = 8192)
```

See `inst/rmarkdown/` for the two fully reproducible tutorial reports.

## Reproducing the official tutorial analyses

The runnable CODEX workflow is intentionally separate from the package API:

```bash
Rscript scripts/run_codex_mouse_spleen.R
Rscript -e 'rmarkdown::render("inst/rmarkdown/codex_mouse_spleen_report.Rmd")'
```

Set `CODEX_H5AD` if the official tutorial file is outside `../data/`. The
workflow records the PCA model, sample-local graphs, HDF5 niche representation,
GMM model, and compact label/metadata table. Generated results are excluded
from the source package; they are needed only to render the report.

The CosMx workflow follows the official tutorial's supplied scVI representation
but deliberately rebuilds sample-local graphs and three-hop niche features in R:

```bash
COSMX_H5AD=/path/to/cosmx_human_nsclc_clustered.h5ad \
  Rscript scripts/run_cosmx_human_nsclc.R
Rscript -e 'rmarkdown::render("inst/rmarkdown/cosmx_human_nsclc_report.Rmd")'
```
