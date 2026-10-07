# rCellCharter

`rCellCharter` is a disk-backed R implementation of the **post-representation**
CellCharter workflow. It accepts a supplied embedding, constructs sample-local
spatial graphs, aggregates exclusive graph-hop features, selects the number of
domains with CellCharter-style AutoK, and fits a full-covariance GMM. It
intentionally does not reimplement scVI or trVAE.

## Why the GMM uses Torch directly

The GMM is an expectation-maximisation (EM) algorithm. Its E-step evaluates
Gaussian probabilities and its M-step has closed-form weighted means and
covariances. These are large matrix operations, so this implementation uses
R Torch tensors and CUDA when available. A neural-network training framework is
not used in the inner loop because EM has closed-form parameter updates.

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

The canonical CODEX report source is in `reports/`; a synchronized packaged
copy is under `inst/rmarkdown/`.

## Reproducing the official tutorial analyses

Download the checksum-pinned public inputs, then run the checkpointed workflow:

```bash
Rscript scripts/download_codex_mouse_spleen.R
Rscript scripts/run_codex_mouse_spleen.R
```

By default, inputs and durable stage checkpoints live under
`results/codex_mouse_spleen/`. Set `CODEX_INPUT_DIR` or `CODEX_RESULTS_DIR` to
use another layout. `CODEX_STAGES` can resume any comma-separated subset of
`representation,graph,aggregation,autok,downstream,report`.

The primary path treats the official trVAE representation as immutable input,
rebuilds sample-local Delaunay graphs with R Qhull, and follows pinned
CellCharter 0.2.0 by applying one global 99th-percentile link cutoff. The
released H5AD instead reflects a separate 99th-percentile cutoff within each
sample. The workflow saves both variants: the global graph drives primary
inference, while the R per-sample reconstruction is verified edge-for-edge
against the published adjacency and reproduces the official 40-dimensional
aggregate. This source-versus-artifact discrepancy is visible in the manifest
and report rather than silently selecting whichever path matches.

AutoK declares candidate `K = 2..20` in advance, fits the expanded `K = 1..21`
grid with at most 10 repeated full-covariance R Torch fits, and uses the
CellCharter 0.2.0 neighboring-K Fowlkes–Mallows rule. The independently selected
K drives all primary maps and statistics. K = 11 is kept only as a separately
labelled tutorial sensitivity analysis. Generated data and checkpoints are
ignored by Git and excluded from the source package.

The permanent K=11 prediction fixture in `inst/extdata/` was generated from
the pinned PyCave model and exact published aggregate with:

```bash
Rscript scripts/build_python_gmm_fixture.R
```

That one-time fixture builder uses the Python standard library only to decode
the official PyTorch archive; Python is not an rCellCharter runtime dependency.

The CosMx workflow follows the official tutorial's supplied scVI representation
but deliberately rebuilds sample-local graphs and three-hop niche features in R:

```bash
COSMX_H5AD=/path/to/cosmx_human_nsclc_clustered.h5ad \
  Rscript scripts/run_cosmx_human_nsclc.R
Rscript -e 'rmarkdown::render("inst/rmarkdown/cosmx_human_nsclc_report.Rmd")'
```
