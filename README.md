# CellCharterR

## Vignettes and reports

- [Beginner SCE-first workflow](https://cperalta-viteri.github.io/CellCharterR/CellCharterR-beginner-workflow.html)
  runs K=2 through K=10 on the packaged 1,000-cell teaching object, on CPU.
- [H5AD download and SCE conversion](https://cperalta-viteri.github.io/CellCharterR/CellCharterR-h5ad-to-sce.html)
  shows checksum-verified download and conversion code without executing the
  full download during package builds.
- [Full CODEX mouse-spleen reproducibility report](https://cperalta-viteri.github.io/CellCharterR/CellCharterR-mouse-spleen-reproducibility.html)
  is the preserved, precomputed full-data report. It is not rerun during build,
  check, examples, or CI.

`CellCharterR` is a disk-backed R implementation of the **post-representation**
CellCharter workflow. It accepts a supplied embedding, constructs sample-local
spatial graphs, aggregates exclusive graph-hop features, selects the number of
domains with CellCharter-style AutoK, and fits a full-covariance GMM. It
intentionally does not reimplement scVI or trVAE.

CellCharterR is an independent R implementation and is not authored, endorsed,
or maintained by the upstream CellCharter project. The package pins and cites
the upstream version used for methodological comparisons.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Make both Bioconductor and CRAN dependencies visible to the installer.
options(repos = BiocManager::repositories())

# R torch is distributed through CRAN, not Bioconductor.
install.packages("torch")

# Install the Linux CPU LibTorch runtime. Restart R if prompted.
torch::install_torch(type = "cpu")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github(
  "cperalta-viteri/CellCharterR",
  build_vignettes = TRUE,
  dependencies = TRUE
)

library(CellCharterR)
torch_gmm_device("cpu")
```

See [INSTALL](INSTALL) for the supported Linux CPU target and smoke test.

## Why the GMM uses Torch directly

The GMM is an expectation-maximisation (EM) algorithm. Its E-step evaluates
Gaussian probabilities and its M-step has closed-form weighted means and
covariances. These are large matrix operations, so this implementation uses
R Torch tensors. The supported package examples and checks run on Linux CPU. A
neural-network training framework is not used in the inner loop because EM has
closed-form parameter updates.

The GMM never materialises the whole feature matrix. Each iteration:

1. reads an HDF5/DelayedArray row chunk;
2. transfers the chunk to CPU/GPU Torch memory;
3. calculates responsibilities and sufficient statistics;
4. releases the chunk and advances to the next one;
5. updates only the small K-by-P means and K P-by-P covariance matrices.

## Minimal example

```r
library(CellCharterR)
data("codex_mouse_spleen_tiny", package = "CellCharterR")

# Build the spatial graph and exclusive one-hop neighborhood representation.
sce <- cellcharter_sce_graph(codex_mouse_spleen_tiny)
sce <- cellcharter_sce_aggregate(sce, n_layers = 1L)
SingleCellExperiment::reducedDimNames(sce)
```

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

The second-review Phase-0 K=4/K=11 deliverables are generated with:

```bash
Rscript scripts/build_second_review_deliverables.R
Rscript -e 'rmarkdown::render("reports/CellCharterR_mouse_spleen_AutoK_K11.Rmd", output_dir="results/CellCharterR_mouse_spleen_v2/report")'
```

They are stored under `results/CellCharterR_mouse_spleen/`. The original
`results/codex_mouse_spleen/` tree is retained as the immutable first-round
checkpoint bundle on which these recomputed downstream analyses depend.

The permanent K=11 prediction fixture in `inst/extdata/` was generated from
the pinned PyCave model and exact published aggregate with:

```bash
Rscript scripts/build_python_gmm_fixture.R
```

That one-time fixture builder uses the Python standard library only to decode
the official PyTorch archive; Python is not a CellCharterR runtime dependency.

The CosMx workflow follows the official tutorial's supplied scVI representation
but deliberately rebuilds sample-local graphs and three-hop niche features in R:

```bash
COSMX_H5AD=/path/to/cosmx_human_nsclc_clustered.h5ad \
  Rscript scripts/run_cosmx_human_nsclc.R
Rscript -e 'rmarkdown::render("inst/rmarkdown/cosmx_human_nsclc_report.Rmd")'
```
