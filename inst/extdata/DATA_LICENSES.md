# Data licenses and attribution

- Original CODEX mouse-spleen data: Goltsev et al. (2018), provided by Mendeley
  Data at DOI `10.17632/zjnpwh8m5b.1`, under
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); publication DOI
  `10.1016/j.cell.2018.07.010`. The tutorial H5AD was retrieved from Figshare
  file 38538101 and verified on 2026-10-12 with SHA-256
  `b386cb72b20d483d3c84675bffcdd22cae1e4db6e4333271bc2c3bdf4ecef946`.
- CellCharter clustered mouse-spleen H5AD and AutoK archive: provided by
  Figshare at DOI `10.6084/m9.figshare.25040192.v2`, under
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The record contains
  file 46832722 (`codex_mouse_spleen_clustered.h5ad`; retrieved and verified
  2026-10-12; SHA-256
  `a23434e23cbfb421b2f51593eb71fe5d642c545aca705e4607a4b2a081047bd6`)
  and file 46833700 (`BALBc-MRL_trVAE_CellCharter_autok.zip`).
- CellCharter method: Varrone et al. (2024), DOI
  `10.1038/s41588-023-01588-4`.

`codex_mouse_spleen_tiny` is a deterministic 1,000-spot CODEX derivative distributed
under the same CC BY 4.0 attribution terms. It is intended for teaching and
software validation, not biological inference.

`mouse_spleen_k4_enrichment.rds` and `mouse_spleen_k4_diff_nhood.rds` are
historical compact derived outputs from that same CODEX data. They are
distributed under the same CC BY 4.0 attribution terms and retain their
CellCharterR/reference provenance within the objects or the archived report.

`mouse_spleen_cpp_k10_labels.rds`, `mouse_spleen_cpp_k10_enrichment.rds`, and
`mouse_spleen_cpp_k10_diff_nhood.rds` are the corresponding compact outputs
from the completed compiled-C++ AutoK primary result. They use the same CODEX
attribution terms and record the backend, ordered-cell identity, graph contract,
and pinned CellCharter reference in the label artifact and rendered report.
