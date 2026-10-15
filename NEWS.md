# CellCharterR 0.99.0

- Added an exported SpatialExperiment-first workflow covering sample-local graph
  construction, neighborhood aggregation, AutoK/GMM clustering, domain-label
  attachment, cell-type enrichment, neighborhood enrichment, and differential
  neighborhood enrichment.
- Added strict ordered-cell identity checks, `colPairs()` graph storage,
  synchronized subsetting tests, and stale-stage detection.
- Added a portable 1,000-segmented-cell CODEX teaching SPE and exact miniature H5AD import
  fixtures under CC BY 4.0 attribution.
- Added a CPU-only beginner vignette over K=2 through K=10, an H5AD-to-SCE
  conversion vignette, and the preserved full mouse-spleen reproducibility
  report as immediately viewable HTML.
- Added a movable, separately distributed full mouse-spleen SPE checkpoint with
  no H5AD or HDF5Array dependency.
- Added `cellcharter_plot_graph_diff()` for deterministic two-panel directed
  differential-neighborhood networks, together with a full mouse-spleen K=4
  illustrative vignette and rendered HTML.
- Incorporated the second-round methodological corrections and tests for graph,
  aggregation, AutoK, full-covariance GMM, enrichment, differential inference,
  provenance, and bounded-memory behavior.
- Made AutoK resume restore committed stability state without replaying every
  historical Fowlkes--Mallows comparison. Representation caching now begins
  only for genuinely unfinished work, and an interrupted final run is
  recoverably quarantined and restarted from its first K.
- Simplified the public workflow to neutral `cellcharter_*()` names. Numerical
  matrix kernels now use internal `.cc_*()` names; the former
  `cellcharter_sce_*()` entry points remain as deprecated compatibility wrappers.
