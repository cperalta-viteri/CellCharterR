# CellCharterR 0.99.0

- Added an exported SCE-first workflow covering sample-local graph
  construction, neighborhood aggregation, AutoK/GMM clustering, domain-label
  attachment, cell-type enrichment, neighborhood enrichment, and differential
  neighborhood enrichment.
- Added strict ordered-cell identity checks, `colPairs()` graph storage,
  synchronized subsetting tests, and stale-stage detection.
- Added a portable 1,000-cell CODEX teaching SCE and exact miniature H5AD import
  fixtures under CC BY 4.0 attribution.
- Added a CPU-only beginner vignette over K=2 through K=10, an H5AD-to-SCE
  conversion vignette, and the preserved full mouse-spleen reproducibility
  report as immediately viewable HTML.
- Added a movable, separately distributed full mouse-spleen SCE checkpoint with
  no H5AD or HDF5Array dependency.
- Incorporated the second-round methodological corrections and tests for graph,
  aggregation, AutoK, full-covariance GMM, enrichment, differential inference,
  provenance, and bounded-memory behavior.
