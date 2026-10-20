# CellCharterR 0.99.0

- Optimized the production C++ CPU GMM path around one shared materialized
  representation, globally centered full-covariance EM, one cached N-by-K
  posterior per concurrent fit, and blocked BLAS `DSYRK` covariance updates.
  Candidate fits now start in descending K order to reduce the scheduling tail
  while results retain canonical caller order. The production path remains
  CPU-only; GPU clustering remains experimental. Set OpenBLAS, OpenMP, and MKL
  thread counts to one when using outer or per-fit parallelism to avoid nested
  oversubscription.
- Full optimized mouse-spleen qualification is still running. The comparable
  pre-optimization compiled-C++ run used 707,474 observations by 40 features,
  `K = 2..20`, five rounds, and 105 fits, and took 2300.1 seconds. Its optimized
  wall time remains `FINAL_OPTIMIZED_SECONDS`; this is a searchable placeholder,
  not a result or speedup claim.
- Made R Torch an optional experimental dependency. Loading CellCharterR or a
  legacy disk-backed SPE now remains a pure Bioconductor/HDF5 operation and
  cannot trigger a Lantern/LibTorch download; Torch entry points fail with an
  actionable installation message when the optional package is absent.
- Avoided repeated hyperslab reads of large HDF5 variable-length string
  datasets during import preflight. This prevents an rhdf5 native crash seen
  while reconstructing the 3.5-million-spot legacy DeepSpot-M SPE.
- Reorganized the public documentation around three workflows: a primary
  precomputed full mouse-spleen downstream tutorial, a large-H5AD disk-backed
  import tutorial, and the preserved experimental full R AutoK reproduction.
- Marked the R Torch GMM/AutoK fitting path experimental for multi-million-spot
  datasets while retaining the supported graph, aggregation, downstream
  statistics, zonation, and plotting APIs. Downstream functions now accept
  cluster labels directly from `colData` or `reducedDims`, without requiring a
  CellCharterR clustering-stage record.
- Added an exported SpatialExperiment-first workflow covering sample-local graph
  construction, neighborhood aggregation, AutoK/GMM clustering, domain-label
  attachment, cell-type enrichment, neighborhood enrichment, and differential
  neighborhood enrichment.
- Added strict ordered-cell identity checks, `colPairs()` graph storage,
  synchronized subsetting tests, and stale-stage detection.
- Added a portable 1,000-spot CODEX teaching SPE and exact miniature H5AD import
  fixtures under CC BY 4.0 attribution.
- Added a CPU-only beginner vignette over K=2 through K=10, an H5AD-to-SCE
  conversion vignette, and the preserved full mouse-spleen reproducibility
  report as immediately viewable HTML.
- Added a movable, separately distributed full mouse-spleen SPE checkpoint with
  no H5AD or HDF5Array dependency.
- Added `cellcharter_plot_graph_diff()` for deterministic two-panel directed
  differential-neighborhood networks, together with a full mouse-spleen K=10
  illustrative vignette and rendered HTML.
- Added block-aware authors-style AUCell zonation with disk-backed score
  access, compact patient-level results, optional spot-level zone labels, and
  plotting for group comparisons, metadata-stratified comparisons, individual
  trajectories, group averages, and BH-adjusted significance annotations.
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
