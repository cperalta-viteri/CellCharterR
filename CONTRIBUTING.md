# Contributing to CellCharterR

Please open an issue before a substantial change so its scientific scope and
reference target can be agreed first. Contributions should be focused, use
explicit package namespaces, preserve `SpatialExperiment` spot identity,
and avoid materialising cohort-scale assays or graphs.

Add deterministic tests for changed behavior. Run the smallest affected tests
first and then `R CMD check --no-manual` on the source tarball. Methodological
changes must identify the pinned CellCharter source behavior, document any
intentional deviation, and include a validation fixture where feasible.

Do not commit downloaded full datasets, generated caches, credentials, or
private data. By contributing, you agree that your contribution is distributed
under the package's MIT license.
