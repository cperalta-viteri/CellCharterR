# Support and maintenance

CellCharterR is an independent implementation and is not supported by the
upstream CellCharter authors. Reproducible bug reports should be filed at the
[`CellCharterR` issue tracker](https://github.com/cperalta-viteri/CellCharterR/issues).
If the tracker is unavailable before the repository is published, contact the
package maintainer listed in `DESCRIPTION`.

A useful report includes the smallest reproducible example, `sessionInfo()`,
the CellCharterR version, CPU thread count, random seed, input dimensions, and
the complete error message. Do not attach private spatial data. For large-data
failures, include checkpoint provenance and measured memory use rather than the
underlying dataset.

The maintainer reviews correctness and security fixes first. Public functions
will not be removed without a documented deprecation period and an entry in
`NEWS.md`. Reference-parity behavior may be retained alongside a recommended
modern alternative when changing it would invalidate a published reproduction.
