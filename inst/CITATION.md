# Citing CellCharterR

If CellCharterR contributes to an analysis, please cite both this R
implementation and the original CellCharter method.

## CellCharterR implementation

Peralta Viteri, C. (2026). *CellCharterR: Disk-Backed R Implementation of
CellCharter Spatial Domains*. R package version 0.99.0.
https://github.com/cperalta-viteri/CellCharterR

```bibtex
@Manual{CellCharterR,
  title  = {CellCharterR: Disk-Backed R Implementation of CellCharter Spatial Domains},
  author = {Christian Peralta Viteri},
  year   = {2026},
  note   = {R package version 0.99.0},
  url    = {https://github.com/cperalta-viteri/CellCharterR}
}
```

After installation, the same package citation can be retrieved in R:

```r
citation("CellCharterR")
```

## Original CellCharter method

Varrone, M. et al. CellCharter reveals spatial cell niches associated with
tissue remodeling and cell plasticity. *Nature Genetics* **56**, 74–84 (2024).
https://doi.org/10.1038/s41588-023-01588-4

```bibtex
@Article{Varrone2024CellCharter,
  title   = {CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity},
  author  = {Varrone, Marco and Tavernari, Daniele and Santamaria-Martínez, Albert and Walsh, Logan A. and Ciriello, Giovanni},
  journal = {Nature Genetics},
  year    = {2024},
  volume  = {56},
  pages   = {74--84},
  doi     = {10.1038/s41588-023-01588-4}
}
```

## Additional citations for the bundled CODEX example

Analyses using the bundled mouse-spleen teaching data or the full reproduction
should also cite:

- Goltsev, Y. et al. Deep profiling of mouse splenic architecture with CODEX
  multiplexed imaging. *Cell* **174**, 968–981.e15 (2018).
  https://doi.org/10.1016/j.cell.2018.07.010
- Original dataset: https://doi.org/10.17632/zjnpwh8m5b.1
- Clustered-data record: https://doi.org/10.6084/m9.figshare.25040192.v2

The clustered-data record is licensed under CC BY 4.0. Consult
`inst/extdata/DATA_LICENSES.md` for the package's data-attribution details.
