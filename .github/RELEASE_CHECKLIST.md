# CellCharterR external release checklist

The source tree is ready for these maintainer-owned external actions. Do not
mark them complete until the public endpoint itself has been checked from an
unauthenticated session.

## 1. Publish the repository and HTML

1. Create the public repository `cperalta-viteri/CellCharterR`.
2. Push the complete reviewed tree to its default branch. The current checkout
   has no Git remote, so this step cannot be performed implicitly.
3. In repository settings, select **GitHub Actions** as the Pages source.
4. Run `Publish vignette HTML` or push a change under `docs/`.
5. Verify HTTP 200 responses for:
   - `https://cperalta-viteri.github.io/CellCharterR/`
   - `https://cperalta-viteri.github.io/CellCharterR/CellCharterR-beginner-workflow.html`
   - `https://cperalta-viteri.github.io/CellCharterR/CellCharterR-h5ad-to-sce.html`
   - `https://cperalta-viteri.github.io/CellCharterR/CellCharterR-mouse-spleen-reproducibility.html`

At the 2026-10-12 checkpoint, the repository and Pages endpoints both returned
HTTP 404. The workflow and rendered files exist locally; publication has not
been claimed.

## 2. Publish the portable full SCE

Upload
`results/CellCharterR_mouse_spleen_v2/portable/CellCharterR_mouse_spleen_sce.rds`
to a durable, versioned public data record such as Figshare or Zenodo. Record:

- filename: `CellCharterR_mouse_spleen_sce.rds`;
- compressed size: 286.0 MiB;
- loaded R object size: approximately 693.7 MiB;
- cells: 707,474;
- SHA-256:
  `d9adda271ad26aca5887a2d4f38feea8bb3a1503b561e34ee9f13df129371219`;
- data licenses and citations from `inst/extdata/DATA_LICENSES.md`; and
- the repository commit and package version that produced the record.

After upload, insert the immutable DOI/download URL into the H5AD-to-SCE
vignette and `BIOCONDUCTOR.md`. Download it into a new temporary directory,
verify SHA-256, and repeat the read-only checks documented in
`PORTABILITY_VALIDATION.md`. Do not rerun the full analysis.

## 3. Complete Bioconductor identity gates

1. Register `christian_94_97@hotmail.com` on the Bioconductor Support Site.
2. From a clean clone of the published repository, run:

       Rscript scripts/build_package_safely.R
       R CMD check --no-manual CellCharterR_0.99.0.tar.gz
       Rscript -e 'BiocCheck::BiocCheck("CellCharterR_0.99.0.tar.gz", new.package = TRUE)'
       Rscript -e 'BiocCheck::BiocCheckGitClone(".")'

   The safe builder stages only standard package inputs before invoking
   `R CMD build`. Use it instead of `R CMD build .` in a working tree that also
   contains ignored large data, because R copies the source tree into `Rtmp`
   before applying `.Rbuildignore`.

3. Confirm the Linux CPU GitHub Actions check and Pages deployment are green.
4. Confirm no unexplained ERROR or WARNING remains before submission.
