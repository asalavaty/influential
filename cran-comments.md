CRAN Comments
================
Adrian (Abbas) Salavaty
15/01/2020 (updated on 28 May, 2026)

## 2.3.1 New version submission

This is a new submission of version 2.3.1 of the package `influential`.

This release substantially improves the `exir` and `ivi` workflows, with
a focus on performance, scalability, and support for modern omics data
structures. The `exir` function now supports data frames, tibbles,
matrices, sparse matrices, and Seurat objects as experimental input. It
also includes updated data preparation routines for bulk and single-cell
omics data, including optional TMM/logCPM normalization using `edgeR`,
pseudo-sampling/pseudo-bulking for large datasets, and conservative
feature filtering prior to association analysis.

The `exir` function has also been optimized across multiple
computational bottlenecks, including PCA, network construction,
neighbourhood scoring, driver/biomarker classification,
mediator-associated driver extraction, and verbose reporting. The `ivi`
workflow and several IVI-related helper functions were optimized while
preserving the original output of the corresponding centrality
calculations. Documentation and vignettes were updated accordingly.

Please see NEWS for complete details.

## Test environments

- local OS X install, R 4.5.0
- windows server (on AppVeyor CI), R 4.1.0
- ubuntu (r-hub CI), devel and release
- win-builder, devel and release

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

There was one NOTE from CRAN incoming checks.

The spell-check NOTE refers to package-specific and standard technical
terms used intentionally in the DESCRIPTION, including ExIR, IVI, omics,
and tibbles.

The Additional_repositories NOTE refers to the Bioconductor repository
used for the `edgeR` dependency. `edgeR` is used for TMM/logCPM
normalization of bulk count-like data and pseudo-bulked single-cell data
in ExIR. The package dependency check passed successfully.

## Downstream dependencies

This package has one reverse dependency: `networksem`.

I have run reverse dependency checks and confirm that there are no
breaking changes.
