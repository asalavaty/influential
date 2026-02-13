influential
================

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# influential 2.3.0 (CRAN version)

- Added first-order and second-order associated drivers of mediators to
  the final mediator result tables.

- Implemented matrix-based linear algebra formulation and C++ code in
  the `fcor` function and consequently the association analysis module
  of the function `exir`, resulting in a highly optimized and
  significantly faster association analysis.

- Removed the data.table dependency, as it is no longer required.

- Enabled parallel multi-core processing across multiple components of
  the function `exir`, including the supervised machine learning module.

- The function `exir` is optimized. It now automatically handles NA
  values of the input experimental data and convert them to zero.

- Debug the function `cent_network.vis`.

# influential 2.2.9 (CRAN version)

- Add the reexports function.

- Update the documentation of the modified functions.

- The syntax of the functions reexported from the igraph package are
  corrected.

- The package version is updated.

# influential 2.2.8 (CRAN version)

- The function `comp_manipulate` is optimized.

- The function `exir` is optimized. Also, the local and online version
  of the ExIR shiny apps are updated according to this optimization
  (calculation of zscores from the raw data rather than range normalized
  data).

- The argument ‘scaled’ of the functions `ivi`, `spreading.score`, and
  `hubness_score` is changed to ‘scale’ and its functionality is
  optimized. Also, the performance and documentations of these functions
  are updated accordingly. Additionally, the local and online version of
  the IVI shiny apps are updated according to these changes.

- SIRIR function, `sirir`, is optimized by providing access to several
  cores. Also, its documentation is updated.

- The function `collective.influence` is updated to match the current
  updates of the `igraph` package.

- The package `readr` is added to the ‘Suggests’ section of the
  Description file.

- The function `runShinyApp` is updated to require the `readr` package.

- The IVI local shiny app is updated to use the `readr` package for
  loading the datasets.

- The IVI local shiny app is debugged.

# influential 2.2.7 (CRAN version)

- The README file is updated.

- DESCRIPTION of the package is updated.

- New packages including *foreach*, and *doParallel* are added to the
  Imports section of the DESCRIPTION.

- ExIR function, `exir`, is optimized by providing access to several
  cores. Also, its documentation is updated.

- IVI function, `IVI`, is optimized by providing access to several
  cores. Also, its documentation is updated.

- LH-index function, `clusterRank`, is optimized by providing access to
  several cores. This also significantly speeds up the `ivi` function.
  Also, its documentation is updated.

- LH-index function, `lh_index`, is optimized by providing access to
  several cores. This also significantly speeds up the `ivi` function.
  Also, its documentation is updated.

- A `verbose` argument is added all centrality measure functions and
  their corresponding documentations are updated as well.

- The fcor function, `fcor`, as well as its documentation are optimized
  and updated.

- The documentation of the SIRIR function, `sirir`, is updated.

- Hubness score function, `hubness.score`, is debugged and optimized.

- Spreading score function and its documentation, `spreading.score`, is
  debugged and optimized.

- IVI function, `ivi`, is debugged and optimized.

- IVI from indices function, `ivi.from.indices`, is debugged and
  optimized.

- Collective Influence function, `collective.influence`, is optimized.

- The documentation of the Collective Influence function,
  `collective.influence`, is updated.

# influential 2.2.6 (CRAN version)

- ExIR Shiny app is updated according to the ExIR function.

- ExIR visualization function, `exir.vis`, is optimized.

- The package vignettes are updated and the `fcor` function described
  and exemplified in the vignettes.

- The ExIR function is optimized (corrected for proper Spearman
  correlation analysis and added Mutual rank as a measure for filtering
  top correlations) and speeded up (using the data.table package).

- A function named `fcor` is added to the package for super-fast
  correlation analysis of large datasets and simultaneous P-value and
  Mutual Rank calculations.

# influential 2.2.5 (CRAN version)

- The package vignettes are updated.

- The `betweenness` function is corrected.

- Documentation of some functions corrected.

- Links to my personal website are corrected.

# influential 2.2.4 (CRAN version)

- The package vignettes are updated.

- The implementations of shiny apps are optimized.

- Shiny apps are debugged.

- A new argument, `label.position`, is added to the `exir.vis` function.

- The `exir.vis` function is debugged and optimized.

- The upload file size limit is removed from the local shiny apps.

- The package website built with pkgdown is added to github.

- The default correlation coefficient is changed from 0.3 to 0.5 in the
  ExIR function and its shiny app.

# influential 2.2.3 (CRAN version)

- The ExIR shiny app is debugged and updated.

- The ExIR function is debugged and updated.

# influential 2.2.2 (CRAN version)

- The DESCRIPTION file is updated.

- The function `ExIR` as well as the ExIR shiny app are debugged.

- Documentations of some functions are updated.

# influential 2.2.1 (CRAN version)

- The DESCRIPTION file is updated.

- The name of the function for running shiny apps is changed to
  `runShinyApp`.

- Shiny apps are updated in the package for local use.

- Documentations of several functions are updated.

- The Vignettes are updated.

- The Read Me file is updated.

# influential 2.2.0 (CRAN version)

- Update package Vignettes.

- Documentation on how to access ExIR shiny app is added to the Read Me
  file, vignettes, and the corresponding functions’ documentations.

- The ExIR shiny app is added to the package.

- A Shiny app is developed for the running the ExIR model and
  visualization of its output.

- Add the `comp_manipulate` function for the simulation of gene knockout
  and up-regulation.

- Debug the function IVI.

- Debug the function ExIR.

- Add dependence to R package janitor for the correction of illegal
  characters in feature names.

- Documentation on how to access IVI shiny app is added to the Read Me
  file, vignettes, and the corresponding functions’ documentations.

- The IVI shiny app is added to the package.

- A Shiny app is developed for the calculation of IVI as well as
  IVI-based network visualization.

- The dependence of the function `ExIR` on the package reshape2 is
  removed.

- The function `ExIR` is debugged.

# influential 2.0.1 (CRAN version)

- The function `ExIR` is debugged.

# influential 2.0.0 (CRAN version)

- The ReadMe file is updated.

- The citation details of the package is updated.

- Vignettes of the package are updated and extended.

- A new function named `exir.vis` is added for the visualization of the
  results of the function `exir.`

- A new function named `cent_network.vis` is added for the visualization
  of a network based on a centrality measure.

- A new function named `graph_from_incidence_matrix` is added (imported
  from the *igraph* package) for the network reconstruction from an
  incidence matrix.

- A new column named **Type** is added to the *Biomarker* table of the
  function `ExIR`.

- The Z-score and statistical significance are added to the results of
  the function `ExIR`.

- The function `ExIR` is updated so that the experimental data are
  ranked prior to correlation analysis. This will result in the
  assessment of the association based on the Spearman method, a more
  robust algorithm in variable conditions and/or non-parametric
  distributions.

- The version of Roxygen2 is updated in the DESCRIPTION file.

- The documentation of the `ExIR` function is updated.

- The `ExIR` function is updated to prevent outputting NULL results.

- The package logo is updated in accordance with the logo design of
  well-known R packages.

- Some documentations in the Read Me file are corrected.

# influential 1.1.2 (CRAN version)

- The `NC` function is improved.

- The documentation of the `exir` function is updated.

# influential 1.1.1 (CRAN version)

- The DESCRIPTION file is updated according to all of the modifications.

- The citation details of the package is updated.

- The function `diff.data.assembly` is added for assembling the
  differential/regression data required as an input for the `exir`
  function.

- The function `exir` is added for the experimental-data-based
  classification and ranking of top candidate features.

- The function `clusterrank` is upgraded to `clusterRank` in order to
  calculate the ClusterRank centrality internally.

- The function `neighborhood.connectivity` is undergone minor
  modifications to enhance the performance speed.

- The function `sif2igraph` is added for automatic importing and
  conversion of a SIF file from the users’ local hard drive, cloud
  space, or internet into a graph with an igraph class.

- The thesis advisors are added to the DESCRIPTION file.

- Some minor corrections are done in the documentation of functions and
  DESCRIPTION file.

# influential 1.0.0 (CRAN version)

- The DESCRIPTION file is updated.

- The Vignettes file is updated.

- The function `sirir` is added for the unsupervised influence ranking
  of network nodes.

- The function `hubness.score` is added for the calculation of the
  Hubness score.

- The function `spreading.score` is added for the calculation of the
  Spreading score.

- The function `ivi` is added for the influential node identification
  from a graph.

- The function `ivi.from.indices` is added for the influential node
  identification from indices.

- The function `ihs` is removed as a refined formula is added for
  influential node identification.

- Updating the normality assessment results of association functions for
  vectors of length \< 4

- The function `clusterrank` is added for the calculation of
  ClusterRank.

- The function `collective.influence` is added for the calculation of
  Collective Influence.

- The function `lh_index` is added for the calculation of local H-index.

- The function `h_index` is added for the calculation of H-index.

- The function `neighborhood.connectivity` is updated so that it will
  not return any NA or NaN values and the output will be a numeric
  vector.

- The formatting of return values in the documentation of each function
  is corrected.

- The `citation` details of the package are updated according to its
  associated published paper.

- The error regarding the use of `NNS` since the second use of
  association assessment functions including `cond.prob.analysis` and
  `double.cent.assess.noRegression` is now corrected.
