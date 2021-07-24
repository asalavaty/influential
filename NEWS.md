influential
================

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# influential 2.2.3.9000 (Developmental version)

-   The package website build with pkgdown is added to github.

-   The default correlation coefficient is changed from 0.3 to 0.5 in
    the ExIR function and its shiny app.

# influential 2.2.3 (CRAN version)

-   The ExIR shiny app is debugged and updated.

-   The ExIR function is debugged and updated.

# influential 2.2.2 (CRAN version)

-   The DESCRIPTION file is updated.

-   The function `ExIR` as well as the ExIR shiny app are debugged.

-   Documentations of some functions are updated.

# influential 2.2.1 (CRAN version)

-   The DESCRIPTION file is updated.

-   The name of the function for running shiny apps is changed to
    `runShinyApp`.

-   Shiny apps are updated in the package for local use.

-   Documentations of several functions are updated.

-   The Vignettes are updated.

-   The Read Me file is updated.

# influential 2.2.0 (CRAN version)

-   Update package Vignettes.

-   Documentation on how to access ExIR shiny app is added to the Read
    Me file, vignettes, and the corresponding functions’ documentations.

-   The ExIR shiny app is added to the package.

-   A Shiny app is developed for the running the ExIR model and
    visualization of its output.

-   Add the `comp_manipulate` function for the simulation of gene
    knockout and up-regulation.

-   Debug the function IVI.

-   Debug the function ExIR.

-   Add dependence to R package janitor for the correction of illegal
    characters in feature names.

-   Documentation on how to access IVI shiny app is added to the Read Me
    file, vignettes, and the corresponding functions’ documentations.

-   The IVI shiny app is added to the package.

-   A Shiny app is developed for the calculation of IVI as well as
    IVI-based network visualization.

-   The dependence of the function `ExIR` on the package reshape2 is
    removed.

-   The function `ExIR` is debugged.

# influential 2.0.1 (CRAN version)

-   The function `ExIR` is debugged.

# influential 2.0.0 (CRAN version)

-   The ReadMe file is updated.

-   The citation details of the package is updated.

-   Vignettes of the package are updated and extended.

-   A new function named `exir.vis` is added for the visualization of
    the results of the function `exir.`

-   A new function named `cent_network.vis` is added for the
    visualization of a network based on a centrality measure.

-   A new function named `graph_from_incidence_matrix` is added
    (imported from the *igraph* package) for the network reconstruction
    from an incidence matrix.

-   A new column named **Type** is added to the *Biomarker* table of the
    function `ExIR`.

-   The Z-score and statistical significance are added to the results of
    the function `ExIR`.

-   The function `ExIR` is updated so that the experimental data are
    ranked prior to correlation analysis. This will result in the
    assessment of the association based on the Spearman method, a more
    robust algorithm in variable conditions and/or non-parametric
    distributions.

-   The version of Roxygen2 is updated in the DESCRIPTION file.

-   The documentation of the `ExIR` function is updated.

-   The `ExIR` function is updated to prevent outputting NULL results.

-   The package logo is updated in accordance with the logo design of
    well-known R packages.

-   Some documentations in the Read Me file are corrected.

# influential 1.1.2 (CRAN version)

-   The `NC` function is improved.

-   The documentation of the `exir` function is updated.

# influential 1.1.1 (CRAN version)

-   The DESCRIPTION file is updated according to all of the
    modifications.

-   The citation details of the package is updated.

-   The function `diff.data.assembly` is added for assembling the
    differential/regression data required as an input for the `exir`
    function.

-   The function `exir` is added for the experimental-data-based
    classification and ranking of top candidate features.

-   The function `clusterrank` is upgraded to `clusterRank` in order to
    calculate the ClusterRank centrality internally.

-   The function `neighborhood.connectivity` is undergone minor
    modifications to enhance the performance speed.

-   The function `sif2igraph` is added for automatic importing and
    conversion of a SIF file from the users’ local hard drive, cloud
    space, or internet into a graph with an igraph class.

-   The thesis advisors are added to the DESCRIPTION file.

-   Some minor corrections are done in the documentation of functions
    and DESCRIPTION file.

# influential 1.0.0 (CRAN version)

-   The DESCRIPTION file is updated.

-   The Vignettes file is updated.

-   The function `sirir` is added for the unsupervised influence ranking
    of network nodes.

-   The function `hubness.score` is added for the calculation of the
    Hubness score.

-   The function `spreading.score` is added for the calculation of the
    Spreading score.

-   The function `ivi` is added for the influential node identification
    from a graph.

-   The function `ivi.from.indices` is added for the influential node
    identification from indices.

-   The function `ihs` is removed as a refined formula is added for
    influential node identification.

-   Updating the normality assessment results of association functions
    for vectors of length &lt; 4

-   The function `clusterrank` is added for the calculation of
    ClusterRank.

-   The function `collective.influence` is added for the calculation of
    Collective Influence.

-   The function `lh_index` is added for the calculation of local
    H-index.

-   The function `h_index` is added for the calculation of H-index.

-   The function `neighborhood.connectivity` is updated so that it will
    not return any NA or NaN values and the output will be a numeric
    vector.

-   The formatting of return values in the documentation of each
    function is corrected.

-   The `citation` details of the package are updated according to its
    associated published paper.

-   The error regarding the use of `NNS` since the second use of
    association assessment functions including `cond.prob.analysis` and
    `double.cent.assess.noRegression` is now corrected.
