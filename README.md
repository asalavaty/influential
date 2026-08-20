
<!-- README.md is generated from README.Rmd. Please edit that file -->

# influential <a href="https://asalavaty.com/widgets/influential"><img src="man/figures/Symbol.png" align="right" height="221" alt="influential logo"/></a>

### Network influence analysis, centrality assessment, and experimental feature prioritization in R

<!-- badges: start -->

[![R-CMD-check](https://github.com/asalavaty/influential/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/asalavaty/influential/actions/workflows/R-CMD-check.yaml)
[![CRAN
Version](https://www.r-pkg.org/badges/version/influential?color=blue)](https://cran.r-project.org/package=influential)
[![Total CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/influential?color=8B0000)](https://cran.r-project.org/package=influential)
[![Monthly CRAN
Downloads](https://cranlogs.r-pkg.org/badges/influential?color=8B0000)](https://cran.r-project.org/package=influential)
[![Integrated Value of Influence:
IVI](https://img.shields.io/badge/Integrated%20Value%20of%20Influence-IVI-006400.svg?labelColor=f5f5f5)](https://doi.org/10.1016/j.patter.2020.100052)
[![SIR-based Influence Ranking:
SIRIR](https://img.shields.io/badge/SIR--based%20Influence%20Ranking-SIRIR-006400.svg?labelColor=f5f5f5)](https://doi.org/10.1016/j.patter.2020.100052)
[![Experimental data-based Integrative Ranking:
ExIR](https://img.shields.io/badge/Experimental%20data--based%20Integrative%20Ranking-ExIR-006400.svg?labelColor=f5f5f5)](https://doi.org/10.1016/j.isci.2026.116303)

<!-- badges: end -->

## Overview

`influential` is an R package for identifying influential nodes in
networks and for classifying and prioritizing candidate features from
experimental data. It brings together **association analysis, network
reconstruction, centrality assessment, influence ranking, visualization,
simulation, and experimental feature prioritization** in a single
toolkit.

Two complementary workflows sit at the center of the package:

- **Integrated Value of Influence (IVI)** integrates local, semi-local,
  and global topological information to identify influential network
  nodes while reducing the limitations of relying on any single
  centrality measure.
- **Experimental data-based Integrative Ranking (ExIR)** combines
  experimental evidence, machine learning, network reconstruction, and
  influence ranking to classify and prioritize candidate **drivers,
  biomarkers, and mediators** from omics data.

The package also provides **Hubness score** for local network power,
**Spreading score** for information-spreading potential, **SIRIR** for
simulation-based influence ranking, computational
knockout/up-regulation, fast correlation analysis, network
reconstruction utilities, centrality measures, and centrality-based
network visualization.

### At a glance

| Capability | What `influential` provides |
|:---|:---|
| **Association analysis** | Fast Pearson/Spearman correlation analysis with optional mutual rank, p-values, and adjusted p-values |
| **Network reconstruction** | Construction of `igraph` networks from data frames, adjacency matrices, incidence matrices, and SIF files |
| **Centrality analysis** | Local, semi-local, and global centrality measures together with association assessment |
| **Network influence** | IVI, Hubness score, Spreading score, and SIRIR |
| **Experimental prioritization** | ExIR-based classification and ranking of drivers, biomarkers, and mediators from omics data |
| **Perturbation & visualization** | In silico knockout/up-regulation and centrality-based network visualization |
| **Interactive analysis** | Browser-based and locally launchable Shiny interfaces for IVI and ExIR |

## Quick links

| Resource | Link |
|:---|:---|
| ✨ **Feature Explorer** | [Explore the features and capabilities of `influential`](https://asalavaty.com/widgets/influential) |
| 📦 **CRAN** | [cran.r-project.org/package=influential](https://cran.r-project.org/package=influential) |
| 📖 **Full vignette** | [Introduction to influential](https://asalavaty.github.io/influential/articles/Vignettes.html) |
| 💻 **GitHub** | [github.com/asalavaty/influential](https://github.com/asalavaty/influential) |
| 🌐 **Interactive web portal** | [Influential Software Package](https://influential.erc.monash.edu/) |
| 🐞 **Issues & feature requests** | [GitHub issue tracker](https://github.com/asalavaty/influential/issues) |

## Installation

Install the current CRAN release:

``` r
install.packages("influential")
```

Or install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github(
  "asalavaty/influential",
  build_vignettes = TRUE
)
```

Then load the package:

``` r
library(influential)
```

## Core methods

### Integrated Value of Influence (IVI)

IVI integrates complementary **local, semi-local, and global**
centrality dimensions to identify influential nodes within a network.

For the methodological details, see:

[**Integrated Value of Influence: An Integrative Method for the
Identification of the Most Influential Nodes within
Networks**](https://doi.org/10.1016/j.patter.2020.100052)

### Experimental data-based Integrative Ranking (ExIR)

ExIR prioritizes candidate features directly from experimental omics
data by integrating multiple levels of evidence with network
reconstruction and influence ranking. Depending on the input data and
analysis settings, ExIR can identify **drivers, biomarkers,
DE-mediators, and nonDE-mediators**.

For the methodological details, see:

[**ExIR enables prioritizing driver and biomarker genes from omics data
in a reference free
manner**](https://doi.org/10.1016/j.isci.2026.116303)

## Documentation

A comprehensive introduction to `influential` and its functions is
available in the package vignette:

**[Read the influential
vignette](https://asalavaty.github.io/influential/articles/Vignettes.html)**

You can also browse installed vignettes directly from R:

``` r
browseVignettes("influential")
```

## Shiny apps

The package provides interactive interfaces for IVI and ExIR through the
[Influential Software Package web
portal](https://influential.erc.monash.edu/).

### IVI Shiny app

The [IVI Shiny App](https://influential.erc.monash.edu/IVI/) supports
calculation of IVI values and IVI-based network visualization.

You can also launch it locally:

``` r
influential::runShinyApp("IVI")
```

### ExIR Shiny app

The [ExIR Shiny App](https://influential.erc.monash.edu/ExIR/) supports
ExIR analysis, result visualization, and downstream exploration.

You can also launch it locally:

``` r
influential::runShinyApp("ExIR")
```

## How to cite `influential`

If you use `influential`, please cite the publication associated with
the method(s) used in your analysis.

For the **Experimental data-based Integrative Ranking (ExIR)** model:

- Salavaty A, Douek AM, Kaslin J, Ramialison M, Currie PD. ExIR enables
  prioritizing driver and biomarker genes from omics data in a reference
  free manner. *iScience*. 2026.06.19. [Read
  online](https://doi.org/10.1016/j.isci.2026.116303).

For the **Integrated Value of Influence (IVI)** and network influence
analysis:

- Salavaty A, Ramialison M, Currie PD. Integrated Value of Influence: An
  Integrative Method for the Identification of the Most Influential
  Nodes within Networks. *Patterns*. 2020.08.14. [Read
  online](https://doi.org/10.1016/j.patter.2020.100052).

Package citation information is also available from R:

``` r
citation("influential")
```

## Author

The `influential` package was developed by [Adrian
Salavaty](https://asalavaty.com/).

### Advisors

- Mirana Ramialison
- Peter D. Currie

## Contributing and support

Bug reports, feature requests, documentation suggestions, and other
contributions are welcome.

Please use the [`influential` GitHub issues
tracker](https://github.com/asalavaty/influential/issues) to report
problems or suggest enhancements.
