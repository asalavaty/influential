
<!-- README.md is generated from README.Rmd. Please edit that file -->

# influential <a href='https://github.com/asalavaty/influential'><img src='man/figures/logo.png' align="right" height="221" /></a>

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/asalavaty/influential?branch=master&svg=true)](https://ci.appveyor.com/project/asalavaty/influential)
[![](https://www.r-pkg.org/badges/version/influential?color=blue)](https://cran.r-project.org/package=influential)
[![](http://cranlogs.r-pkg.org/badges/grand-total/influential?color=green)](https://cran.r-project.org/package=influential)
[![Rdoc](http://www.rdocumentation.org/badges/version/influential)](http://www.rdocumentation.org/packages/influential)
[![](https://img.shields.io/badge/Integrated%20Value%20of%20Influence-IVI-blue.svg)](https://doi.org/10.1016/j.patter.2020.100052)
[![](https://img.shields.io/badge/SIR--based%20Influence%20Ranking-SIRIR-green.svg)](https://doi.org/10.1016/j.patter.2020.100052)
[![](https://img.shields.io/badge/Experimental%20data--based%20Integrative%20Ranking-ExIR-blue.svg)](https://github.com/asalavaty/influential)
<!-- badges: end -->

## Overview

The goal of `influential` is to help identification of the most
`influential` nodes in a network as well as the classification and
ranking of top candidate features. This package contains functions for
the classification and ranking of features, reconstruction of networks
from adjacency matrices and data frames, analysis of the topology of the
network and calculation of centrality measures as well as a novel and
powerful `influential` node ranking. The **Experimental data-based
Integrative Ranking (ExIR)** is a sophisticated model for classification
and ranking of the top candidate features based on only the experimental
data. The first integrative method, namely the **Integrated Value of
Influence (IVI)**, that captures all topological dimensions of the
network for the identification of network most `influential` nodes is
also provided as a function. Also, neighborhood connectivity, H-index,
local H-index, and collective influence (CI), all of which required
centrality measures for the calculation of **IVI**, are for the first
time provided in an R package. Additionally, a function is provided for
running **SIRIR** model, which is the combination of leave-one-out cross
validation technique and the conventional SIR model, on a network to
unsupervisedly rank the true influence of vertices. Furthermore, some
functions have been provided for the assessment of dependence and
correlation of two network centrality measures as well as the
conditional probability of deviation from their corresponding means in
opposite directions.

Check out [**our paper**](https://doi.org/10.1016/j.patter.2020.100052)
for a more complete description of the IVI formula and all of its
underpinning methods and analyses.

## Author

The `influential` package was written by [Adrian (Abbas)
Salavaty](https://www.AbbasSalavaty.com)

## Advisors

Mirana Ramialison and Peter D. Currie

## How to Install

You can install the official [CRAN
release](https://cran.r-project.org/package=influential) of the
`influential` with the following code:

``` r
install.packages("influential")
```

Or the development version from GitHub:

``` r
## install.packages("devtools")
devtools::install_github("asalavaty/influential")
```

## Vignettes

  - [Detailed description of functions and their
    outputs](https://cran.r-project.org/web/packages/influential/vignettes/Vignettes.html)

  - [Vignettes of the developmental version of the
    package](https://github.com/asalavaty/influential/blob/master/vignettes/Vignettes.md)

You may browse Vignettes from within R using the following code.

``` r
browseVignettes("influential")
```

## How to cite `influential`

To cite `influential`, please cite its associated paper:

  - Integrated Value of Influence: An Integrative Method for the
    Identification of the Most Influential Nodes within Networks. Abbas
    Salavaty, Mirana Ramialison, Peter D Currie. *Patterns*. 2020.08.14
    ([Read online](https://doi.org/10.1016/j.patter.2020.100052)).

You can also refer to the package’s citation information using the
`citation()` function.

``` r
citation("influential")
```

## How to contribute

Please don’t hesitate to report any bugs/issues and request for
enhancement or any other contributions. To submit a bug report or
enhancement request, please use the [`influential` GitHub issues
tracker](https://github.com/asalavaty/influential/issues).
