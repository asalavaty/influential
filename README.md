
<!-- README.md is generated from README.Rmd. Please edit that file -->

# influential

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/asalavaty/influential?branch=master&svg=true)](https://ci.appveyor.com/project/asalavaty/influential)
[![](https://www.r-pkg.org/badges/version/influential?color=blue)](https://cran.r-project.org/package=influential)
[![](http://cranlogs.r-pkg.org/badges/grand-total/influential?color=green)](https://cran.r-project.org/package=influential)
[![](https://img.shields.io/badge/First%20integrative%20method%20for-Hub%20identification-blue.svg)](https://dx.doi.org/10.2139/ssrn.3565980)
<!-- badges: end -->

![The influential R package
logo](https://github.com/asalavaty/influential/blob/master/logo.png)

## Overview

The goal of `influential` is to help identification of the most
influential nodes in a network. This package contains functions for
reconstruction of networks from adjacency matrices and data frames,
analysis of the topology of the network and calculation of centrality
measures as well as a novel and powerful influential node ranking. The
first integrative method, namely the **Integrated Vector of Influence
(IVI)**, that captures all topological dimensions of the network for the
identification of network most influential nodes is also provided as a
function. Also, neighborhood connectivity, H-index, local H-index, and
collective influence (CI), all of which required centrality measures for
the calculation of **IVI**, are for the first time provided in an R
package. Furthermore, some functions have been provided for the
assessment of dependence and correlation of two network centrality
measures as well as the conditional probability of deviation from their
corresponding means in opposite directions.

Check out [**our paper**](https://doi.org/10.1101/2020.02.17.953430) for
a more complete description of the IVI formula and all of its
underpinning methods and analyses.

## Author

The `influential` package was written by [Adrian (Abbas)
Salavaty](https://www.AbbasSalavaty.com)

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

[Detailed description of the functions and their
outputs](https://github.com/asalavaty/influential/blob/master/vignettes/Vignettes.md)

You may browse Vignettes from within R using the following code.

``` r
browseVignettes("influential")
```

## An Example for Calculation of IVI

This is a basic example which shows you how to solve a common problem:

``` r
library(influential)

MyData <- centrality.measures # A data frame of centrality measures

# This function calculates the Integrated Vector of Influence (IVI)
My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,       # Calculation of IVI
                                   CR = centrality.measures$CR,
                                   NC = centrality.measures$NC,
                                   LH_index = centrality.measures$LH_index,
                                   BC = centrality.measures$BC,
                                   CI = centrality.measures$CI)

print(head(My.vertices.IVI))
#> [1] 24.670056  8.344337 18.621049  1.017768 29.437028 33.512598
```

## How to cite `influential`

To cite `influential`, please cite the [**associated
paper**](https://dx.doi.org/10.2139/ssrn.3565980). You can also refer to
the package’s citation information using the citation() function.

``` r
citation("influential")
```

## How to contribute

Please don’t hesitate to report any bugs/issues and request for
enhancement or any other contributions. To submit a bug report or
enhancement request, please use the [`influential` GitHub issues
tracker](https://github.com/asalavaty/influential/issues).
