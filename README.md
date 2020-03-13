
<!-- README.md is generated from README.Rmd. Please edit that file -->

# influential

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/asalavaty/influential?branch=master&svg=true)](https://ci.appveyor.com/project/asalavaty/influential)
[![](https://www.r-pkg.org/badges/version/influential?color=blue)](https://cran.r-project.org/package=influential)
[![](http://cranlogs.r-pkg.org/badges/grand-total/influential?color=green)](https://cran.r-project.org/package=influential)
[![](https://img.shields.io/badge/First%20integrative%20method%20for-Hub%20identification-blue.svg)](https://doi.org/10.1101/2020.02.17.953430)
<!-- badges: end -->

![The influential R package
logo](https://github.com/asalavaty/influential/blob/master/logo.png)

## Overview

The goal of `influential` is to help identification of the most
influential nodes (hubs) in a network. This package contains functions
for reconstruction of networks from adjacency matrices and data frames,
analysis of the topology of the network and calculation of centrality
measures. The first integrative method (i.e. the *integrated hubness
score (IHS)* method) for the identification of network hubs is also
provided as a function, which is the main purpose of the package. Also,
neighborhood connectivity, one of the required centrality measures for
the calculation of IHS, is for the first time calculable in the R
environment. Furthermore, some functions have been provided for the
assessment of dependence and correlation of two network centrality
measures as well as the conditional probability of deviation from their
corresponding means in opposite directions.

Check out [**our paper**](https://doi.org/10.1101/2020.02.17.953430) for
a more complete description of the IHS formula and all of its
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

## An Example for Calculation of IHS

This is a basic example which shows you how to solve a common problem:

``` r
library(influential)

MyData <- centrality.measures # A data frame of centrality measures

# This function calculates the Integrated Hubness Score (IHS)
My.vertices.IHS <- ihs(DC = centrality.measures$Degree,
                       BC = centrality.measures$BetweennessCentrality,
                       NC = centrality.measures$NeighborhoodConnectivity)

print(head(My.vertices.IHS))
#> [1] 196.2039215   2.9822273   0.1572078   6.1757221   0.3199993   0.5095222
```

## How to cite `influential`

To cite `influential`, please cite the [**associated
paper**](https://doi.org/10.1101/2020.02.17.953430). You can also refer
to the package’s citation information using the citation() function.

``` r
citation("influential")
```

## How to contribute

Please don’t hesitate to report any bugs/issues and request for
enhancement or any other contributions. To submit a bug report or
enhancement request, please use the [`influential` GitHub issues
tracker](https://github.com/asalavaty/influential/issues).
