Adrian (Abbas) Salavaty
16/01/2020 (updated on 17 January, 2020)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# influential

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/asalavaty/influential?branch=master&svg=true)](https://ci.appveyor.com/project/asalavaty/influential)
<!-- badges: end -->

The goal of influential is to help identification of the most
influential nodes (hubs) in a network. This package contains functions
for reconstruction of networks from adjacency matrices and data frames,
analysis of the topology of the network and calculation of centrality
measures, and identification of the most influential nodes (network
hubs). Also, some functions have been provided for the assessment of
dependence and correlation of two network centrality measures as well as
the conditional probability of deviation from their corresponding means
in opposite directions.

## Installation

You can install the released version of influential from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("influential")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(influential)
MyData <- centrality.measures
My.vertices.IHS <- ihs(DC = centrality.measures$Degree,
                       BC = centrality.measures$BetweennessCentrality,
                       NC = centrality.measures$NeighborhoodConnectivity)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
