Vignettes of the `influential`
================

## Overview

The influential package contains several functions that could be
categorized into four groups according to their purpose:

  - Network reconstruction
  - Calculation of centrality measures
  - Assessment of the association of centrality measures
  - Identification of the most `influential` network nodes (hubs)

-----

``` r
library(influential)
```

## Network reconstruction

Two functions have been obtained from the `igraph` package for the
reconstruction of networks.

  - ### From a data frame
    
    In the data frame the first and second columns should be composed of
    source and target nodes.
    
    A sample appropriate data frame is brought below:

| lncRNA      | Coexpressed.Gene |
| :---------- | :--------------- |
| ADAMTS9-AS2 | A2M              |
| ADAMTS9-AS2 | ABCA6            |
| ADAMTS9-AS2 | ABCA8            |
| ADAMTS9-AS2 | ABCA9            |
| ADAMTS9-AS2 | ABI3BP           |
| ADAMTS9-AS2 | AC093110.3       |

This is a co-expression dataset obtained from
[PMID: 31211495](https://www.ncbi.nlm.nih.gov/pubmed/31211495).

``` r
MyData <- coexpression.data
My_graph <- graph_from_data_frame(d=MyData)
```

If you look at the class of `My_graph` you should see that it has an
`igraph` class:

``` r
class(My_graph)
#> [1] "igraph"
```

-----

  - ### From an adjacency matrix

A sample appropriate data frame is brought
below:

|             | ADAMTS9-AS2 | C8orf34-AS1 | CADM3-AS1 | FAM83A-AS1 | FENDRR | LANCL1-AS1 | LINC00092 | LINC00467 | LINC00857 | LINC00891 |
| ----------- | ----------: | ----------: | --------: | ---------: | -----: | ---------: | --------: | --------: | --------: | --------: |
| ADAMTS9-AS2 |           0 |           0 |         0 |          0 |      0 |          1 |         0 |         0 |         0 |         1 |
| C8orf34-AS1 |           0 |           0 |         0 |          0 |      0 |          0 |         0 |         0 |         0 |         1 |
| CADM3-AS1   |           0 |           0 |         0 |          0 |      1 |          0 |         0 |         0 |         0 |         0 |
| FAM83A-AS1  |           0 |           0 |         0 |          0 |      0 |          0 |         0 |         0 |         1 |         0 |
| FENDRR      |           0 |           0 |         0 |          0 |      0 |          0 |         0 |         0 |         0 |         0 |
| LANCL1-AS1  |           0 |           0 |         0 |          0 |      1 |          0 |         0 |         0 |         0 |         0 |

``` r
MyData <- coexpression.adjacency
My_graph <- graph_from_adjacency_matrix(MyData)
```

-----

## Calculation of centrality measures

  - ### Degree centrality

  - ### Betweenness centrality

  - ### Neighborhood connectivity

-----

## Assessment of the association of centrality measures

## Identification of the most `influential` network nodes (hubs)
