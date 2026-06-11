# SIF to igraph

This function imports and converts a SIF file from your local hard
drive, cloud space, or internet into a graph with an igraph class, which
can then be used for the identification of most influential nodes via
the ivi function, for instance.

## Usage

``` r
sif2igraph(Path, directed = FALSE)
```

## Arguments

- Path:

  A string or character vector indicating the path to the desired SIF
  file. The SIF file could be on your local hard drive, cloud space, or
  on the internet.

- directed:

  Logical scalar, whether or not to create a directed graph.

## Value

An igraph graph object.

## Examples

``` r
if (FALSE) { # \dontrun{
MyGraph <- sif2igraph(Path = "/Users/User1/Desktop/mygraph.sif", directed=FALSE)
} # }
```
