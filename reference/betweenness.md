# Vertex betweenness centrality

This function and all of its descriptions have been obtained from the
igraph package.

## Usage

``` r
betweenness(
  graph,
  v = V(graph),
  directed = TRUE,
  weights = NULL,
  normalized = FALSE,
  ...
)
```

## Arguments

- graph:

  The graph to analyze (an igraph graph).

- v:

  The vertices for which the vertex betweenness will be calculated.

- directed:

  Logical, whether directed paths should be considered while determining
  the shortest paths.

- weights:

  Optional positive weight vector for calculating weighted betweenness.
  If the graph has a weight edge attribute, then this is used by
  default. Weights are used to calculate weighted shortest paths, so
  they are interpreted as distances.

- normalized:

  Logical scalar, whether to normalize the betweenness scores. If TRUE,
  then the results are normalized.

- ...:

  Additional arguments according to the original
  [`betweenness`](https://r.igraph.org/reference/betweenness.html)
  function in the package igraph.

## Value

A numeric vector with the betweenness score for each vertex in v.

## See also

[`ivi`](https://asalavaty.github.io/influential/reference/ivi.md),
[`cent_network.vis`](https://asalavaty.github.io/influential/reference/cent_network.vis.md),
and [`betweenness`](https://r.igraph.org/reference/betweenness.html) for
a complete description on this function

Other centrality functions:
[`clusterRank()`](https://asalavaty.github.io/influential/reference/clusterRank.md),
[`collective.influence()`](https://asalavaty.github.io/influential/reference/collective.influence.md),
[`h_index()`](https://asalavaty.github.io/influential/reference/h_index.md),
[`lh_index()`](https://asalavaty.github.io/influential/reference/lh_index.md),
[`neighborhood.connectivity()`](https://asalavaty.github.io/influential/reference/neighborhood.connectivity.md),
[`sirir()`](https://asalavaty.github.io/influential/reference/sirir.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- coexpression.data
My_graph <- graph_from_data_frame(MyData)
GraphVertices <- V(My_graph)
My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,
                                    directed = FALSE,
                                    normalized = FALSE)
} # }
```
