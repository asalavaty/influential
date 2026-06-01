# ClusterRank (CR)

This function calculates the ClusterRank of input vertices and works
with both directed and undirected networks. This function and all of its
descriptions have been adapted from the centiserve package with some
minor modifications. ClusterRank is a local ranking algorithm which
takes into account not only the number of neighbors and the neighbors’
influences, but also the clustering coefficient.

## Usage

``` r
clusterRank(
  graph,
  vids = V(graph),
  directed = FALSE,
  loops = TRUE,
  ncores = "default",
  verbose = FALSE
)
```

## Arguments

- graph:

  The input graph as igraph object

- vids:

  Vertex sequence, the vertices for which the centrality values are
  returned. Default is all vertices.

- directed:

  Logical scalar, whether to directed graph is analyzed. This argument
  is ignored for undirected graphs.

- loops:

  Logical; whether the loop edges are also counted.

- ncores:

  Integer; the number of cores to be used for parallel processing. If
  ncores == "default" (default), the number of cores to be used will be
  the max(number of available cores) - 1. We recommend leaving ncores
  argument as is (ncores = "default").

- verbose:

  Logical; whether the accomplishment of different stages of the
  algorithm should be printed (default is FALSE).

## Value

A numeric vector contaning the ClusterRank centrality scores for the
selected vertices.

## See also

[`ivi`](https://asalavaty.github.io/influential/reference/ivi.md),
[`cent_network.vis`](https://asalavaty.github.io/influential/reference/cent_network.vis.md)

Other centrality functions:
[`betweenness()`](https://asalavaty.github.io/influential/reference/betweenness.md),
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
cr <- clusterRank(graph = My_graph, vids = GraphVertices, 
directed = FALSE, loops = TRUE, ncores = 1)
} # }
```
