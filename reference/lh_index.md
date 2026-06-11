# local H-index (LH-index)

This function calculates the local H-index of input vertices and works
with both directed and undirected networks.

## Usage

``` r
lh_index(
  graph,
  vertices = V(graph),
  mode = "all",
  ncores = "default",
  verbose = FALSE
)
```

## Arguments

- graph:

  A graph (network) of the igraph class.

- vertices:

  A vector of desired vertices, which could be obtained by the V
  function.

- mode:

  The mode of local H-index depending on the directedness of the graph.
  If the graph is undirected, the mode "all" should be specified.
  Otherwise, for the calculation of local H-index based on incoming
  connections select "in" and for the outgoing connections select "out".
  Also, if all of the connections are desired, specify the "all" mode.
  Default mode is set to "all".

- ncores:

  Integer; the number of cores to be used for parallel processing. If
  ncores == "default" (default), the number of cores to be used will be
  the max(number of available cores) - 1. We recommend leaving ncores
  argument as is (ncores = "default").

- verbose:

  Logical; whether the accomplishment of different stages of the
  algorithm should be printed (default is FALSE).

## Value

A vector including the local H-index of each vertex inputted.

## See also

[`ivi`](https://asalavaty.github.io/influential/reference/ivi.md),
[`cent_network.vis`](https://asalavaty.github.io/influential/reference/cent_network.vis.md)

Other centrality functions:
[`betweenness()`](https://asalavaty.github.io/influential/reference/betweenness.md),
[`clusterRank()`](https://asalavaty.github.io/influential/reference/clusterRank.md),
[`collective.influence()`](https://asalavaty.github.io/influential/reference/collective.influence.md),
[`h_index()`](https://asalavaty.github.io/influential/reference/h_index.md),
[`neighborhood.connectivity()`](https://asalavaty.github.io/influential/reference/neighborhood.connectivity.md),
[`sirir()`](https://asalavaty.github.io/influential/reference/sirir.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- coexpression.data
My_graph <- graph_from_data_frame(MyData)
GraphVertices <- V(My_graph)
lh.index <- lh_index(graph = My_graph, vertices = GraphVertices, mode = "all", ncores = 1)
} # }
```
