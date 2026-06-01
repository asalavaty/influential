# H-index

This function calculates the H-index of input vertices and works with
both directed and undirected networks.

## Usage

``` r
h_index(graph, vertices = V(graph), mode = "all", verbose = FALSE)
```

## Arguments

- graph:

  A graph (network) of the igraph class.

- vertices:

  A vector of desired vertices, which could be obtained by the V
  function.

- mode:

  The mode of H-index depending on the directedness of the graph. If the
  graph is undirected, the mode "all" should be specified. Otherwise,
  for the calculation of H-index based on incoming connections select
  "in" and for the outgoing connections select "out". Also, if all of
  the connections are desired, specify the "all" mode. Default mode is
  set to "all".

- verbose:

  Logical; whether the accomplishment of different stages of the
  algorithm should be printed (default is FALSE).

## Value

A vector including the H-index of each vertex inputted.

## See also

[`ivi`](https://asalavaty.github.io/influential/reference/ivi.md),
[`cent_network.vis`](https://asalavaty.github.io/influential/reference/cent_network.vis.md)

Other centrality functions:
[`betweenness()`](https://asalavaty.github.io/influential/reference/betweenness.md),
[`clusterRank()`](https://asalavaty.github.io/influential/reference/clusterRank.md),
[`collective.influence()`](https://asalavaty.github.io/influential/reference/collective.influence.md),
[`lh_index()`](https://asalavaty.github.io/influential/reference/lh_index.md),
[`neighborhood.connectivity()`](https://asalavaty.github.io/influential/reference/neighborhood.connectivity.md),
[`sirir()`](https://asalavaty.github.io/influential/reference/sirir.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- coexpression.data
My_graph <- graph_from_data_frame(MyData)
GraphVertices <- V(My_graph)
h.index <- h_index(graph = My_graph, vertices = GraphVertices, mode = "all")
} # }
```
