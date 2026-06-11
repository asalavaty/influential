# Collective Influence (CI)

This function calculates the collective influence of input vertices and
works with both directed and undirected networks. This function and its
descriptions are obtained from
https://github.com/ronammar/collective_influence with minor
modifications. Collective Influence as described by Morone & Makse
(2015). In simple terms, it is the product of the reduced degree
(degree - 1) of a node and the total (sum of) reduced degrees of all
nodes at a distance d from the node.

## Usage

``` r
collective.influence(
  graph,
  vertices = V(graph),
  mode = "all",
  d = 3,
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

  The mode of collective influence depending on the directedness of the
  graph. If the graph is undirected, the mode "all" should be specified.
  Otherwise, for the calculation of collective influence based on
  incoming connections select "in" and for the outgoing connections
  select "out". Also, if all of the connections are desired, specify the
  "all" mode. Default mode is set to "all".

- d:

  The distance, expressed in number of steps from a given node
  (default=3). Distance must be \> 0. According to Morone & Makse
  (https://doi.org/10.1038/nature14604), optimal results can be reached
  at d=3,4, but this depends on the size/"radius" of the network. NOTE:
  the distance d is not inclusive. This means that nodes at a distance
  of 3 from our node-of-interest do not include nodes at distances 1
  and 2. Only 3.

- verbose:

  Logical; whether the accomplishment of different stages of the
  algorithm should be printed (default is FALSE).

## Value

A vector of collective influence for each vertex of the graph
corresponding to the order of vertices output by V(graph).

## See also

[`ivi`](https://asalavaty.github.io/influential/reference/ivi.md),
[`cent_network.vis`](https://asalavaty.github.io/influential/reference/cent_network.vis.md)

Other centrality functions:
[`betweenness()`](https://asalavaty.github.io/influential/reference/betweenness.md),
[`clusterRank()`](https://asalavaty.github.io/influential/reference/clusterRank.md),
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
ci <- collective.influence(graph = My_graph, vertices = GraphVertices, mode = "all", d=3)
} # }
```
