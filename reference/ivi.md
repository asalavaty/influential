# Integrated Value of Influence (IVI)

This function calculates the IVI of the desired nodes from a graph. \#'
A shiny app has also been developed for the calculation of IVI as well
as IVI-based network visualization, which is accessible using the
\`influential::runShinyApp("IVI")\` command. You can also access the
shiny app online at https://influential.erc.monash.edu/.

## Usage

``` r
ivi(
  graph,
  vertices = V(graph),
  weights = NULL,
  directed = FALSE,
  mode = "all",
  loops = TRUE,
  d = 3,
  scale = "range",
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

- weights:

  Optional positive weight vector for calculating weighted betweenness
  centrality of nodes as a requirement for calculation of IVI. If the
  graph has a weight edge attribute, then this is used by default.
  Weights are used to calculate weighted shortest paths, so they are
  interpreted as distances.

- directed:

  Logical scalar, whether to directed graph is analyzed. This argument
  is ignored for undirected graphs.

- mode:

  The mode of IVI depending on the directedness of the graph. If the
  graph is undirected, the mode "all" should be specified. Otherwise,
  for the calculation of IVI based on incoming connections select "in"
  and for the outgoing connections select "out". Also, if all of the
  connections are desired, specify the "all" mode. Default mode is set
  to "all".

- loops:

  Logical; whether the loop edges are also counted.

- d:

  The distance, expressed in number of steps from a given node
  (default=3). Distance must be \> 0. According to Morone & Makse
  (https://doi.org/10.1038/nature14604), optimal results can be reached
  at d=3,4, but this depends on the size/"radius" of the network. NOTE:
  the distance d is not inclusive. This means that nodes at a distance
  of 3 from our node-of-interest do not include nodes at distances 1
  and 2. Only 3.

- scale:

  Character string; the method used for scaling/normalizing the results.
  Options include 'range' (normalization within a 1-100 range),
  'z-scale' (standardization using the z-score), and 'none' (no data
  scaling). The default selection is 'range'. Opting for the 'range'
  method is suitable when exploring a single network, allowing you to
  observe the complete spectrum and distribution of node influences. In
  this case, there is no intention to establish a specific threshold for
  the outcomes. However, it is possible to identify and present the top
  influential nodes based on their rankings. Conversely, the 'z-scale'
  option proves advantageous if the aim is to compare node influences
  across multiple networks or if there is a desire to establish a
  threshold (usually z-score \> 1.645) for generating a list of the most
  influential nodes without manual intervention.

- ncores:

  Integer; the number of cores to be used for parallel processing. If
  ncores == "default" (default), the number of cores to be used will be
  the max(number of available cores) - 1. We recommend leaving ncores
  argument as is (ncores = "default").

- verbose:

  Logical; whether the accomplishment of different stages of the
  algorithm should be printed (default is FALSE).

## Value

A numeric vector with the IVI values based on the provided centrality
measures.

## See also

[`cent_network.vis`](https://asalavaty.github.io/influential/reference/cent_network.vis.md)

Other integrative ranking functions:
[`comp_manipulate()`](https://asalavaty.github.io/influential/reference/comp_manipulate.md),
[`exir()`](https://asalavaty.github.io/influential/reference/exir.md),
[`hubness.score()`](https://asalavaty.github.io/influential/reference/hubness.score.md),
[`ivi.from.indices()`](https://asalavaty.github.io/influential/reference/ivi.from.indices.md),
[`spreading.score()`](https://asalavaty.github.io/influential/reference/spreading.score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- coexpression.data
My_graph <- graph_from_data_frame(MyData)
GraphVertices <- V(My_graph)
My.vertices.IVI <- ivi(graph = My_graph, vertices = GraphVertices,
                       weights = NULL, directed = FALSE, mode = "all",
                       loops = TRUE, d = 3, scale = "range")
} # }
```
