# Hubness score

This function calculates the Hubness score of the desired nodes from a
graph. Hubness score reflects the power of each node in its surrounding
environment and is one of the major components of the IVI.

## Usage

``` r
hubness.score(
  graph,
  vertices = V(graph),
  directed = FALSE,
  mode = "all",
  loops = TRUE,
  scale = "range",
  verbose = FALSE
)
```

## Arguments

- graph:

  A graph (network) of the igraph class.

- vertices:

  A vector of desired vertices, which could be obtained by the V
  function.

- directed:

  Logical scalar, whether to directed graph is analyzed. This argument
  is ignored for undirected graphs.

- mode:

  The mode of Hubness score depending on the directedness of the graph.
  If the graph is undirected, the mode "all" should be specified.
  Otherwise, for the calculation of Hubness score based on incoming
  connections select "in" and for the outgoing connections select "out".
  Also, if all of the connections are desired, specify the "all" mode.
  Default mode is set to "all".

- loops:

  Logical; whether the loop edges are also counted.

- scale:

  Character string; the method used for scaling/normalizing the results.
  Options include 'range' (normalization within a 1-100 range),
  'z-scale' (standardization using the z-score), and 'none' (no data
  scaling). The default selection is 'range'. Opting for the 'range'
  method is suitable when exploring a single network, allowing you to
  observe the complete spectrum and distribution of node influences. In
  this case, there is no intention to establish a specific threshold for
  the outcomes. However, it is possible to identify and present the top
  hub nodes based on their rankings. Conversely, the 'z-scale' option
  proves advantageous if the aim is to compare node influences across
  multiple networks or if there is a desire to establish a threshold
  (usually z-score \> 1.645) for generating a list of the most hub nodes
  without manual intervention.

- verbose:

  Logical; whether the accomplishment of different stages of the
  algorithm should be printed (default is FALSE).

## Value

A numeric vector with the Hubness scores.

## See also

[`cent_network.vis`](https://asalavaty.github.io/influential/reference/cent_network.vis.md)

Other integrative ranking functions:
[`comp_manipulate()`](https://asalavaty.github.io/influential/reference/comp_manipulate.md),
[`exir()`](https://asalavaty.github.io/influential/reference/exir.md),
[`ivi()`](https://asalavaty.github.io/influential/reference/ivi.md),
[`ivi.from.indices()`](https://asalavaty.github.io/influential/reference/ivi.from.indices.md),
[`spreading.score()`](https://asalavaty.github.io/influential/reference/spreading.score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- coexpression.data
My_graph <- graph_from_data_frame(MyData)
GraphVertices <- V(My_graph)
Hubness.score <- hubness.score(graph = My_graph, vertices = GraphVertices,
                               directed = FALSE, mode = "all",
                               loops = TRUE, scale = "range")
} # }
```
