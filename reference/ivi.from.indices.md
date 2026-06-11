# Integrated Value of Influence (IVI)

This function calculates the IVI of the desired nodes from previously
calculated centrality measures. This function is not dependent to other
packages and the required centrality measures, namely degree centrality,
ClusterRank, betweenness centrality, Collective Influence, local
H-index, and neighborhood connectivity could have been calculated by any
means beforehand. A shiny app has also been developed for the
calculation of IVI as well as IVI-based network visualization, which is
accessible using the \`influential::runShinyApp("IVI")\` command. You
can also access the shiny app online at
https://influential.erc.monash.edu/.

## Usage

``` r
ivi.from.indices(
  DC,
  CR,
  LH_index,
  NC,
  BC,
  CI,
  scale = "range",
  verbose = FALSE
)
```

## Arguments

- DC:

  A vector containing the values of degree centrality of the desired
  vertices.

- CR:

  A vector containing the values of ClusterRank of the desired vertices.

- LH_index:

  A vector containing the values of local H-index of the desired
  vertices.

- NC:

  A vector containing the values of neighborhood connectivity of the
  desired vertices.

- BC:

  A vector containing the values of betweenness centrality of the
  desired vertices.

- CI:

  A vector containing the values of Collective Influence of the desired
  vertices.

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
[`ivi()`](https://asalavaty.github.io/influential/reference/ivi.md),
[`spreading.score()`](https://asalavaty.github.io/influential/reference/spreading.score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- centrality.measures
My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,
                                    CR = centrality.measures$CR,
                                    NC = centrality.measures$NC,
                                    LH_index = centrality.measures$LH_index,
                                    BC = centrality.measures$BC,
                                    CI = centrality.measures$CI)
                                    } # }
```
