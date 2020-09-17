  - [Overview](#overview)
  - [Network reconstruction](#network-reconstruction)
      - [From a data frame](#from-a-data-frame)
      - [From an adjacency matrix](#from-an-adjacency-matrix)
      - [From an incidence matrix](#from-an-incidence-matrix)
      - [From a SIF file](#from-a-sif-file)
  - [Calculation of centrality
    measures](#calculation-of-centrality-measures)
      - [Network vertices](#network-vertices)
      - [Degree centrality](#degree-centrality)
      - [Betweenness centrality](#betweenness-centrality)
      - [Neighborhood connectivity](#neighborhood-connectivity)
      - [H-index](#h-index)
      - [Local H-index](#local-h-index)
      - [Collective Influence](#collective-influence)
      - [ClusterRank](#clusterrank)
  - [Assessment of the association of centrality
    measures](#assessment-of-the-association-of-centrality-measures)
      - [Conditional probability of deviation from
        means](#conditional-probability-of-deviation-from-means)
      - [Nature of association (considering dependent and
        independent)](#nature-of-association-considering-dependent-and-independent)
      - [Nature of association (without considering dependence
        direction)](#nature-of-association-without-considering-dependence-direction)
  - [Identification of the most `influential` network nodes](#IVI)
      - [Integrated Value of Influence (IVI) from centrality
        measures](#integrated-value-of-influence-ivi-from-centrality-measures)
      - [Integrated Value of Influence (IVI) from a
        graph](#integrated-value-of-influence-ivi-from-a-graph)
  - [Identification of the most important network
    spreaders](#identification-of-the-most-important-network-spreaders)
  - [Identification of the most important network
    hubs](#identification-of-the-most-important-network-hubs)
  - [Ranking the influence of nodes on the topology of a network based
    on the `SIRIR`
    model](#ranking-the-influence-of-nodes-on-the-topology-of-a-network-based-on-the-sirir-model)
  - [Experimental data-based classification and ranking of top candidate
    features based on the `ExIR`
    model](#experimental-data-based-classification-and-ranking-of-top-candidate-features-based-on-the-exir-model)
      - [Assembling the Diff\_data](#assembling-the-diff_data)
      - [Preparing the Exptl\_data](#preparing-the-exptl_data)
      - [Running the `ExIR` model](#running-the-exir-model)
  - [Visualization](#visualization)
      - [Network visualization](#network-visualization)
      - [ExIR visualization](#exir-visualization)
  - [References](#references)

# Overview

<a href='https://github.com/asalavaty/influential'><img src='../man/figures/logo.png' align="right" height="198" /></a>

`influential` is an R package mainly for the identification of the most
influential nodes in a network as well as the classification and ranking
of top candidate features.The `influential` package contains several
functions that could be categorized into six groups according to their
purpose:

  - Network reconstruction
  - Calculation of centrality measures
  - Assessment of the association of centrality measures
  - Identification of the most `influential` network nodes
  - Experimental data-based classification and ranking of features
  - Visualization

The sections below introduce these six categories.

-----

``` r
library(influential)
```

# Network reconstruction

Three functions have been obtained from the `igraph`\[1\] R package for
the reconstruction of networks.

## From a data frame

In the data frame the first and second columns should be composed of
source and target nodes.  
A sample appropriate data frame is brought below:

<table>
<thead>
<tr class="header">
<th style="text-align: left;">lncRNA</th>
<th style="text-align: left;">Coexpressed.Gene</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">ADAMTS9-AS2</td>
<td style="text-align: left;">A2M</td>
</tr>
<tr class="even">
<td style="text-align: left;">ADAMTS9-AS2</td>
<td style="text-align: left;">ABCA6</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ADAMTS9-AS2</td>
<td style="text-align: left;">ABCA8</td>
</tr>
<tr class="even">
<td style="text-align: left;">ADAMTS9-AS2</td>
<td style="text-align: left;">ABCA9</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ADAMTS9-AS2</td>
<td style="text-align: left;">ABI3BP</td>
</tr>
<tr class="even">
<td style="text-align: left;">ADAMTS9-AS2</td>
<td style="text-align: left;">AC093110.3</td>
</tr>
</tbody>
</table>

This is a co-expression dataset obtained from a paper by Salavaty et
al.\[2\]

``` r
# Preparing the data
MyData <- coexpression.data

# Reconstructing the graph
My_graph <- graph_from_data_frame(d=MyData)
```

If you look at the class of `My_graph` you should see that it has an
`igraph` class:

``` r
class(My_graph)
#> [1] "igraph"
```

## From an adjacency matrix

A sample appropriate adjacency matrix is brought below:

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">LINC00891</th>
<th style="text-align: right;">LINC00968</th>
<th style="text-align: right;">LINC00987</th>
<th style="text-align: right;">LINC01506</th>
<th style="text-align: right;">MAFG-AS1</th>
<th style="text-align: right;">MIR497HG</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">LINC00891</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">LINC00968</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">LINC00987</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">LINC01506</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">MAFG-AS1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">MIR497HG</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
</tbody>
</table>

  - Note that the matrix has the same number of rows and columns.

<!-- end list -->

``` r
# Preparing the data
MyData <- coexpression.adjacency        

# Reconstructing the graph
My_graph <- graph_from_adjacency_matrix(MyData)        
```

## From an incidence matrix

A sample appropriate incidence matrix is brought below:

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">Gene_1</th>
<th style="text-align: right;">Gene_2</th>
<th style="text-align: right;">Gene_3</th>
<th style="text-align: right;">Gene_4</th>
<th style="text-align: right;">Gene_5</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">cell_1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td style="text-align: left;">cell_2</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">cell_3</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">cell_4</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0</td>
</tr>
</tbody>
</table>

``` r
# Reconstructing the graph
My_graph <- graph_from_adjacency_matrix(MyData)        
```

## From a SIF file

SIF is the common output format of the Cytoscape software.

``` r
# Reconstructing the graph
My_graph <- sif2igraph(Path = "Sample SIF.sif")        

class(My_graph)
#> [1] "igraph"
```

<a href="#top">Back to top</a>

# Calculation of centrality measures

To calculate the centrality of nodes within a network several different
options are available. The following sections describe how to obtain the
names of network nodes and use different functions to calculate the
centrality of nodes within a network. Although several centrality
functions are provided, we recommend the [IVI](#IVI) for the
identification of the most `influential` nodes within a network.

## Network vertices

Network vertices (nodes) are required in order to calculate their
centrality measures. Thus, before calculation of network centrality
measures we need to obtain the name of required network vertices. To
this end, we use the `V` function, which is obtained from the `igraph`
package. However, you may provide a character vector of the name of your
desired nodes manually.

  - Note in many of the centrality index functions the entire network
    nodes are assessed if no vector of desired vertices is provided.

<!-- end list -->

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
My_graph_vertices <- V(My_graph)        

head(My_graph_vertices)
#> + 6/794 vertices, named, from 775cff6:
#> [1] ADAMTS9-AS2 C8orf34-AS1 CADM3-AS1   FAM83A-AS1  FENDRR      LANCL1-AS1
```

## Degree centrality

Degree centrality is the most commonly used local centrality measure
which could be calculated via the `degree` function obtained from the
`igraph` package.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculating degree centrality
My_graph_degree <- degree(My_graph, v = GraphVertices, normalized = FALSE) 

head(My_graph_degree)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>         172         121         168          26         189         176
```

Degree centrality could be also calculated for *directed* graphs via
specifying the `mode` parameter.

## Betweenness centrality

Betweenness centrality, like degree centrality, is one of the most
commonly used centrality measures but is representative of the global
centrality of a node. This centrality metric could also be calculated
using a function obtained from the `igraph` package.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculating betweenness centrality
My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,    
                                    directed = FALSE, normalized = FALSE)

head(My_graph_betweenness)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>   21719.857   28185.199   26946.625    2940.467   33333.369   21830.511
```

Betweenness centrality could be also calculated for *directed* and/or
*weighted* graphs via specifying the `directed` and `weights`
parameters, respectively.

## Neighborhood connectivity

Neighborhood connectivity is one of the other important centrality
measures that reflect the semi-local centrality of a node. This
centrality measure was first represented in a Science paper\[3\] in 2002
and is for the first time calculable in R environment via the
`influential` package.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculating neighborhood connectivity
neighrhood.co <- neighborhood.connectivity(graph = My_graph,    
                                           vertices = GraphVertices,
                                           mode = "all")

head(neighrhood.co)
#>  ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>   11.290698    4.983471    7.970238    3.000000   15.153439   13.465909
```

Neighborhood connectivity could be also calculated for *directed* graphs
via specifying the `mode` parameter.

## H-index

H-index is H-index is another semi-local centrality measure that was
inspired from its application in assessing the impact of researchers and
is for the first time calculable in R environment via the `influential`
package.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculating H-index
h.index <- h_index(graph = My_graph,    
                   vertices = GraphVertices,
                   mode = "all")

head(h.index)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>          11           9          11           2          12          12
```

H-index could be also calculated for *directed* graphs via specifying
the `mode` parameter.

## Local H-index

Local H-index (LH-index) is a semi-local centrality measure and an
improved version of H-index centrality that leverages the H-index to the
second order neighbors of a node and is for the first time calculable in
R environment via the `influential` package.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculating Local H-index
lh.index <- lh_index(graph = My_graph,    
                   vertices = GraphVertices,
                   mode = "all")

head(lh.index)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>        1165         446         994          34        1289        1265
```

Local H-index could be also calculated for *directed* graphs via
specifying the `mode` parameter.

## Collective Influence

Collective Influence (CI) is a global centrality measure that calculates
the product of the reduced degree (degree - 1) of a node and the total
reduced degree of all nodes at a distance d from the node. This
centrality measure is for the first time provided in an R package.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculating Collective Influence
ci <- collective.influence(graph = My_graph,    
                          vertices = GraphVertices,
                          mode = "all", d=3)

head(ci)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>        9918       70560       39078         675       10716        7350
```

Collective Influence could be also calculated for *directed* graphs via
specifying the `mode` parameter.

## ClusterRank

ClusterRank is a local centrality measure that makes a connection
between local and semi-local characteristics of a node and at the same
time removes the negative effects of local clustering.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculating ClusterRank
cr <- clusterRank(graph = My_graph,    
                  vids = GraphVertices,
                  directed = FALSE, loops = TRUE)

head(cr)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>   63.459812    5.185675   21.111776    1.280000  135.098278   81.255195
```

ClusterRank could be also calculated for *directed* graphs via
specifying the `directed` parameter.

<a href="#top">Back to top</a>

# Assessment of the association of centrality measures

## Conditional probability of deviation from means

The function `cond.prob.analysis` assesses the conditional probability
of deviation of two centrality measures (or any other two continuous
variables) from their corresponding means in opposite directions.

``` r
# Preparing the data
MyData <- centrality.measures        

# Assessing the conditional probability
My.conditional.prob <- cond.prob.analysis(data = MyData,       
                                          nodes.colname = rownames(MyData),
                                          Desired.colname = "BC",
                                          Condition.colname = "NC")

print(My.conditional.prob)
#> $ConditionalProbability
#> [1] 51.61871
#> 
#> $ConditionalProbability_split.half.sample
#> [1] 51.73611
```

  - As you can see in the results, the whole data is also randomly
    splitted into half in order to further test the validity of
    conditional probability assessments.
  - *The higher the conditional probability the more two centrality
    measures behave in contrary manners*.

## Nature of association (considering dependent and independent)

The function `double.cent.assess` could be used to automatically assess
both the distribution mode of centrality measures (two continuous
variables) and the nature of their association. The analyses done
through this formula are as follows:

1.  **Normality assessment**:
      - Variables with **lower than** 5000 observations: *Shapiro-Wilk
        test*
      - Variables with **over** 5000 observations: *Anderson-Darling
        test* <br><br>
2.  **Assessment of non-linear/non-monotonic correlation**:
      - *Non-linearity assessment*: Fitting a generalized additive model
        (GAM) with integrated smoothness approximations using the `mgcv`
        package <br><br>
      - *Non-monotonicity assessment*: Comparing the squared
        coefficients of the correlation based on Spearman’s rank
        correlation analysis and ranked regression test with non-linear
        splines.
          - Squared coefficient of Spearman’s rank correlation **\>**
            R-squared ranked regression with non-linear splines:
            *Monotonic*
          - Squared coefficient of Spearman’s rank correlation **\<**
            R-squared ranked regression with non-linear splines:
            *Non-monotonic* <br><br>
3.  **Dependence assessment**:
      - *Hoeffding’s independence test*: Hoeffding’s test of
        independence is a test based on the population measure of
        deviation from independence which computes a D Statistics
        ranging from -0.5 to 1: Greater D values indicate a higher
        dependence between variables.
      - *Descriptive non-linear non-parametric dependence test*: This
        assessment is based on non-linear non-parametric statistics
        (NNS) which outputs a dependence value ranging from 0 to 1. For
        further details please refer to the NNS R package\[4\]: Greater
        values indicate a higher dependence between variables. <br><br>
4.  **Correlation assessment**: As the correlation between most of the
    centrality measures follows a non-monotonic form, this part of the
    assessment is done based on the NNS statistics which itself
    calculates the correlation based on partial moments and outputs a
    correlation value ranging from -1 to 1. For further details please
    refer to the NNS R package. <br><br>
5.  **Assessment of conditional probability of deviation from means**
    This step assesses the conditional probability of deviation of two
    centrality measures (or any other two continuous variables) from
    their corresponding means in opposite directions.
      - The independent centrality measure (variable) is considered as
        the condition variable and the other as the desired one.
      - As you will see in the results, the whole data is also randomly
        splitted into half in order to further test the validity of
        conditional probability assessments.
      - *The higher the conditional probability the more two centrality
        measures behave in contrary manners*.

<!-- end list -->

``` r
# Preparing the data
MyData <- centrality.measures        

# Association assessment
My.metrics.assessment <- double.cent.assess(data = MyData,       
                                            nodes.colname = rownames(MyData),
                                            dependent.colname = "BC",
                                            independent.colname = "NC")

print(My.metrics.assessment)
#> $Summary_statistics
#>         BC NC
#> Min.              0.000000000                   1.2000
#> 1st Qu.           0.000000000                  66.0000
#> Median            0.000000000                 156.0000
#> Mean              0.005813357                 132.3443
#> 3rd Qu.           0.000340000                 179.3214
#> Max.              0.529464720                 192.0000
#> 
#> $Normality_results
#>                               p.value
#> BC    1.415450e-50
#> NC 9.411737e-30
#> 
#> $Dependent_Normality
#> [1] "Non-normally distributed"
#> 
#> $Independent_Normality
#> [1] "Non-normally distributed"
#> 
#> $GAM_nonlinear.nonmonotonic.results
#>      edf  p-value 
#> 8.992406 0.000000 
#> 
#> $Association_type
#> [1] "nonlinear-nonmonotonic"
#> 
#> $HoeffdingD_Statistic
#>         D_statistic P_value
#> Results  0.01770279   1e-08
#> 
#> $Dependence_Significance
#>                       Hoeffding
#> Results Significantly dependent
#> 
#> $NNS_dep_results
#>         Correlation Dependence
#> Results  -0.7948106  0.8647164
#> 
#> $ConditionalProbability
#> [1] 55.35386
#> 
#> $ConditionalProbability_split.half.sample
#> [1] 55.90331
```

**Note**: It should also be noted that as a single regression line does
not fit all models with a certain degree of freedom, based on the size
and correlation mode of the variables provided, this function might
return an error due to incapability of running step 2. In this case, you
may follow each step manually or as an alternative run the other
function named `double.cent.assess.noRegression` which does not perform
any regression test and consequently it is not required to determine the
dependent and independent variables.

## Nature of association (without considering dependence direction)

The function `double.cent.assess.noRegression` could be used to
automatically assess both the distribution mode of centrality measures
(two continuous variables) and the nature of their association. The
analyses done through this formula are as follows:

1.  **Normality assessment**:
      - Variables with **lower than** 5000 observations: *Shapiro-Wilk
        test*
      - Variables with **over** 5000 observations: *Anderson–Darling
        test* <br><br>
2.  **Dependence assessment**:
      - *Hoeffding’s independence test*: Hoeffding’s test of
        independence is a test based on the population measure of
        deviation from independence which computes a D Statistics
        ranging from -0.5 to 1: Greater D values indicate a higher
        dependence between variables.
      - *Descriptive non-linear non-parametric dependence test*: This
        assessment is based on non-linear non-parametric statistics
        (NNS) which outputs a dependence value ranging from 0 to 1. For
        further details please refer to the NNS R package: Greater
        values indicate a higher dependence between variables. <br><br>
3.  **Correlation assessment**: As the correlation between most of the
    centrality measures follows a non-monotonic form, this part of the
    assessment is done based on the NNS statistics which itself
    calculates the correlation based on partial moments and outputs a
    correlation value ranging from -1 to 1. For further details please
    refer to the NNS R package. <br><br>
4.  **Assessment of conditional probability of deviation from means**
    This step assesses the conditional probability of deviation of two
    centrality measures (or any other two continuous variables) from
    their corresponding means in opposite directions.
      - The `centrality2` variable is considered as the condition
        variable and the other (`centrality1`) as the desired one.
      - As you will see in the results, the whole data is also randomly
        splitted into half in order to further test the validity of
        conditional probability assessments.
      - *The higher the conditional probability the more two centrality
        measures behave in contrary manners*.

<!-- end list -->

``` r
# Preparing the data
MyData <- centrality.measures        

# Association assessment
My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,       
                                                         nodes.colname = rownames(MyData),
                                                         centrality1.colname = "BC",
                                                         centrality2.colname = "NC")

print(My.metrics.assessment)
#> $Summary_statistics
#>         BC NC
#> Min.              0.000000000                   1.2000
#> 1st Qu.           0.000000000                  66.0000
#> Median            0.000000000                 156.0000
#> Mean              0.005813357                 132.3443
#> 3rd Qu.           0.000340000                 179.3214
#> Max.              0.529464720                 192.0000
#> 
#> $Normality_results
#>                               p.value
#> BC    1.415450e-50
#> NC 9.411737e-30
#> 
#> $Centrality1_Normality
#> [1] "Non-normally distributed"
#> 
#> $Centrality2_Normality
#> [1] "Non-normally distributed"
#> 
#> $HoeffdingD_Statistic
#>         D_statistic P_value
#> Results  0.01770279   1e-08
#> 
#> $Dependence_Significance
#>                       Hoeffding
#> Results Significantly dependent
#> 
#> $NNS_dep_results
#>         Correlation Dependence
#> Results  -0.7948106  0.8647164
#> 
#> $ConditionalProbability
#> [1] 55.35386
#> 
#> $ConditionalProbability_split.half.sample
#> [1] 55.68163
```

<a href="#top">Back to top</a>

# Identification of the most `influential` network nodes

**IVI** : `IVI` is the first integrative method for the identification
of network most influential nodes in a way that captures all network
topological dimensions. The `IVI` formula integrates the most important
local (i.e. degree centrality and ClusterRank), semi-local
(i.e. neighborhood connectivity and local H-index) and global
(i.e. betweenness centrality and collective influence) centrality
measures in such a way that both synergize their effects and remove
their biases.

## Integrated Value of Influence (IVI) from centrality measures

``` r
# Preparing the data
MyData <- centrality.measures        

# Calculation of IVI
My.vertices.IVI <- ivi.from.indices(DC = MyData$DC,       
                                   CR = MyData$CR,
                                   NC = MyData$NC,
                                   LH_index = MyData$LH_index,
                                   BC = MyData$BC,
                                   CI = MyData$CI)

head(My.vertices.IVI)
#> [1] 24.670056  8.344337 18.621049  1.017768 29.437028 33.512598
```

## Integrated Value of Influence (IVI) from a graph

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculation of IVI
My.vertices.IVI <- ivi(graph = My_graph, vertices = GraphVertices, 
                       weights = NULL, directed = FALSE, mode = "all",
                       loops = TRUE, d = 3, scaled = TRUE)

head(My.vertices.IVI)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>    39.53878    19.94999    38.20524     1.12371   100.00000    47.49356
```

IVI could be also calculated for *directed* and/or *weighted* graphs via
specifying the `directed`, `mode`, and `weights` parameters.

Check out [our paper](https://doi.org/10.1016/j.patter.2020.100052)\[5\]
for a more complete description of the IVI formula and all of its
underpinning methods and analyses.

The following tutorial video demonstrates how to simply calculate the
IVI value of all of the nodes within a network.

<iframe width="100%" height="400" src="https://www.youtube.com/embed/nF9qe4dvuJc
" frameborder="0" allowfullscreen>

</iframe>

<a href="#top">Back to top</a>

# Identification of the most important network spreaders

Sometimes we seek to identify not necessarily the most influential nodes
but the nodes with most potential in spreading of information throughout
the network.

**Spreading score** : `spreading.score` is an integrative score made up
of four different centrality measures including ClusterRank,
neighborhood connectivity, betweenness centrality, and collective
influence. Also, Spreading score reflects the spreading potential of
each node within a network and is one of the major components of the
`IVI`.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculation of Spreading score
Spreading.score <- spreading.score(graph = My_graph,     
                                   vertices = GraphVertices, 
                                   weights = NULL, directed = FALSE, mode = "all",
                                   loops = TRUE, d = 3, scaled = TRUE)

head(Spreading.score)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>   42.932497   38.094111   45.114648    1.587262  100.000000   49.193292 
```

Spreading score could be also calculated for *directed* and/or
*weighted* graphs via specifying the `directed`, `mode`, and `weights`
parameters.

<a href="#top">Back to top</a>

# Identification of the most important network hubs

In some cases we want to identify not the nodes with the most
sovereignty in their surrounding local environments.

**Hubness score** : `hubness.score` is an integrative score made up of
two different centrality measures including degree centrality and local
H-index. Also, Hubness score reflects the power of each node in its
surrounding environment and is one of the major components of the `IVI`.

``` r
# Preparing the data
MyData <- coexpression.data        

# Reconstructing the graph
My_graph <- graph_from_data_frame(MyData)        

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculation of Hubness score
Hubness.score <- hubness.score(graph = My_graph,     
                                   vertices = GraphVertices, 
                                   directed = FALSE, mode = "all",
                                   loops = TRUE, scaled = TRUE)

head(Hubness.score)
#> ADAMTS9-AS2 C8orf34-AS1   CADM3-AS1  FAM83A-AS1      FENDRR  LANCL1-AS1 
#>   84.299719   46.741660   77.441514    8.437142   92.870451   88.734131
```

Spreading score could be also calculated for *directed* graphs via
specifying the `directed` and `mode` parameters.

<a href="#top">Back to top</a>

# Ranking the influence of nodes on the topology of a network based on the `SIRIR` model

**SIRIR** : `SIRIR` is achieved by the integration
susceptible-infected-recovered (SIR) model with the leave-one-out cross
validation technique and ranks network nodes based on their true
universal influence on the network topology and spread of information.
One of the applications of this function is the assessment of
performance of a novel algorithm in identification of network
influential nodes.

``` r
# Reconstructing the graph
My_graph <-  sif2igraph(Path = "Sample SIF.sif")       

# Extracting the vertices
GraphVertices <- V(My_graph)        

# Calculation of influence rank
Influence.Ranks <- sirir(graph = My_graph,     
                                   vertices = GraphVertices, 
                                   beta = 0.5, gamma = 1, no.sim = 10, seed = 1234)
```

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">difference.value</th>
<th style="text-align: right;">rank</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">MRAP</td>
<td style="text-align: right;">49.7</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td style="text-align: left;">FOXM1</td>
<td style="text-align: right;">49.5</td>
<td style="text-align: right;">2</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ATAD2</td>
<td style="text-align: right;">49.5</td>
<td style="text-align: right;">2</td>
</tr>
<tr class="even">
<td style="text-align: left;">POSTN</td>
<td style="text-align: right;">49.4</td>
<td style="text-align: right;">4</td>
</tr>
<tr class="odd">
<td style="text-align: left;">CDC7</td>
<td style="text-align: right;">49.3</td>
<td style="text-align: right;">5</td>
</tr>
<tr class="even">
<td style="text-align: left;">ZWINT</td>
<td style="text-align: right;">42.1</td>
<td style="text-align: right;">6</td>
</tr>
<tr class="odd">
<td style="text-align: left;">MKI67</td>
<td style="text-align: right;">41.9</td>
<td style="text-align: right;">7</td>
</tr>
<tr class="even">
<td style="text-align: left;">FN1</td>
<td style="text-align: right;">41.9</td>
<td style="text-align: right;">7</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ASPM</td>
<td style="text-align: right;">41.8</td>
<td style="text-align: right;">9</td>
</tr>
<tr class="even">
<td style="text-align: left;">ANLN</td>
<td style="text-align: right;">41.8</td>
<td style="text-align: right;">9</td>
</tr>
</tbody>
</table>

<a href="#top">Back to top</a>

# Experimental data-based classification and ranking of top candidate features based on the `ExIR` model

**ExIR** : `ExIR` is a model for the classification and ranking of top
candidate features. The input data could come from any type of
experiment such as transcriptomics and proteomics. This model is based
on multi-level filtration and scoring based on several supervised and
unsupervised analyses followed by the classification and integrative
ranking of top candidate features. Using this function and depending on
the input data and specified arguments, the user can get one to four
tables including:

  - **Drivers**: Prioritized drivers are supposed to have the highest
    impact on the progression of a biological process or disease under
    investigation.
  - **Biomarkers**: Prioritized biomarkers are supposed to have the
    highest sensitivity to different conditions under investigation and
    the severity of each condition.
  - **DE-mediators**: Prioritized DE-mediators are those features that
    are differentially expressed/abundant but in a fluctuating manner
    and play mediatory roles between drivers.
  - **nonDE-mediators**: Prioritized nonDE-mediators are those features
    that are not differentially expressed/abundant but have associations
    with and play mediatory roles between drivers.

First, prepare your data. Suppose we have the data for time-course
transcriptomics and we have previously performed differential expression
analysis for each step-wise pair of time-points. Also, we have performed
trajectory analysis to identify the genes that have significant
alterations across all time-points.

``` r
# Prepare sample data
gene.names <- paste("gene", c(1:1000), sep = "_")

set.seed(60)
tp2.vs.tp1.DEGs <- data.frame(logFC = runif(n = 90, min = -5, max = 5),
                              FDR = runif(n = 90, min = 0.0001, max = 0.049))

set.seed(60)
rownames(tp2.vs.tp1.DEGs) <- sample(gene.names, size = 90)

set.seed(70)
tp3.vs.tp2.DEGs <- data.frame(logFC = runif(n = 121, min = -3, max = 6),
                              FDR = runif(n = 121, min = 0.0011, max = 0.039))

set.seed(70)
rownames(tp3.vs.tp2.DEGs) <- sample(gene.names, size = 121)

set.seed(80)
regression.data <- data.frame(R_squared = runif(n = 65, min = 0.1, max = 0.85))

set.seed(80)
rownames(regression.data) <- sample(gene.names, size = 65)
```

## Assembling the Diff\_data

Use the function `diff_data.assembly` to automatically generate the
Diff\_data table for the `ExIR` model.

``` r
my_Diff_data <- diff_data.assembly(tp2.vs.tp1.DEGs,
                                   tp3.vs.tp2.DEGs,
                                   regression.data)

knitr::kable(my_Diff_data[c(1:10),])
```

## Preparing the Exptl\_data

Now, prepare a sample normalized experimental data matrix

``` r
set.seed(60)
MyExptl_data <- matrix(data = runif(n = 50000, min = 2, max = 100), 
                       nrow = 50, ncol = 1000, 
                       dimnames = list(c(paste("cancer_sample", c(1:25), sep = "_"),
                                         paste("normal_sample", c(1:25), sep = "_")), 
                                       gene.names))

# Log transform the data to bring them closer to normal distribution
MyExptl_data <- log2(MyExptl_data)  

knitr::kable(MyExptl_data[c(1:5),c(1:5)])
```

Now add the “condition” column to the Exptl\_data table.

``` r
MyExptl_data <- as.data.frame(MyExptl_data)
MyExptl_data$condition <- c(rep("C", 25), rep("N", 25))
```

## Running the `ExIR` model

Finally, prepare the other required input data for the `ExIR` model.

``` r

#The table of differential/regression previously prepared
my_Diff_data

#The column indices of differential values in the Diff_data table
Diff_value <- c(1,3) 

#The column indices of regression values in the Diff_data table
Regr_value <- 5

#The column indices of significance (P-value/FDR) values in 
# the Diff_data table
Sig_value <- c(2,4) 

#The matrix/data frame of normalized experimental
# data previously prepared
MyExptl_data

#The name of the column delineating the conditions of
# samples in the Exptl_data matrix
Condition_colname <- "condition"

#The desired list of features
set.seed(60)
MyDesired_list <- sample(gene.names, size = 200)  #Optional

#Running the ExIR model
My.exir <- exir(Desired_list = MyDesired_list,
Diff_data = my_Diff_data, Diff_value = Diff_value,
Regr_value = Regr_value, Sig_value = Sig_value,
Exptl_data = MyExptl_data, Condition_colname = Condition_colname,
seed = 60, verbose = FALSE)

names(My.exir)
#> [1] "Driver table"          "DE-mediator table"     "nonDE-mediator table"
#> [4] "Biomarker table"

class(My.exir)
#> [1] "ExIR_Result"
```

Have a look at the heads of the outputs of ExIR:

  - **Drivers**

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">Score</th>
<th style="text-align: right;">Z.score</th>
<th style="text-align: right;">Rank</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">P.adj</th>
<th style="text-align: left;">Type</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">gene_738</td>
<td style="text-align: right;">100.00000</td>
<td style="text-align: right;">5.380024</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3.723800e-08</td>
<td style="text-align: right;">0.0000040</td>
<td style="text-align: left;">Decelerator</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_520</td>
<td style="text-align: right;">70.84240</td>
<td style="text-align: right;">3.569259</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">1.789961e-04</td>
<td style="text-align: right;">0.0095763</td>
<td style="text-align: left;">Accelerator</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_78</td>
<td style="text-align: right;">55.87662</td>
<td style="text-align: right;">2.639844</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">4.147204e-03</td>
<td style="text-align: right;">0.1479169</td>
<td style="text-align: left;">Decelerator</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_53</td>
<td style="text-align: right;">52.47212</td>
<td style="text-align: right;">2.428416</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">7.582469e-03</td>
<td style="text-align: right;">0.2028311</td>
<td style="text-align: left;">Accelerator</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_886</td>
<td style="text-align: right;">50.63736</td>
<td style="text-align: right;">2.314472</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">1.032092e-02</td>
<td style="text-align: right;">0.2208677</td>
<td style="text-align: left;">Accelerator</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_883</td>
<td style="text-align: right;">45.78593</td>
<td style="text-align: right;">2.013185</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2.204757e-02</td>
<td style="text-align: right;">0.3527402</td>
<td style="text-align: left;">Accelerator</td>
</tr>
</tbody>
</table>

-----

  - **Biomarkers**

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">Score</th>
<th style="text-align: right;">Z.score</th>
<th style="text-align: right;">Rank</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">P.adj</th>
<th style="text-align: left;">Type</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">gene_711</td>
<td style="text-align: right;">100.000000</td>
<td style="text-align: right;">9.1261490</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3.548963e-20</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: left;">Up-regulated</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_245</td>
<td style="text-align: right;">37.107082</td>
<td style="text-align: right;">3.1933225</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">7.032290e-04</td>
<td style="text-align: right;">0.0376227</td>
<td style="text-align: left;">Down-regulated</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_705</td>
<td style="text-align: right;">24.598890</td>
<td style="text-align: right;">2.0133974</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">2.203642e-02</td>
<td style="text-align: right;">0.5842342</td>
<td style="text-align: left;">Up-regulated</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_910</td>
<td style="text-align: right;">23.430181</td>
<td style="text-align: right;">1.9031504</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2.851046e-02</td>
<td style="text-align: right;">0.5842342</td>
<td style="text-align: left;">Down-regulated</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_895</td>
<td style="text-align: right;">13.919517</td>
<td style="text-align: right;">1.0059887</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">1.572105e-01</td>
<td style="text-align: right;">0.5842342</td>
<td style="text-align: left;">Down-regulated</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_800</td>
<td style="text-align: right;">5.843283</td>
<td style="text-align: right;">0.2441399</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">4.035613e-01</td>
<td style="text-align: right;">0.5842342</td>
<td style="text-align: left;">Down-regulated</td>
</tr>
</tbody>
</table>

-----

  - **DE-mediators**

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">Score</th>
<th style="text-align: right;">Z.score</th>
<th style="text-align: right;">Rank</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">P.adj</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">gene_315</td>
<td style="text-align: right;">100.000000</td>
<td style="text-align: right;">1.6261374</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.05196021</td>
<td style="text-align: right;">0.3117613</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_963</td>
<td style="text-align: right;">68.097002</td>
<td style="text-align: right;">0.8442286</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">0.19927085</td>
<td style="text-align: right;">0.5978125</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_49</td>
<td style="text-align: right;">19.104450</td>
<td style="text-align: right;">-0.3565272</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">0.63927712</td>
<td style="text-align: right;">0.7882165</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_682</td>
<td style="text-align: right;">10.102668</td>
<td style="text-align: right;">-0.5771514</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">0.71808142</td>
<td style="text-align: right;">0.7882165</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_493</td>
<td style="text-align: right;">3.603502</td>
<td style="text-align: right;">-0.7364391</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">0.76926825</td>
<td style="text-align: right;">0.7882165</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_898</td>
<td style="text-align: right;">1.000000</td>
<td style="text-align: right;">-0.8002482</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">0.78821650</td>
<td style="text-align: right;">0.7882165</td>
</tr>
</tbody>
</table>

-----

  - **nonDE-mediators**

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">Score</th>
<th style="text-align: right;">Z.score</th>
<th style="text-align: right;">Rank</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">P.adj</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">gene_117</td>
<td style="text-align: right;">100.00000</td>
<td style="text-align: right;">5.236991</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">8.160804e-08</td>
<td style="text-align: right;">0.0000071</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_516</td>
<td style="text-align: right;">57.59469</td>
<td style="text-align: right;">2.693473</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">3.535593e-03</td>
<td style="text-align: right;">0.1468437</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_262</td>
<td style="text-align: right;">54.38275</td>
<td style="text-align: right;">2.500817</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">6.195351e-03</td>
<td style="text-align: right;">0.1468437</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_974</td>
<td style="text-align: right;">53.87269</td>
<td style="text-align: right;">2.470223</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">6.751434e-03</td>
<td style="text-align: right;">0.1468437</td>
</tr>
<tr class="odd">
<td style="text-align: left;">gene_441</td>
<td style="text-align: right;">46.60681</td>
<td style="text-align: right;">2.034408</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">2.095524e-02</td>
<td style="text-align: right;">0.3646212</td>
</tr>
<tr class="even">
<td style="text-align: left;">gene_533</td>
<td style="text-align: right;">40.91995</td>
<td style="text-align: right;">1.693304</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">4.519882e-02</td>
<td style="text-align: right;">0.6553829</td>
</tr>
</tbody>
</table>

<a href="#top">Back to top</a>

# Visualization

Two visualization functions have been developed for the illustration of
the results generated by other functions of the `influential` R package.
These functions including `cent_network.vis` and `exir.vis` and their
use are ellaborated in detail in the following sections.

## Network visualization

The `cent_network.vis` is a function for the visualization of a network
based on applying a centrality measure to the size and color of network
nodes. The centrality of network nodes could be calculated by any means
and based on any centrality index. Here, we demonstrate the
visualization of a network according to [IVI](#IVI) values.

``` r
# Reconstructing the graph
set.seed(70)
My_graph <-  igraph::sample_gnm(n = 50, m = 120, directed = TRUE)

# Calculating the IVI values
My_graph_IVI <- ivi(My_graph, directed = TRUE)

# Visualizing the graph based on IVI values
My_graph_IVI_Vis <- cent_network.vis(graph = My_graph,
                                     cent.metric = My_graph_IVI,
                                     directed = TRUE,
                                     plot.title = "IVI-based Network",
                                     legend.title = "IVI value")

My_graph_IVI_Vis
```

![](../man/figures/IVI_Net_Vis.png)

## ExIR visualization

<a href="#top">Back to top</a>

# References

1.  Csardi G., Nepusz T. *The igraph software package for complex
    network research*. InterJournal. 2006; (1695).

2.  Salavaty A, Rezvani Z, Najafi A. *Survival analysis and functional
    annotation of long non-coding RNAs in lung adenocarcinoma*. J Cell
    Mol Med. 2019;23:5600–5617.
    ([PMID: 31211495](https://www.ncbi.nlm.nih.gov/pubmed/31211495)).

3.  Maslov S., Sneppen K. *Specificity and stability in topology of
    protein networks*. Science. 2002; 296: 910-913
    ([PMID:11988575](https://www.ncbi.nlm.nih.gov/pubmed/11988575))

4.  NNS: Nonlinear Nonparametric Statistics.
    [CRAN](https://cran.r-project.org/package=NNS)

5.  Salavaty A, Ramialison M, Currie PD. *Integrated Value of Influence:
    An Integrative Method for the Identification of the Most Influential
    Nodes within Networks*. Patterns. 2020.08.14 ([Read
    online](https://doi.org/10.1016/j.patter.2020.100052))
