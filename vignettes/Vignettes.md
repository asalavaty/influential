Vignettes of the `influential`
================

## Overview

The influential package contains several functions that could be
categorized into four groups according to their purpose:

  - Network reconstruction
  - Calculation of centrality measures
  - Assessment of the association of centrality measures
  - Identification of the most `influential` network nodes

-----

``` r
library(influential)
```

## Network reconstruction

Two functions have been obtained from the `igraph` package for the
reconstruction of networks.

  - ### From a data frame
    
    In the data frame the first and second columns should be composed of
    source and target nodes.
    
    A sample appropriate data frame is brought below:

| lncRNA      | Coexpressed.Gene |
| :---------- | :--------------- |
| ADAMTS9-AS2 | A2M              |
| ADAMTS9-AS2 | ABCA6            |
| ADAMTS9-AS2 | ABCA8            |
| ADAMTS9-AS2 | ABCA9            |
| ADAMTS9-AS2 | ABI3BP           |
| ADAMTS9-AS2 | AC093110.3       |

This is a co-expression dataset obtained from
[PMID: 31211495](https://www.ncbi.nlm.nih.gov/pubmed/31211495).

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(d=MyData)        # Reconstructing the graph
```

If you look at the class of `My_graph` you should see that it has an
`igraph` class:

``` r
class(My_graph)
#> [1] "igraph"
```

-----

  - ### From an adjacency matrix

A sample appropriate data frame is brought
below:

|             | ADAMTS9-AS2 | C8orf34-AS1 | CADM3-AS1 | FAM83A-AS1 | FENDRR | LANCL1-AS1 | LINC00092 | LINC00467 | LINC00857 | LINC00891 |
| ----------- | ----------: | ----------: | --------: | ---------: | -----: | ---------: | --------: | --------: | --------: | --------: |
| ADAMTS9-AS2 |           0 |           0 |         0 |          0 |      0 |          1 |         0 |         0 |         0 |         1 |
| C8orf34-AS1 |           0 |           0 |         0 |          0 |      0 |          0 |         0 |         0 |         0 |         1 |
| CADM3-AS1   |           0 |           0 |         0 |          0 |      1 |          0 |         0 |         0 |         0 |         0 |
| FAM83A-AS1  |           0 |           0 |         0 |          0 |      0 |          0 |         0 |         0 |         1 |         0 |
| FENDRR      |           0 |           0 |         0 |          0 |      0 |          0 |         0 |         0 |         0 |         0 |
| LANCL1-AS1  |           0 |           0 |         0 |          0 |      1 |          0 |         0 |         0 |         0 |         0 |

``` r
MyData <- coexpression.adjacency        # Preparing the data

My_graph <- graph_from_adjacency_matrix(MyData)        # Reconstructing the graph
```

-----

  - ### Network vertices

Network vertices (nodes) are required in order to calculate their
centrality measures. Thus, before calculation of network centrality
measures we need to obtain the name of required network vertices. To
this end, we use the `V` function, which is obtained from the `igraph`
package. However, you may provide a character vector of the name of your
desired nodes manually.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

My_graph_vertices <- V(My_graph)        # Extracting the vertices
```

-----

## Calculation of centrality measures

  - ### Degree centrality
    
    Degree centrality is the most commonly used local centrality measure
    which could be calculated via the `degree` function obtained from
    the `igraph` package.

<!-- end list -->

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

My_graph_degree <- degree(My_graph, v = GraphVertices, normalized = FALSE) # Calculating degree centrality
```

Degree centrality could be also calculated for *directed* graphs via
specifying the `mode` parameter.

  - ### Betweenness centrality
    
    Betweenness centrality, like degree centrality, is one of the most
    commonly used centrality measures but is representative of the
    global centrality of a node. This centrality metric could also be
    calculated using a function obtained from the `igraph` package.

<!-- end list -->

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,    # Calculating betweenness centrality
                                    directed = FALSE, normalized = FALSE)
```

Betweenness centrality could be also calculated for *directed* and/or
*weighted* graphs via specifying the `directed` and `weights`
parameters, respectively.

  - ### Neighborhood connectivity

Neighborhood connectivity is one of the other important centrality
measures that reflect the semi-local centrality of a node. This
centrality measure was first represented in a [Science
paper](https://www.ncbi.nlm.nih.gov/pubmed/11988575) in 2002 and is for
the first time calculable in R environment via the `influential`
package.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

neighrhood.co <- neighborhood.connectivity(graph = My_graph,    # Calculating neighborhood connectivity
                                           vertices = GraphVertices,
                                           mode = "all")
```

Neighborhood connectivity could be also calculated for *directed* graphs
via specifying the `mode` parameter.

  - ### H-index

H-index is H-index is another semi-local centrality measure that was
inspired from its application in assessing the impact of researchers and
is for the first time calculable in R environment via the `influential`
package.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

h.index <- h_index(graph = My_graph,    # Calculating H-index
                   vertices = GraphVertices,
                   mode = "all")
```

H-index could be also calculated for *directed* graphs via specifying
the `mode` parameter.

  - ### Local H-index

Local H-index (LH-index) is a semi-local centrality measure and an
improved version of H-index centrality that leverages the H-index to the
second order neighbors of a node and is for the first time calculable in
R environment via the `influential` package.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

lh.index <- lh_index(graph = My_graph,    # Calculating Local H-index
                   vertices = GraphVertices,
                   mode = "all")
```

Local H-index could be also calculated for *directed* graphs via
specifying the `mode` parameter.

  - ### Collective Influence

Collective Influence (CI) is a global centrality measure that calculates
the product of the reduced degree (degree - 1) of a node and the total
reduced degree of all nodes at a distance d from the node. This
centrality measure is for the first time provided in an R package.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

ci <- collective.influence(graph = My_graph,    # Calculating Collective Influence
                          vertices = GraphVertices,
                          mode = "all", d=3)
```

Collective Influence could be also calculated for *directed* graphs via
specifying the `mode` parameter.

  - ### ClusterRank

ClusterRank is a local centrality measure that makes a connection
between local and semi-local characteristics of a node and at the same
time removes the negative effects of local clustering.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

cr <- clusterRank(graph = My_graph,    # Calculating ClusterRank
                  vids = GraphVertices,
                  directed = FALSE, loops = TRUE)
```

ClusterRank could be also calculated for *directed* graphs via
specifying the `directed` parameter.

-----

## Assessment of the association of centrality measures

  - ### Conditional probability of deviation from means

The function `cond.prob.analysis` assesses the conditional probability
of deviation of two centrality measures (or any other two continuous
variables) from their corresponding means in opposite directions.

``` r
MyData <- centrality.measures        # Preparing the data

My.conditional.prob <- cond.prob.analysis(data = MyData,       # Assessing the conditional probability
                                          nodes.colname = rownames(MyData),
                                          Desired.colname = "BC",
                                          Condition.colname = "NC")

print(My.conditional.prob)
#> $ConditionalProbability
#> [1] 51.61871
#> 
#> $ConditionalProbability_split.half.sample
#> [1] 51.33333
```

  - As you can see in the results, the whole data is also randomly
    splitted into half in order to further test the validity of
    conditional probability assessments.
  - *The higher the conditional probability the more two centrality
    measures behave in contrary manners*.

-----

  - ### Nature of association (considering dependent and independent)

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
        further details please refer to [NNS: Nonlinear Nonparametric
        Statistics](https://cran.r-project.org/package=NNS): Greater
        values indicate a higher dependence between variables. <br><br>
4.  **Correlation assessment**: As the correlation between most of the
    centrality measures follows a non-monotonic form, this part of the
    assessment is done based on the non-linear non-parametric statistics
    (NNS) which itself calculates the correlation based on partial
    moments and outputs a correlation value ranging from -1 to 1. For
    further details please refer to [NNS: Nonlinear Nonparametric
    Statistics](https://cran.r-project.org/package=NNS). <br><br>
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
MyData <- centrality.measures        # Preparing the data

My.metrics.assessment <- double.cent.assess(data = MyData,       # Association assessment
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

-----

  - ### Nature of association (without considering dependence direction)

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
        further details please refer to [NNS: Nonlinear Nonparametric
        Statistics](https://cran.r-project.org/package=NNS): Greater
        values indicate a higher dependence between variables. <br><br>
3.  **Correlation assessment**: As the correlation between most of the
    centrality measures follows a non-monotonic form, this part of the
    assessment is done based on the non-linear non-parametric statistics
    (NNS) which itself calculates the correlation based on partial
    moments and outputs a correlation value ranging from -1 to 1. For
    further details please refer to [NNS: Nonlinear Nonparametric
    Statistics](https://cran.r-project.org/package=NNS). <br><br>
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
MyData <- centrality.measures        # Preparing the data

My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,       # Association assessment
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

-----

## Identification of the most `influential` network nodes

**IVI** : `IVI` is the first integrative method for the identification
of network most influential nodes in a way that captures all network
topological dimensions. The `IVI` formula integrates the most important
local (i.e. degree centrality and ClusterRank), semi-local
(i.e. neighborhood connectivity and local H-index) and global
(i.e. betweenness centrality and collective influence) centrality
measures in such a way that both synergize their effects and remove
their biases.

  - ### Integrated Value of Influence (IVI) from centrality measures

<!-- end list -->

``` r
MyData <- centrality.measures        # Preparing the data

My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,       # Calculation of IVI
                                   CR = centrality.measures$CR,
                                   NC = centrality.measures$NC,
                                   LH_index = centrality.measures$LH_index,
                                   BC = centrality.measures$BC,
                                   CI = centrality.measures$CI)
```

  - ### Integrated Value of Influence (IVI) from graph

<!-- end list -->

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

My.vertices.IVI <- ivi(graph = My_graph, vertices = GraphVertices, # Calculation of IVI
                       weights = NULL, directed = FALSE, mode = "all",
                       loops = TRUE, d = 3, scaled = TRUE)
```

IVI could be also calculated for *directed* and/or *weighted* graphs via
specifying the `directed`, `mode`, and `weights` parameters.

-----

## Identification of the most important network spreaders

**Spreading score** : `spreading.score` is an integrative score made up
of four different centrality measures including ClusterRank,
neighborhood connectivity, betweenness centrality, and collective
influence. Also, Spreading score reflects the spreading potential of
each node within a network and is one of the major components of the
`IVI`.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

Spreading.score <- spreading.score(graph = My_graph,     # Calculation of Spreading score
                                   vertices = GraphVertices, 
                                   weights = NULL, directed = FALSE, mode = "all",
                                   loops = TRUE, d = 3, scaled = TRUE)
```

Spreading score could be also calculated for *directed* and/or
*weighted* graphs via specifying the `directed`, `mode`, and `weights`
parameters.

-----

## Identification of the most important network hubs

**Hubness score** : `hubness.score` is an integrative score made up of
two different centrality measures including degree centrality and local
H-index. Also, Hubness score reflects the power of each node in its
surrounding environment and is one of the major components of the `IVI`.

``` r
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

Hubness.score <- hubness.score(graph = My_graph,     # Calculation of Hubness score
                                   vertices = GraphVertices, 
                                   directed = FALSE, mode = "all",
                                   loops = TRUE, scaled = TRUE)
```

Spreading score could be also calculated for *directed* graphs via
specifying the `directed` and `mode` parameters.

-----

## Ranking the influence of network nodes based on the `SIRIR` model

**SIRIR** : `SIRIR` is achieved by the integration
susceptible-infected-recovered (SIR) model with the leave-one-out cross
validation technique and ranks network nodes based on their true
universal influence. One of the applications of this function is the
assessment of performance of a novel algorithm in identification of
network influential nodes.

``` r
set.seed(1234)
My_graph <- igraph::sample_gnp(n=50, p=0.05)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

Influence.Ranks <- sirir(graph = My_graph,     # Calculation of influence rank
                                   vertices = GraphVertices, 
                                   beta = 0.5, gamma = 1, no.sim = 10, seed = 1234)

knitr::kable(Influence.Ranks[c(order(Influence.Ranks$rank)[1:10]),])
```

|    | difference.value | rank |
| -- | ---------------: | ---: |
| 6  |              9.0 |    1 |
| 1  |              8.9 |    2 |
| 2  |              8.9 |    2 |
| 8  |              8.9 |    2 |
| 10 |              8.7 |    5 |
| 24 |              8.7 |    5 |
| 18 |              8.6 |    7 |
| 19 |              8.6 |    7 |
| 20 |              8.6 |    7 |
| 21 |              8.6 |    7 |

-----

## Experimental-data-based classification and ranking of top candidate features based on the `ExIR` model

**ExIR** : `ExIR` is a model for the classification and ranking of top
candidate features. The input data could come from any type of
experiment such as transcriptomics and proteomics. This model is based
on multi-level filteration and scoring based on several supervised and
unsupervised analyses followed by the classification and integrative
ranking of top candidate features. Using this function and depending on
the input data and specified arguments, the user can get 1 to four
tables.

``` r

#Preparing the required arguments for the `exir` function
MyDesired_list <- Desiredlist
MyDiff_data <- Diffdata
Diff_value <- c(1,3,5)
Regr_value <- 7
Sig_value <- c(2,4,6,8)
MyExptl_data <- Exptldata
Condition_colname <- "condition"

#Running the ExIR model
My.exir <- exir(Desired_list = MyDesired_list,
Diff_data = MyDiff_data, Diff_value = Diff_value,
Regr_value = Regr_value, Sig_value = Sig_value,
Exptl_data = MyExptl_data, Condition_colname = Condition_colname,
verbose = FALSE)

names(My.exir)
#> [1] "Driver table"          "DE-mediator table"     "nonDE-mediators table"
#> [4] "Biomarker table"
```
