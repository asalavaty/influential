Vignettes of the `influential`
================

## Overview

The influential package contains several functions that could be
categorized into four groups according to their purpose:

  - Network reconstruction
  - Calculation of centrality measures
  - Assessment of the association of centrality measures
  - Identification of the most `influential` network nodes (hubs)

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
package.

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
    
    Betweeneess centrality, like degree centrality, is one of the most
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

Degree centrality could be also calculated for *directed* and/or
*weighted* graphs via specifying the `directed` and `weights`
parameters, respectively.

  - ### Neighborhood connectivity

Neighborhood connectivity is one of the other important centrality
measures that reflect the (semi-) local centrality of a node. This
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

-----

## Assessment of the association of centrality measures

  - ### Conditional probability of deviation from means

The function `cond.prob.analysis` assesses the conditional probability
of deviation of two centrality measures (or any other two continuous
variables) from their corresponding means in oposite directions.

``` r
MyData <- centrality.measures        # Preparing the data

My.conditional.prob <- cond.prob.analysis(data = MyData,       # Assessing the conditional probability
                                          nodes.colname = "name",
                                          Desired.colname = "BetweennessCentrality",
                                          Condition.colname = "NeighborhoodConnectivity")

print(My.conditional.prob)
#> $ConditionalProbability
#> [1] 55.35386
#> 
#> $ConditionalProbability_split.half.sample
#> [1] 55.16751
```

  - As you can see in the results, the whole data is also randomly
    splitted into half in order to further test the validity of
    conditional probability assessments.
  - *The higher the conditional probability the more two centrality
    measures behave in contrary manners*.

-----

  - ### Nature of association (considering dependent and independent)

The function `double.cent.assess` could be used to automaticelly assess
both the distribution mode of centrality measures (two continuous
variables) and the nature of their association. The analyses done
through this formula are as follows:

1.  **Normality assessment**:
      - Variables with **lower than** 5000 observations: *Shapiro-Wilk
        test*
      - Variables with **over** 5000 observations: *Anderson-Darlin
        test* <br><br>
2.  **Assessment of non-linear/non-monotonic correlation**:
      - *Non-linearity assessment*: Fitting a generalized additive model
        (GAM) with integrated smoothness approximations using the `mgcv`
        package <br><br>
      - *Non-monotonicity assessment*: Comaring the squared coefficients
        of the correlation based on Spearman’s rank correlation analysis
        and ranked regression test with non-linear splines.
          - Squared coefficient of Spearman’s rank correlation **\>**
            R-squared ranked regression with non-linear splines:
            *Monotonic*
          - Squared coefficient of Spearman’s rank correlation **\<**
            R-squared ranked regression with non-linear splines:
            *Non-monotonic* <br><br>
3.  **Dependemnce assessment**:
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
    their corresponding means in oposite directions.
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
                                            nodes.colname = "name",
                                            dependent.colname = "BetweennessCentrality",
                                            independent.colname = "NeighborhoodConnectivity")

print(My.metrics.assessment)
#> $Summary_statistics
#>         BetweennessCentrality NeighborhoodConnectivity
#> Min.              0.000000000                   1.2000
#> 1st Qu.           0.000000000                  66.0000
#> Median            0.000000000                 156.0000
#> Mean              0.005813357                 132.3443
#> 3rd Qu.           0.000340000                 179.3214
#> Max.              0.529464720                 192.0000
#> 
#> $Normality_results
#>                               p.value
#> BetweennessCentrality    1.415450e-50
#> NeighborhoodConnectivity 9.411737e-30
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
automaticelly assess both the distribution mode of centrality measures
(two continuous variables) and the nature of their association. The
analyses done through this formula are as follows:

1.  **Normality assessment**:
      - Variables with **lower than** 5000 observations: *Shapiro-Wilk
        test*
      - Variables with **over** 5000 observations: *Anderson-Darlin
        test* <br><br>
2.  **Dependemnce assessment**:
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
    their corresponding means in oposite directions.
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
                                                         nodes.colname = "name",
                                                         centrality1.colname = "BetweennessCentrality",
                                                         centrality2.colname = "NeighborhoodConnectivity")

print(My.metrics.assessment)
#> $Summary_statistics
#>         BetweennessCentrality NeighborhoodConnectivity
#> Min.              0.000000000                   1.2000
#> 1st Qu.           0.000000000                  66.0000
#> Median            0.000000000                 156.0000
#> Mean              0.005813357                 132.3443
#> 3rd Qu.           0.000340000                 179.3214
#> Max.              0.529464720                 192.0000
#> 
#> $Normality_results
#>                               p.value
#> BetweennessCentrality    1.415450e-50
#> NeighborhoodConnectivity 9.411737e-30
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

## Identification of the most `influential` network nodes (hubs)

**IHS** : `IHS` is the first integrative method for the identification
of network hubs. The `IHS` formula integrates three network centrality
meassure including degree centrality, betweenness centrality, and
neighborhood connectivity in such a way that both synergize their
effects and remove their biases.

``` r
MyData <- centrality.measures        # Preparing the data

My.vertices.IHS <- ihs(DC = centrality.measures$Degree,       # Calculation of IHS
                       BC = centrality.measures$BetweennessCentrality,
                       NC = centrality.measures$NeighborhoodConnectivity)

print(head(My.vertices.IHS))
#> [1] 196.2039215   2.9822273   0.1572078   6.1757221   0.3199993   0.5095222
```
