## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tibble.print_min = 4L, tibble.print_max = 4L)
set.seed(60)

## ----setup--------------------------------------------------------------------
library(influential)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(coexpression.data))

## ----g_dataframe--------------------------------------------------------------
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(d=MyData)        # Reconstructing the graph

## -----------------------------------------------------------------------------
class(My_graph)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(coexpression.adjacency)[,1:10])

## ----g_adj--------------------------------------------------------------------
MyData <- coexpression.adjacency        # Preparing the data

My_graph <- graph_from_adjacency_matrix(MyData)        # Reconstructing the graph

## ----Vertices-----------------------------------------------------------------
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

My_graph_vertices <- V(My_graph)        # Extracting the vertices

## ----DC-----------------------------------------------------------------------
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

My_graph_degree <- degree(My_graph, v = GraphVertices, normalized = FALSE) # Calculating degree centrality

## ----BC-----------------------------------------------------------------------
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,    # Calculating betweenness centrality
                                    directed = FALSE, normalized = FALSE)

## ----NC-----------------------------------------------------------------------
MyData <- coexpression.data        # Preparing the data

My_graph <- graph_from_data_frame(MyData)        # Reconstructing the graph

GraphVertices <- V(My_graph)        # Extracting the vertices

neighrhood.co <- neighborhood.connectivity(graph = My_graph,    # Calculating neighborhood connectivity
                                           vertices = GraphVertices,
                                           mode = "all")

## ----cond.prob----------------------------------------------------------------
MyData <- centrality.measures        # Preparing the data

My.conditional.prob <- cond.prob.analysis(data = MyData,       # Assessing the conditional probability
                                          nodes.colname = "name",
                                          Desired.colname = "BetweennessCentrality",
                                          Condition.colname = "NeighborhoodConnectivity")

print(My.conditional.prob)

## ----double.cent.assess, eval=FALSE-------------------------------------------
#  MyData <- centrality.measures        # Preparing the data
#  
#  My.metrics.assessment <- double.cent.assess(data = MyData,       # Association assessment
#                                              nodes.colname = "name",
#                                              dependent.colname = "BetweennessCentrality",
#                                              independent.colname = "NeighborhoodConnectivity")
#  
#  print(My.metrics.assessment)
#  #> $Summary_statistics
#  #>         BetweennessCentrality NeighborhoodConnectivity
#  #> Min.              0.000000000                   1.2000
#  #> 1st Qu.           0.000000000                  66.0000
#  #> Median            0.000000000                 156.0000
#  #> Mean              0.005813357                 132.3443
#  #> 3rd Qu.           0.000340000                 179.3214
#  #> Max.              0.529464720                 192.0000
#  #>
#  #> $Normality_results
#  #>                               p.value
#  #> BetweennessCentrality    1.415450e-50
#  #> NeighborhoodConnectivity 9.411737e-30
#  #>
#  #> $Dependent_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $Independent_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $GAM_nonlinear.nonmonotonic.results
#  #>      edf  p-value
#  #> 8.992406 0.000000
#  #>
#  #> $Association_type
#  #> [1] "nonlinear-nonmonotonic"
#  #>
#  #> $HoeffdingD_Statistic
#  #>         D_statistic P_value
#  #> Results  0.01770279   1e-08
#  #>
#  #> $Dependence_Significance
#  #>                       Hoeffding
#  #> Results Significantly dependent
#  #>
#  #> $NNS_dep_results
#  #>         Correlation Dependence
#  #> Results  -0.7948106  0.8647164
#  #>
#  #> $ConditionalProbability
#  #> [1] 55.35386
#  #>
#  #> $ConditionalProbability_split.half.sample
#  #> [1] 55.90331

## ----double.cent.assess.noRegr., eval=FALSE-----------------------------------
#  MyData <- centrality.measures        # Preparing the data
#  
#  My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,       # Association assessment
#                                                           nodes.colname = "name",
#                                                           centrality1.colname = "BetweennessCentrality",
#                                                           centrality2.colname = "NeighborhoodConnectivity")
#  
#  print(My.metrics.assessment)
#  #> $Summary_statistics
#  #>         BetweennessCentrality NeighborhoodConnectivity
#  #> Min.              0.000000000                   1.2000
#  #> 1st Qu.           0.000000000                  66.0000
#  #> Median            0.000000000                 156.0000
#  #> Mean              0.005813357                 132.3443
#  #> 3rd Qu.           0.000340000                 179.3214
#  #> Max.              0.529464720                 192.0000
#  #>
#  #> $Normality_results
#  #>                               p.value
#  #> BetweennessCentrality    1.415450e-50
#  #> NeighborhoodConnectivity 9.411737e-30
#  #>
#  #> $Centrality1_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $Centrality2_Normality
#  #> [1] "Non-normally distributed"
#  #>
#  #> $HoeffdingD_Statistic
#  #>         D_statistic P_value
#  #> Results  0.01770279   1e-08
#  #>
#  #> $Dependence_Significance
#  #>                       Hoeffding
#  #> Results Significantly dependent
#  #>
#  #> $NNS_dep_results
#  #>         Correlation Dependence
#  #> Results  -0.7948106  0.8647164
#  #>
#  #> $ConditionalProbability
#  #> [1] 55.35386
#  #>
#  #> $ConditionalProbability_split.half.sample
#  #> [1] 55.68163

## ----ihs----------------------------------------------------------------------
MyData <- centrality.measures        # Preparing the data

My.vertices.IHS <- ihs(DC = centrality.measures$Degree,       # Calculation of IHS
                       BC = centrality.measures$BetweennessCentrality,
                       NC = centrality.measures$NeighborhoodConnectivity)

print(head(My.vertices.IHS))

