#=============================================================================
#
#    Code chunk 0: Influential package description
#
#=============================================================================

#' @keywords internal
#' @title Influential package
#' @description
#' The goal of \emph{\strong{\code{influential}}} is to help identification of the most influential nodes in a network
#' as well as the classification and ranking of top candidate features.
#' This package contains functions for the classification and ranking of features,
#' reconstruction of networks from adjacency matrices and
#' data frames, analysis of the topology of the network and calculation of centrality measures
#' as well as a novel and powerful influential node ranking.
#' The \strong{Experimental data-based Integrative Ranking (ExIR)} is a sophisticated model
#' for classification and ranking of the top candidate features based on only the experimental data.
#' The first integrative method, namely the \strong{Integrated Value of Influence (IVI)},
#' that captures all topological dimensions of the network for
#' the identification of network most influential nodes is also provided as
#' a function. Also, neighborhood connectivity, H-index, local H-index, and collective
#' influence (CI), all of which required centrality measures for the calculation of IVI,
#' are for the first time provided in an R package. Additionally, a function is provided
#' for running \strong{SIRIR} model, which is the combination of leave-one-out cross validation
#' technique and the conventional SIR model, on a network to unsupervisedly rank the true
#' influence of vertices.Furthermore, some functions have been provided for the
#' assessment of dependence and correlation of two network centrality measures as well
#' as the conditional probability of deviation from their corresponding
#' means in opposite directions.
#'
#' You may check the latest developmental version of the \emph{influential} package on its
#' \href{https://github.com/asalavaty/influential}{GitHub repository}
#'
#' @details
#' \itemize{
#'   \item Package: influential
#'   \item Type: Package
#'   \item Version: 1.1.1
#'   \item Date: 24-06-2020
#'   \item License: GPL-3
#' }
#'
#' @author
#' Author: Adrian (Abbas) Salavaty
#'
#' Advisors: Mirana Ramialison and Peter D. Currie
#'
#'
#' Maintainer: Adrian (Abbas) Salavaty \email{abbas.salavaty@@gmail.com}
#'
#'
#' You may find more information on my personal website at \href{https://www.abbassalavaty.com/}{www.AbbasSalavaty.com}
#'
#' @references
#' \itemize{
#'   \item Fred Viole and David Nawrocki (2013, ISBN:1490523995).
#'   \item Csardi G, Nepusz T (2006). “The igraph software package for complex network research.”
#' InterJournal, Complex Systems, 1695. \url{http://igraph.org}.
#' }
#'
#' \strong{Note:} Adopted algorithms and sources are referenced in function document.
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

#=============================================================================
#
#    Code chunk 1: Calculation of neighborhood connectivity
#
#=============================================================================

#' Neighborhood connectivity
#'
#' This function calculates the neighborhood connectivity of input vertices and
#' works with both directed and undirected networks.
#' @param graph A graph (network) of the igraph class.
#' @param vertices A vector of desired vertices, which could be obtained by the V function.
#' @param mode The mode of neighborhood connectivity depending on the directedness of the graph.
#' If the graph is undirected, the mode "all" should be specified.
#' Otherwise, for the calculation of neighborhood connectivity based on
#' incoming connections select "in" and for the outgoing connections select "out".
#' Also, if all of the connections are desired, specify the "all" mode. Default mode is set to "all".
#' @return A vector including the neighborhood connectivity score of each vertex inputted.
#' @aliases NC
#' @keywords neighborhood_connectivity
#' @family centrality functions
#' @export
#' @examples
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' neighrhood.co <- neighborhood.connectivity(graph = My_graph,
#'                                            vertices = GraphVertices,
#'                                            mode = "all")
neighborhood.connectivity <- function(graph, vertices = V(graph), mode = "all") {

  # Getting the names of vertices with the order in use
  node.names <- as.character(igraph::as_ids(igraph::V(graph = graph)))

  # Getting the first neighbors of each node
  node.neighbors <- sapply(as.list(node.names),
                           FUN = function(i) as.character(igraph::as_ids(igraph::neighbors(graph = graph,
                                                                                           v = i,
                                                                                           mode = mode))))

  # Getting the neighborhood size of each node
    first.neighbors.size <- sapply(node.neighbors,
                                   function(s) igraph::neighborhood.size(graph = graph,
                                                                         nodes = s,
                                                                         mode = mode,
                                                                         order = 1) - 1)

  first.neighbors.size.sum <- sapply(first.neighbors.size, sum)

  # Calculation of neighborhood connectivity
  temp.nc <- vector(mode = "numeric", length = length(vertices))

  for (i in 1:length(vertices)) {
    temp.nc[i] <- first.neighbors.size.sum[i]/length(node.neighbors[[i]])
  }

  nc.table <- data.frame(Neighborhood_connectivity = temp.nc)

  if (length(vertices) > 1) {
    nc.table <- nc.table
  } else if (length(vertices) == 1) {
    nc.table <- data.frame(Neighborhood_connectivity =
                             sum(nc.table$Neighborhood_connectivity)/
                             nrow(nc.table))
  }

  rownames(nc.table) <- node.names

  nc.table$Neighborhood_connectivity[c(which(is.nan(nc.table$Neighborhood_connectivity)),
                                       which(is.na(nc.table$Neighborhood_connectivity)))] <- 0

  nc.table <- nc.table[,1]

  return(nc.table)

}

#=============================================================================
#
#    Code chunk 2: Calculation of H-index
#
#=============================================================================

#' H-index
#'
#' This function calculates the H-index of input vertices and
#' works with both directed and undirected networks.
#' @param graph A graph (network) of the igraph class.
#' @param vertices A vector of desired vertices, which could be obtained by the V function.
#' @param mode The mode of H-index depending on the directedness of the graph.
#' If the graph is undirected, the mode "all" should be specified.
#' Otherwise, for the calculation of H-index based on
#' incoming connections select "in" and for the outgoing connections select "out".
#' Also, if all of the connections are desired, specify the "all" mode. Default mode is set to "all".
#' @return A vector including the H-index of each vertex inputted.
#' @aliases h.index
#' @keywords h_index
#' @family centrality functions
#' @seealso \code{\link[influential]{lh_index}}
#' @export
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' h.index <- h_index(graph = My_graph, vertices = GraphVertices, mode = "all")
#' }
#' @importFrom utils tail
  h_index <- function(graph, vertices = V(graph), mode = "all") {

  # Getting the first neighbors of each node
  first.neighbors <- igraph::neighborhood(graph, nodes = vertices, mode = mode)

  # Getting the neighbors of each vertex
  node.neighbors <- sapply(first.neighbors, function(n) rownames(as.matrix(n[[]][-1])), simplify = F)

  # Getting the neighborhood size of each node
  first.neighbors.size <- lapply(node.neighbors, function(s) igraph::neighborhood.size(graph, s,
                                                                               mode = mode, order = 1) - 1)

  # Calculation of H-index
  hindex <- vector(mode = "integer", length = length(vertices))

  for (i in 1:length(vertices)) {

    temp.neighbors.size <- unlist(first.neighbors.size[i])

    temp.neighbors.size <- temp.neighbors.size[order(temp.neighbors.size,
                                                     decreasing = TRUE)]
    if(length(temp.neighbors.size) == 0) { hindex[i] <- 0

    } else if(max(temp.neighbors.size) == 0) {
      hindex[i] <- 0
    } else {hindex[i] <- utils::tail(which(temp.neighbors.size >=
                                      seq_along(temp.neighbors.size)), 1)
    }

    rm(temp.neighbors.size)
  }
  return(hindex)
  }

  #=============================================================================
  #
  #    Code chunk 3: Calculation of local H-index (LH-index)
  #
  #=============================================================================

  #' local H-index (LH-index)
  #'
  #' This function calculates the local H-index of input vertices and
  #' works with both directed and undirected networks.
  #' @param graph A graph (network) of the igraph class.
  #' @param vertices A vector of desired vertices, which could be obtained by the V function.
  #' @param mode The mode of local H-index depending on the directedness of the graph.
  #' If the graph is undirected, the mode "all" should be specified.
  #' Otherwise, for the calculation of local H-index based on
  #' incoming connections select "in" and for the outgoing connections select "out".
  #' Also, if all of the connections are desired, specify the "all" mode. Default mode is set to "all".
  #' @return A vector including the local H-index of each vertex inputted.
  #' @aliases lh.index
  #' @keywords lh_index
  #' @family centrality functions
  #' @export
  #' @examples
  #' \dontrun{
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' lh.index <- lh_index(graph = My_graph, vertices = GraphVertices, mode = "all")
  #' }
  lh_index <- function(graph, vertices = V(graph), mode = "all") {

    # Getting the first neighbors of each node
    first.neighbors <- igraph::neighborhood(graph, nodes = vertices, mode = mode)

    # Calculation of local H-index (LH-index)
    lhindex <- vector(mode = "integer", length = length(vertices))

    for (i in 1:length(vertices)) {

      lhindex[i] <- sum(h_index(graph = graph,
                                vertices = unlist(first.neighbors[i]),
                                mode = mode))
    }
    return(lhindex)
  }

  #=============================================================================
  #
  #    Code chunk 4: Calculation of Collective Influence (CI)
  #
  #=============================================================================

  #' Collective Influence (CI)
  #'
  #' This function calculates the collective influence of input vertices and
  #' works with both directed and undirected networks. This function and its descriptions are
  #' obtained from https://github.com/ronammar/collective_influence with minor modifications.
  #' Collective Influence as described by Morone & Makse (2015). In simple terms,
  #' it is the product of the reduced degree (degree - 1) of a node and the total reduced
  #' degree of all nodes at a distance d from the node.
  #' @param graph A graph (network) of the igraph class.
  #' @param vertices A vector of desired vertices, which could be obtained by the V function.
  #' @param mode The mode of collective influence depending on the directedness of the graph.
  #' If the graph is undirected, the mode "all" should be specified.
  #' Otherwise, for the calculation of collective influence based on
  #' incoming connections select "in" and for the outgoing connections select "out".
  #' Also, if all of the connections are desired, specify the "all" mode. Default mode is set to "all".
  #' @param d The distance, expressed in number of steps from a given node (default=3). Distance
  #' must be > 0. According to Morone & Makse (https://doi.org/10.1038/nature14604), optimal
  #' results can be reached at d=3,4, but this depends on the size/"radius" of the network.
  #' NOTE: the distance d is not inclusive. This means that nodes at a distance of 3 from
  #' our node-of-interest do not include nodes at distances 1 and 2. Only 3.
  #' @return A vector of collective influence for each vertex of the graph corresponding to
  #' the order of vertices output by V(graph).
  #' @aliases CI
  #' @keywords collective.influence
  #' @family centrality functions
  #' @export
  #' @examples
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' ci <- collective.influence(graph = My_graph, vertices = GraphVertices, mode = "all", d=3)
  collective.influence <- function(graph, vertices = V(graph), mode = "all", d=3) {

    ci <- vector(mode="numeric", length=length(vertices))  # collective influence output

    reduced.degrees <- degree(graph = graph,
                              v = vertices,
                              mode = mode) - 1

    # Only identify nodes at distance d
    nodes.at.distance <- igraph::neighborhood(graph = graph, nodes = vertices,
                                              mode = mode, order=d, mindist=d)

    for (i in 1:length(nodes.at.distance)) {
      rd <- reduced.degrees[i]  # i is the index of the node
      rd.list <- reduced.degrees[igraph::as_ids(nodes.at.distance[[i]])]
      # Setting 0 as default in case the graph doesn't have 2nd-order neighbours
      rd.neighbours <- ifelse(length(rd.list) > 0, sum(rd.list), 0)
      ci[i] <- rd * rd.neighbours
    }

    return(ci)
  }

  #=============================================================================
  #
  #    Code chunk 5: Calculation of ClusterRank
  #
  #=============================================================================

  #' ClusterRank (CR)
  #'
  #' This function calculates the ClusterRank of input vertices and
  #' works with both directed and undirected networks.
  #' This function and all of its descriptions have been adapted from the centiserve package with
  #' some minor modifications. ClusterRank is a local ranking algorithm which takes into account not only
  #' the number of neighbors and the neighbors’ influences, but also the clustering coefficient.
  #' @param graph The input graph as igraph object
  #' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices.
  #' @param directed Logical scalar, whether to directed graph is analyzed. This argument is ignored for undirected graphs.
  #' @param loops Logical; whether the loop edges are also counted.
  #' @return A numeric vector contaning the ClusterRank centrality scores for the selected vertices.
  #' @aliases CR
  #' @keywords clusterRank
  #' @family centrality functions
  #' @export
  #' @examples
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' cr <- clusterRank(graph = My_graph, vids = GraphVertices, directed = FALSE, loops = TRUE)
  clusterRank <- function(graph, vids = V(graph),
                        directed = FALSE, loops = TRUE) {

  vertex.transitivity <- vector(mode = "numeric")

  if(directed) {

    cl.Rank.mode <- "out"

    for(i in V(graph)) {
      vertex.neighborhood <- igraph::neighborhood(graph = graph,
                                          order = 1, nodes=i,
                                          mode = cl.Rank.mode)[[1]][-1]
      if(length(vertex.neighborhood) < 2){
        vertex.transitivity <- base::append(vertex.transitivity, NaN)
      } else {
        indc.subgraph <- igraph::induced.subgraph(graph = graph, vertex.neighborhood)
        vertex.transitivity <- base::append(vertex.transitivity,
                                            igraph::ecount(indc.subgraph)/(igraph::vcount(indc.subgraph)*(igraph::vcount(indc.subgraph)-1)))
      }
    }

  } else {

    cl.Rank.mode <- "all"

    vertex.transitivity <- igraph::transitivity(graph = graph, type = "local")

  }

  if(class(vids) == "igraph.vs") {
    vertices.index <- stats::na.omit(match(vids, V(graph)))
  } else {
    vertices.index <- stats::na.omit(match(vids, igraph::as_ids(V(graph))))
  }

  cl.Rank <- vector(mode = "numeric")

  for (i in V(graph)[vertices.index]) {
    if (is.nan(vertex.transitivity[i])) {
      cl.Rank <- append(cl.Rank, NaN)
    }
    else {

      selected.v.neighborhood <- igraph::neighborhood(graph = graph,
                                                      order = 1, nodes = i,
                                                      mode = cl.Rank.mode)[[1]][-1]
      temp.cl.Rank <- 0
      for (j in selected.v.neighborhood) {
        temp.cl.Rank <- temp.cl.Rank + igraph::degree(graph = graph,
                                                      v = j, mode = cl.Rank.mode,
                                                      loops = loops) + 1
      }
      cl.Rank <- append(cl.Rank, temp.cl.Rank * vertex.transitivity[i])
    }
  }

  if (igraph::is.named(graph)) {
    names(cl.Rank) <- igraph::V(graph)$name[vertices.index]
  }

  return(cl.Rank)
  }

#=============================================================================
#
#    Code chunk 6: Calculation of conditional probability of deviation from
#                  means in opposite directions
#
#=============================================================================

#' Conditional probability of deviation from means
#'
#' This function calculates the conditional probability of deviation of two
#' centrality measures (or any two other continuous variables) from their corresponding
#' means in opposite directions.
#' @param data A data frame containing the values of two continuous variables and the name of
#' observations (nodes).
#' @param nodes.colname The character format (quoted) name of the column containing
#' the name of observations (nodes).
#' @param Desired.colname The character format (quoted) name of the column containing
#' the values of the desired variable.
#' @param Condition.colname The character format (quoted) name of the column containing
#' the values of the condition variable.
#' @return A list of two objects including the conditional probability of deviation of two
#' centrality measures (or any two other continuous variables) from their corresponding
#' means in opposite directions based on both the entire network and the split-half random
#' sample of network nodes.
#' @aliases CPA
#' @keywords conditional_probability association_assessment
#' @family centrality association assessment functions
#' @export
#' @examples
#' MyData <- centrality.measures
#' My.conditional.prob <- cond.prob.analysis(data = MyData,
#'                                           nodes.colname = rownames(MyData),
#'                                           Desired.colname = "BC",
#'                                           Condition.colname = "NC")
cond.prob.analysis <- function(data, nodes.colname, Desired.colname, Condition.colname) {

  #filtering the data to find those nodes meeting the conditions
  ncpositive <- data[data[, Condition.colname] >
                       mean(data[, Condition.colname]),]
  ncpositive.bcnegative <- ncpositive[ncpositive[, Desired.colname] <
                                        mean(data[, Desired.colname]),]
  ncpositive.bcnegative.prob <- (nrow(ncpositive.bcnegative)/nrow(ncpositive))*100

  ncnegative <- data[data[, Condition.colname] <
                       mean(data[, Condition.colname]),]
  ncnegative.bcpositive <- ncnegative[ncnegative[, Desired.colname] >
                                        mean(data[, Desired.colname]),]
  ncnegative.bcpositive.prob <- (nrow(ncnegative.bcpositive)/nrow(ncnegative))*100

  #calculation of conditional probability
  final.cond.prob <- sum(ncpositive.bcnegative.prob, ncnegative.bcpositive.prob)/2

  ##Reliability analysis based on split-half random sampling

  #split-half random sampling

  sample.data <- data[sample(1:nrow(data), size = round(nrow(data)/2), replace = FALSE),]

  #filtering the data to find those nodes meeting the conditions
  sample.ncpositive <- sample.data[sample.data[, Condition.colname] >
                                     mean(sample.data[, Condition.colname]),]
  sample.ncpositive.bcnegative <- sample.ncpositive[sample.ncpositive[, Desired.colname] <
                                                      mean(sample.data[, Desired.colname]),]
  sample.ncpositive.bcnegative.prob <- (nrow(sample.ncpositive.bcnegative)/nrow(sample.ncpositive))*100

  sample.ncnegative <- sample.data[sample.data[, Condition.colname] <
                                     mean(sample.data[, Condition.colname]),]
  sample.ncnegative.bcpositive <- sample.ncnegative[sample.ncnegative[, Desired.colname] >
                                                      mean(sample.data[, Desired.colname]),]
  sample.ncnegative.bcpositive.prob <- (nrow(sample.ncnegative.bcpositive)/nrow(sample.ncnegative))*100

  #calculation of conditional probability
  sample.final.cond.prob <- sum(sample.ncpositive.bcnegative.prob, sample.ncnegative.bcpositive.prob)/2

  results <- list(ConditionalProbability = final.cond.prob,
                  ConditionalProbability_split.half.sample = sample.final.cond.prob)

  return(results)
}


#=============================================================================
#
#    Code chunk 7: Assessment of innate features and associations of two network
#                  centrality measures, one independent and one dependent
#
#=============================================================================

#' Assessment of innate features and associations of two network centrality measures (dependent and independent)
#'
#' This function assesses innate features and the association of two centrality measures
#' (or any two other continuous variables) from the aspect of distribution mode, dependence,
#' linearity, monotonicity, partial-moments based correlation, and conditional probability of
#' deviating from corresponding means in opposite direction. This function assumes one
#' variable as dependent and the other as independent for regression analyses. The non-linear nature of
#' the association of two centrality measures is evaluated based on generalized additive models (GAM).
#' The monotonicity of the association is evaluated based on comparing the squared coefficient of
#' Spearman correlation and R-squared of rank regression analysis.
#' Also, the correlation between two variables is assessed via non-linear non-parametric statistics (NNS).
#' For the conditional probability assessment, the independent variable is considered as the condition variable.
#' @param data A data frame containing the values of two continuous variables and the name of
#' observations (nodes).
#' @param nodes.colname The character format (quoted) name of the column containing
#' the name of observations (nodes).
#' @param dependent.colname The character format (quoted) name of the column containing
#' the values of the dependent variable.
#' @param independent.colname The character format (quoted) name of the column containing
#' the values of the independent variable.
#' @param plot logical; FALSE (default) Plots quadrant means of NNS correlation analysis.
#' @return A list of 11 objects including:
#'
#'     - Summary of the basic statistics of two centrality measures (or any two other continuous variables).
#'
#'     - The results of normality assessment of two variable (p-value > 0.05 imply that the variable is normally distributed).
#'
#'     - Description of the normality assessment of the dependent variable.
#'
#'     - Description of the normality assessment of the independent variable.
#'
#'     - Results of the generalized additive modeling (GAM) of the data.
#'
#'     - The association type based on simultaneous consideration of normality assessment,
#' GAM Computation with smoothness estimation, Spearman correlation, and ranked regression analysis of splines.
#'
#'     - The Hoeffding's D Statistic of dependence (ranging from -0.5 to 1).
#'
#'     - Description of the dependence significance.
#'
#'     - Correlation between variables based on the NNS method.
#'
#'     - The last two objects are the conditional probability of deviation of two
#' centrality measures from their corresponding means in opposite directions based
#' on both the entire network and the split-half random sample of network nodes.
#' @aliases DCA
#' @keywords association_assessment dependence_assessment
#' @family centrality association assessment functions
#' @seealso \code{\link[nortest]{ad.test}} for Anderson-Darling test for normality,
#' \code{\link[mgcv]{gam}} for Generalized additive models with integrated smoothness estimation,
#' \code{\link[stats]{lm}} for Fitting Linear Models,
#' \code{\link[Hmisc]{hoeffd}} for Matrix of Hoeffding's D Statistics, and
#' \code{\link[NNS]{NNS.dep}} for NNS Dependence
#' @export
#' @examples
#' \dontrun{
#' MyData <- centrality.measures
#' My.metrics.assessment <- double.cent.assess(data = MyData,
#'                                             nodes.colname = rownames(MyData),
#'                                             dependent.colname = "BC",
#'                                             independent.colname = "NC")
#' }
double.cent.assess <- function(data, nodes.colname, dependent.colname, independent.colname, plot = FALSE) {

  if("parallel" %in% (.packages())) {
    detach("package:parallel", unload = TRUE)
    base::attachNamespace("parallel")
    parallel::detectCores(logical = TRUE)
  } else {
    base::attachNamespace("parallel")
    parallel::detectCores(logical = TRUE)
  }

  #Checking the availability of required packages

  if (nrow(data) >= 5000) { if(!requireNamespace(c("nortest", "Hmisc", "mgcv", "NNS"), quietly = TRUE)) {
    stop("The packages \"nortest\" \"Hmisc\", \"mgcv\" and \"NNS\" are required for this function to work.
    Please install the required packages before using this function.

  You can install the packages via one of the following options:

         install.packages(\"Package Name\")

         Or

         install.packages(\"BiocManager\")
         BiocManager::install(\"Package Name\")",
         call. = FALSE)
  }
  }

  if(!requireNamespace(c("Hmisc", "mgcv", "NNS"), quietly = TRUE)) {
    stop("The packages \"Hmisc\", \"mgcv\" and \"NNS\" are required for this function to work.
    Please install the required packages before using this function.

  You can install the packages via one of the following options:

         install.packages(\"Package Name\")

         Or

         install.packages(\"BiocManager\")
         BiocManager::install(\"Package Name\")",
         call. = FALSE)
  }

  #checking the normality of data
  summary.stat <- apply(data[, c(dependent.colname, independent.colname)], 2, summary)

  if(length(unique(data[,dependent.colname])) < 3 &
     length(unique(data[,independent.colname])) >= 3) {
    if(nrow(data) < 5000) {
      normality <- data.frame(p.value = c(NA, stats::shapiro.test(data[, independent.colname])$p.value))
    } else if (nrow(data) >= 5000) {
      normality <- data.frame(p.value = c(NA, nortest::ad.test(data[, independent.colname])$p.value))
    }
  } else if (length(unique(data[,dependent.colname])) >= 3 &
             length(unique(data[,independent.colname])) < 3) {
    if(nrow(data) < 5000) {
      normality <- data.frame(p.value = c(stats::shapiro.test(data[, dependent.colname])$p.value, NA))
    } else if (nrow(data) >= 5000) {
      normality <- data.frame(p.value = c(nortest::ad.test(data[, dependent.colname])$p.value, NA))
    }
  } else if(length(unique(data[,dependent.colname])) < 3 &
            length(unique(data[,independent.colname])) < 3) {
    normality <- data.frame(p.value = c(NA, NA))
  } else {

    if(nrow(data) < 5000) {
      normality <- apply(data[, c(dependent.colname, independent.colname)], 2, stats::shapiro.test)
    } else if(nrow(data) >= 5000) {
      normality <- apply(data[, c(dependent.colname, independent.colname)], 2, nortest::ad.test)
    }
    normality <- as.data.frame(sapply(normality, function(m) m[]$p.value))
    colnames(normality) <- "p.value"
  }

  if(is.na(normality[1,1])) {dependent.normality <- NA} else if (normality[1,1] < 0.05) {
    dependent.normality <- "Non-normally distributed"
  } else {dependent.normality <- "Normally distributed"}

  if(is.na(normality[2,1])) {independent.normality <- NA} else if(normality[2,1] < 0.05) {
    independent.normality <- "Non-normally distributed"
  } else {independent.normality <- "Normally distributed"}

  #Assessment of non-linear/non-monotonic correlation of dependent and independent variables
  nl.assess <- summary(mgcv::gam(data[, dependent.colname] ~ s(data[, independent.colname])))
  nl.assess <- nl.assess$s.table[,c(1,4)]

  #Assessment of non-monotonic vs non-linear monotonic correlation of dependent and independent variables
  squared.pearson <- stats::cor(data[, dependent.colname], data[, independent.colname])^2
  squared.spearman <- stats::cor(rank(data[, dependent.colname]), rank(data[, independent.colname]))^2
  squared.regression <- summary(stats::lm(rank(data[, dependent.colname]) ~
                                     splines::ns(rank(data[, independent.colname]),
                                        df = 6)))$r.squared
  if(nl.assess[1] > 1 & squared.spearman < squared.regression) {
    association.type <- "nonlinear-nonmonotonic"
  } else if(nl.assess[1] > 1 & squared.spearman > squared.regression) {
    association.type <- "nonlinear-monotonic"
  } else if(nl.assess[1] <= 1 & squared.spearman < squared.pearson) {
    association.type <- "linear-monotonic"
  }


  #calculation of Hoeffding’s D Statistics (Hoeffding Dependence Coefficient)
  hoeffd <- data.frame(D_statistic = as.data.frame(Hmisc::hoeffd(x = data[, independent.colname],
                                                          y = data[, dependent.colname])[1])[1,2],
                       P_value = as.data.frame(Hmisc::hoeffd(x = data[, independent.colname],
                                                      y = data[, dependent.colname])[3])[1,2], row.names = "Results")

  if(hoeffd[1,2] < 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Significantly dependent",
                                          row.names = "Results")
  } else if(hoeffd[1,2] >= 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Not significantly dependent",
                                          row.names = "Results")
  } else if(hoeffd[1,2] < 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Significantly dependent",
                                          row.names = "Results")
  } else if (hoeffd[1,2] >= 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Not significantly dependent",
                                          row.names = "Results")
  }


  ##assessment of descriptive non-linear non-parametric correlation/dependence between
  #dependent and independent variables
  if(association.type == "nonlinear-nonmonotonic") {
    if(plot == TRUE) {
    #prepare a PDF devide to save the plot in
    grDevices::pdf(file = paste("NNS_scatter.plot", "pdf", sep = "."),
        width = 26, height = 12) }
    nl.cor.dep <- NNS::NNS.dep(x = data[, independent.colname], y = data[, dependent.colname],
                          print.map = plot, order = 2)
    if(plot == TRUE) {
    grDevices::dev.off() }
    nl.cor.dep <- data.frame(Correlation = unlist(nl.cor.dep)[1],
                             Dependence = unlist(nl.cor.dep)[2], row.names = "Results")
  } else if(association.type == "nonlinear-monotonic") {
    if(plot == TRUE) {
    #prepare a PDF devide to save the plot in
    grDevices::pdf(file = paste("NNS_scatter.plot", "pdf", sep = "."),
        width = 26, height = 12) }
    nl.cor.dep <- NNS::NNS.dep(x = data[, independent.colname], y = data[, dependent.colname],
                          print.map = plot, order = 2)
    if(plot == TRUE) {
    grDevices::dev.off() }
    nl.cor.dep <- data.frame(Correlation = unlist(nl.cor.dep)[1],
                             Dependence = unlist(nl.cor.dep)[2], row.names = "Results")
  } else {nl.cor.dep <- "The association is linear!"}

  ##assessment of conditional probability of deviation of BC and NC from their means in opposite directions
  #filtering the data to find those nodes meeting the conditions
  ncpositive <- data[data[, independent.colname] >
                       mean(data[, independent.colname]),]
  ncpositive.bcnegative <- ncpositive[ncpositive[, dependent.colname] <
                                        mean(data[, dependent.colname]),]
  ncpositive.bcnegative.prob <- (nrow(ncpositive.bcnegative)/nrow(ncpositive))*100

  ncnegative <- data[data[, independent.colname] <
                       mean(data[, independent.colname]),]
  ncnegative.bcpositive <- ncnegative[ncnegative[, dependent.colname] >
                                        mean(data[, dependent.colname]),]
  ncnegative.bcpositive.prob <- (nrow(ncnegative.bcpositive)/nrow(ncnegative))*100

  #calculation of conditional probability
  final.cond.prob <- sum(ncpositive.bcnegative.prob, ncnegative.bcpositive.prob)/2

  ##Reliability analysis based on split-half random sampling

  #split-half random sampling

  sample.data <- data[sample(1:nrow(data), size = round(nrow(data)/2), replace = FALSE),]

  #filtering the data to find those nodes meeting the conditions
  sample.ncpositive <- sample.data[sample.data[, independent.colname] >
                                     mean(sample.data[, independent.colname]),]
  sample.ncpositive.bcnegative <- sample.ncpositive[sample.ncpositive[, dependent.colname] <
                                                      mean(sample.data[, dependent.colname]),]
  sample.ncpositive.bcnegative.prob <- (nrow(sample.ncpositive.bcnegative)/nrow(sample.ncpositive))*100

  sample.ncnegative <- sample.data[sample.data[, independent.colname] <
                                     mean(sample.data[, independent.colname]),]
  sample.ncnegative.bcpositive <- sample.ncnegative[sample.ncnegative[, dependent.colname] >
                                                      mean(sample.data[, dependent.colname]),]
  sample.ncnegative.bcpositive.prob <- (nrow(sample.ncnegative.bcpositive)/nrow(sample.ncnegative))*100

  #calculation of conditional probability
  sample.final.cond.prob <- sum(sample.ncpositive.bcnegative.prob, sample.ncnegative.bcpositive.prob)/2

  results <- list(Summary_statistics = summary.stat,
                  Normality_results = normality,
                  Dependent_Normality = dependent.normality,
                  Independent_Normality = independent.normality,
                  GAM_nonlinear.nonmonotonic.results = nl.assess,
                  Association_type = association.type,
                  HoeffdingD_Statistic = hoeffd,
                  Dependence_Significance = dependence.significance,
                  NNS_dep_results = nl.cor.dep,
                  ConditionalProbability = final.cond.prob,
                  ConditionalProbability_split.half.sample = sample.final.cond.prob)

  return(results)

}

#=============================================================================
#
#    Code chunk 8: Assessment of innate features and associations of two network centrality
#                  measures, without considering dependent and independent ones
#
#=============================================================================


#' Assessment of innate features and associations of two network centrality measures
#'
#' This function assesses innate features and the association of two centrality measures
#' (or any two other continuous variables) from the aspect of distribution mode, dependence,
#' linearity, partial-moments based correlation, and conditional probability of
#' deviating from corresponding means in opposite direction (centrality2 is used as the condition variable).
#' This function doesn't consider which variable is dependent and which one is
#' independent and no regression analysis is done.
#' Also, the correlation between two variables is assessed via non-linear non-parametric statistics (NNS).
#' For the conditional probability assessment, the centrality2 variable is considered
#' as the condition variable.
#' @param data A data frame containing the values of two continuous variables and the name of
#' observations (nodes).
#' @param nodes.colname The character format (quoted) name of the column containing
#' the name of observations (nodes).
#' @param centrality1.colname The character format (quoted) name of the column containing
#' the values of the Centrality_1 variable.
#' @param centrality2.colname The character format (quoted) name of the column containing
#' the values of the Centrality_2 variable.
#' @return A list of nine objects including:
#'
#' - Summary of the basic statistics of two centrality measures (or any two other continuous variables).
#'
#' - The results of normality assessment of two variable (p-value > 0.05 imply that the variable is normally distributed).
#'
#' - Description of the normality assessment of the centrality1 (first variable).
#'
#' - Description of the normality assessment of the centrality2 (second variable).
#'
#' - The Hoeffding's D Statistic of dependence (ranging from -0.5 to 1).
#'
#' - Description of the dependence significance.
#'
#' - Correlation between variables based on the NNS method.
#'
#' - The last two objects are the conditional probability of deviation of two
#' centrality measures from their corresponding means in opposite directions based
#' on both the entire network and the split-half random sample of network nodes.
#' @aliases DCANR
#' @keywords association_assessment dependence_assessment
#' @family centrality association assessment functions
#' @seealso \code{\link[nortest]{ad.test}} for Anderson-Darling test for normality,
#' \code{\link[Hmisc]{hoeffd}} for Matrix of Hoeffding's D Statistics, and
#' \code{\link[NNS]{NNS.dep}} for NNS Dependence
#' @export
#' @examples
#' \dontrun{
#' MyData <- centrality.measures
#' My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,
#'                                             nodes.colname = rownames(MyData),
#'                                             centrality1.colname = "BC",
#'                                             centrality2.colname = "NC")
#' }
double.cent.assess.noRegression <- function(data, nodes.colname,
                                            centrality1.colname,
                                            centrality2.colname) {

  if("parallel" %in% (.packages())) {
    detach("package:parallel", unload = TRUE)
    base::attachNamespace("parallel")
    parallel::detectCores(logical = TRUE)
  } else {
    base::attachNamespace("parallel")
    parallel::detectCores(logical = TRUE)
  }

  #Checking the availability of required packages

  if (nrow(data) >= 5000) { if(!requireNamespace(c("nortest", "Hmisc", "NNS"), quietly = TRUE)) {
    stop("The packages \"nortest\", \"Hmisc\" and \"NNS\" are required for this function to work.
    Please install the required packages before using this function.

  You can install the packages via one of the following options:

         install.packages(\"Package Name\")

         Or

         install.packages(\"BiocManager\")
         BiocManager::install(\"Package Name\")",
         call. = FALSE)
  }
  }

  if(!requireNamespace(c("Hmisc", "mgcv", "NNS"), quietly = TRUE)) {
    stop("The packages \"Hmisc\", \"mgcv\" and \"NNS\" are required for this function to work.
    Please install the required packages before using this function.

  You can install the packages via one of the following options:

         install.packages(\"Package Name\")

         Or

         install.packages(\"BiocManager\")
         BiocManager::install(\"Package Name\")",
         call. = FALSE)
  }

  #checking the normality of data
  summary.stat <- apply(data[, c(centrality1.colname, centrality2.colname)], 2, summary)

  if(length(unique(data[,centrality1.colname])) < 3 &
     length(unique(data[,centrality2.colname])) >= 3) {
    if(nrow(data) < 5000) {
      normality <- data.frame(p.value = c(NA, stats::shapiro.test(data[, centrality2.colname])$p.value))
    } else if (nrow(data) >= 5000) {
      normality <- data.frame(p.value = c(NA, nortest::ad.test(data[, centrality2.colname])$p.value))
    }
  } else if (length(unique(data[,centrality1.colname])) >= 3 &
             length(unique(data[,centrality2.colname])) < 3) {
    if(nrow(data) < 5000) {
      normality <- data.frame(p.value = c(stats::shapiro.test(data[, centrality1.colname])$p.value, NA))
    } else if (nrow(data) >= 5000) {
      normality <- data.frame(p.value = c(nortest::ad.test(data[, centrality1.colname])$p.value, NA))
    }
  } else if(length(unique(data[,centrality1.colname])) < 3 &
            length(unique(data[,centrality2.colname])) < 3) {
    normality <- data.frame(p.value = c(NA, NA))
  } else {
    if(nrow(data) < 5000) {
      normality <- apply(data[, c(centrality1.colname, centrality2.colname)], 2, stats::shapiro.test)
    } else if(nrow(data) >= 5000) {
      normality <- apply(data[, c(centrality1.colname, centrality2.colname)], 2, nortest::ad.test)
    }
    normality <- as.data.frame(sapply(normality, function(m) m[]$p.value))
    colnames(normality) <- "p.value"
  }

  if(is.na(normality[1,1])) {dependent.normality <- NA} else if(normality[1,1] < 0.05) {
    dependent.normality <- "Non-normally distributed"
  } else {dependent.normality <- "Normally distributed"}

  if(is.na(normality[2,1])) {independent.normality <- NA} else if(normality[2,1] < 0.05) {
    independent.normality <- "Non-normally distributed"
  } else {independent.normality <- "Normally distributed"}

  #calculation of Hoeffding’s D Statistics (Hoeffding Dependence Coefficient)
  hoeffd <- data.frame(D_statistic = as.data.frame(Hmisc::hoeffd(x = data[, centrality2.colname],
                                                          y = data[, centrality1.colname])[1])[1,2],
                       P_value = as.data.frame(Hmisc::hoeffd(x = data[, centrality2.colname],
                                                      y = data[, centrality1.colname])[3])[1,2], row.names = "Results")

  if(hoeffd[1,2] < 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Significantly dependent",
                                          row.names = "Results")
  } else if(hoeffd[1,2] >= 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Not significantly dependent",
                                          row.names = "Results")
  } else if(hoeffd[1,2] < 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Significantly dependent",
                                          row.names = "Results")
  } else if (hoeffd[1,2] >= 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Not significantly dependent",
                                          row.names = "Results")
  }


  ##assessment of descriptive non-linear non-parametric correlation/dependence between
  #dependent and independent variables

  nl.cor.dep <- NNS::NNS.dep(x = data[, centrality2.colname], y = data[, centrality1.colname],
                        print.map = FALSE, order = 2)

  nl.cor.dep <- data.frame(Correlation = unlist(nl.cor.dep)[1],
                           Dependence = unlist(nl.cor.dep)[2], row.names = "Results")

  ##assessment of conditional probability of deviation of BC and NC from their means in opposite directions
  #filtering the data to find those nodes meeting the conditions
  ncpositive <- data[data[, centrality2.colname] >
                       mean(data[, centrality2.colname]),]
  ncpositive.bcnegative <- ncpositive[ncpositive[, centrality1.colname] <
                                        mean(data[, centrality1.colname]),]
  ncpositive.bcnegative.prob <- (nrow(ncpositive.bcnegative)/nrow(ncpositive))*100

  ncnegative <- data[data[, centrality2.colname] <
                       mean(data[, centrality2.colname]),]
  ncnegative.bcpositive <- ncnegative[ncnegative[, centrality1.colname] >
                                        mean(data[, centrality1.colname]),]
  ncnegative.bcpositive.prob <- (nrow(ncnegative.bcpositive)/nrow(ncnegative))*100

  #calculation of conditional probability
  final.cond.prob <- sum(ncpositive.bcnegative.prob, ncnegative.bcpositive.prob)/2

  ##Reliability analysis based on split-half random sampling

  #split-half random sampling

  sample.data <- data[sample(1:nrow(data), size = round(nrow(data)/2), replace = FALSE),]

  #filtering the data to find those nodes meeting the conditions
  sample.ncpositive <- sample.data[sample.data[, centrality2.colname] >
                                     mean(sample.data[, centrality2.colname]),]
  sample.ncpositive.bcnegative <- sample.ncpositive[sample.ncpositive[, centrality1.colname] <
                                                      mean(sample.data[, centrality1.colname]),]
  sample.ncpositive.bcnegative.prob <- (nrow(sample.ncpositive.bcnegative)/nrow(sample.ncpositive))*100

  sample.ncnegative <- sample.data[sample.data[, centrality2.colname] <
                                     mean(sample.data[, centrality2.colname]),]
  sample.ncnegative.bcpositive <- sample.ncnegative[sample.ncnegative[, centrality1.colname] >
                                                      mean(sample.data[, centrality1.colname]),]
  sample.ncnegative.bcpositive.prob <- (nrow(sample.ncnegative.bcpositive)/nrow(sample.ncnegative))*100

  #calculation of conditional probability
  sample.final.cond.prob <- sum(sample.ncpositive.bcnegative.prob, sample.ncnegative.bcpositive.prob)/2

  results <- list(Summary_statistics = summary.stat,
                  Normality_results = normality,
                  Centrality1_Normality = dependent.normality,
                  Centrality2_Normality = independent.normality,
                  HoeffdingD_Statistic = hoeffd,
                  Dependence_Significance = dependence.significance,
                  NNS_dep_results = nl.cor.dep,
                  ConditionalProbability = final.cond.prob,
                  ConditionalProbability_split.half.sample = sample.final.cond.prob)

  return(results)

}

#=============================================================================
#
#    Code chunk 9: Calculation of IVI from centrality measures
#
#=============================================================================

#' Integrated Value of Influence (IVI)
#'
#' This function calculates the IVI of the desired nodes from previously calculated centrality
#' measures. This function is not dependent to other packages and the required centrality
#' measures, namely degree centrality, ClusterRank, betweenness centrality, Collective Influence,
#' local H-index, and neighborhood connectivity could have been calculated by any means beforehand.
#' @param DC A vector containing the values of degree centrality of the desired vertices.
#' @param CR A vector containing the values of ClusterRank of the desired vertices.
#' @param LH_index A vector containing the values of local H-index of the desired vertices.
#' @param NC A vector containing the values of neighborhood connectivity of the desired vertices.
#' @param BC A vector containing the values of betweenness centrality of the desired vertices.
#' @param CI A vector containing the values of Collective Influence of the desired vertices.
#' @param scaled Logical; whether the end result should be 1-100 range normalized or not (default is TRUE).
#' @return A numeric vector with the IVI values based on the provided centrality measures.
#' @aliases IVI.FI
#' @keywords ivi.from.indices
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{ivi}},
#' \code{\link[influential]{exir}}
#' @export
#' @examples
#' MyData <- centrality.measures
#' My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,
#'                                     CR = centrality.measures$CR,
#'                                     NC = centrality.measures$NC,
#'                                     LH_index = centrality.measures$LH_index,
#'                                     BC = centrality.measures$BC,
#'                                     CI = centrality.measures$CI)
ivi.from.indices <- function(DC, CR, LH_index, NC, BC, CI, scaled = TRUE) {

  #Generating temporary measures

  temp.DC <- DC
  temp.CR <- CR
  temp.LH_index <- LH_index
  temp.NC <- NC
  temp.BC <- BC
  temp.CI <- CI

  #Removing the NAN and NA values

  temp.DC[c(which(is.nan(temp.DC)), which(is.na(temp.DC)))] <- 0
  temp.CR[c(which(is.nan(temp.CR)), which(is.na(temp.CR)))] <- 0
  temp.LH_index[c(which(is.nan(temp.LH_index)), which(is.na(temp.LH_index)))] <- 0
  temp.NC[c(which(is.nan(temp.NC)), which(is.na(temp.NC)))] <- 0
  temp.BC[c(which(is.nan(temp.BC)), which(is.na(temp.BC)))] <- 0
  temp.CI[c(which(is.nan(temp.CI)), which(is.na(temp.CI)))] <- 0

  #1-100 normalization of centrality measures

  if(any(temp.DC > 0)) {
    temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
  }

  if(any(temp.CR > 0)) {
    temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
  }

  if(any(temp.LH_index > 0)) {
    temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
  }

  if(any(temp.NC > 0)) {
    temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
  }

  if(any(temp.BC > 0)) {
    temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
  }

  if(any(temp.CI > 0)) {
    temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
  }

  #Calculation of IVI

  spreading.rank <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))

  hubness.rank <- (temp.DC+temp.LH_index)

  temp.ivi <- (hubness.rank)*(spreading.rank)

  #1-100 normalization of IVI

  if(scaled == TRUE) {

    temp.ivi <- 1+(((temp.ivi-min(temp.ivi))*(100-1))/(max(temp.ivi)-min(temp.ivi)))

  }

  return(temp.ivi)
}

#=============================================================================
#
#    Code chunk 10: Calculation of IVI from graph
#
#=============================================================================

#' Integrated Value of Influence (IVI)
#'
#' This function calculates the IVI of the desired nodes from a graph.
#' @param graph A graph (network) of the igraph class.
#' @param vertices A vector of desired vertices, which could be obtained by the V function.
#' @param weights Optional positive weight vector for calculating weighted betweenness centrality
#' of nodes as a requirement for calculation of IVI. If the graph has a weight edge attribute,
#' then this is used by default. Weights are used to calculate weighted shortest paths,
#' so they are interpreted as distances.
#' @param directed Logical scalar, whether to directed graph is analyzed. This argument
#' is ignored for undirected graphs.
#' @param mode The mode of IVI depending on the directedness of the graph.
#' If the graph is undirected, the mode "all" should be specified.
#' Otherwise, for the calculation of IVI based on
#' incoming connections select "in" and for the outgoing connections select "out".
#' Also, if all of the connections are desired, specify the "all" mode. Default mode is set to "all".
#' @param loops Logical; whether the loop edges are also counted.
#' @param d The distance, expressed in number of steps from a given node (default=3). Distance
#' must be > 0. According to Morone & Makse (https://doi.org/10.1038/nature14604), optimal
#' results can be reached at d=3,4, but this depends on the size/"radius" of the network.
#' NOTE: the distance d is not inclusive. This means that nodes at a distance of 3 from
#' our node-of-interest do not include nodes at distances 1 and 2. Only 3.
#' @param scaled Logical; whether the end result should be 1-100 range normalized or not (default is TRUE).
#' @return A numeric vector with the IVI values based on the provided centrality measures.
#' @aliases IVI
#' @keywords IVI integrated_value_of_influence
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{ivi.from.indices}},
#' \code{\link[influential]{exir}}
#' @export
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' My.vertices.IVI <- ivi(graph = My_graph, vertices = GraphVertices,
#'                        weights = NULL, directed = FALSE, mode = "all",
#'                        loops = TRUE, d = 3, scaled = TRUE)
#' }
ivi <- function(graph, vertices = V(graph), weights = NULL, directed = FALSE,
                mode = "all", loops = TRUE, d = 3, scaled = TRUE) {

  #Calculation of required centrality measures

  DC <- igraph::degree(graph = graph, v = vertices, mode = mode, loops = loops)
  CR <- clusterRank(graph = graph, vids = vertices, directed = directed, loops = loops)
  LH_index <- lh_index(graph = graph, vertices = vertices, mode = mode)
  NC <- neighborhood.connectivity(graph = graph, vertices = vertices, mode = mode)
  BC <- betweenness(graph = graph, v = vertices, directed = directed, weights = weights)
  CI <- collective.influence(graph = graph, vertices = vertices, mode = mode, d = d)

  #Generating temporary measures

  temp.DC <- DC
  temp.CR <- CR
  temp.LH_index <- LH_index
  temp.NC <- unlist(NC)
  temp.BC <- BC
  temp.CI <- CI

  #Removing the NAN and NA values

  temp.DC[c(which(is.nan(temp.DC)), which(is.na(temp.DC)))] <- 0
  temp.CR[c(which(is.nan(temp.CR)), which(is.na(temp.CR)))] <- 0
  temp.LH_index[c(which(is.nan(temp.LH_index)), which(is.na(temp.LH_index)))] <- 0
  temp.NC[c(which(is.nan(temp.NC)), which(is.na(temp.NC)))] <- 0
  temp.BC[c(which(is.nan(temp.BC)), which(is.na(temp.BC)))] <- 0
  temp.CI[c(which(is.nan(temp.CI)), which(is.na(temp.CI)))] <- 0

  #1-100 normalization of centrality measures

  if(any(temp.DC > 0)) {
    temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
  }

  if(any(temp.CR > 0)) {
    temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
  }

  if(any(temp.LH_index > 0)) {
    temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
  }

  if(any(temp.NC > 0)) {
    temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
  }

  if(any(temp.BC > 0)) {
    temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
  }

  if(any(temp.CI > 0)) {
    temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
  }

  #Calculation of IVI

  spreading.rank <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))

  hubness.rank <- (temp.DC+temp.LH_index)

  temp.ivi <- (hubness.rank)*(spreading.rank)

  #1-100 normalization of IVI

  if(scaled == TRUE) {

    temp.ivi <- 1+(((temp.ivi-min(temp.ivi))*(100-1))/(max(temp.ivi)-min(temp.ivi)))

  }

  return(temp.ivi)
}

#=============================================================================
#
#    Code chunk 11: Calculation of Spreading score
#
#=============================================================================

#' Spreading score
#'
#' This function calculates the Spreading score of the desired nodes from a graph.
#' Spreading score reflects the spreading potential of each node within a network and is
#' one of the major components of the IVI.
#' @param graph A graph (network) of the igraph class.
#' @param vertices A vector of desired vertices, which could be obtained by the V function.
#' @param weights Optional positive weight vector for calculating weighted betweenness centrality
#' of nodes as a requirement for calculation of IVI. If the graph has a weight edge attribute,
#' then this is used by default. Weights are used to calculate weighted shortest paths,
#' so they are interpreted as distances.
#' @param directed Logical scalar, whether to directed graph is analyzed. This argument
#' is ignored for undirected graphs.
#' @param mode The mode of Spreading score depending on the directedness of the graph.
#' If the graph is undirected, the mode "all" should be specified.
#' Otherwise, for the calculation of Spreading score based on
#' incoming connections select "in" and for the outgoing connections select "out".
#' Also, if all of the connections are desired, specify the "all" mode. Default mode is set to "all".
#' @param loops Logical; whether the loop edges are also counted.
#' @param d The distance, expressed in number of steps from a given node (default=3). Distance
#' must be > 0. According to Morone & Makse (https://doi.org/10.1038/nature14604), optimal
#' results can be reached at d=3,4, but this depends on the size/"radius" of the network.
#' NOTE: the distance d is not inclusive. This means that nodes at a distance of 3 from
#' our node-of-interest do not include nodes at distances 1 and 2. Only 3.
#' @param scaled Logical; whether the end result should be 1-100 range normalized or not (default is TRUE).
#' @return A numeric vector with Spreading scores.
#' @keywords spreading.score
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{hubness.score}}
#' @export
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' Spreading.score <- spreading.score(graph = My_graph, vertices = GraphVertices,
#'                                    weights = NULL, directed = FALSE, mode = "all",
#'                                    loops = TRUE, d = 3, scaled = TRUE)
#' }
spreading.score <- function(graph, vertices = V(graph), weights = NULL, directed = FALSE,
                            mode = "all", loops = TRUE, d = 3, scaled = TRUE) {

  #Calculation of required centrality measures

  CR <- clusterRank(graph = graph, vids = vertices, directed = directed, loops = loops)
  CR[which(is.nan(CR))] <- 0
  NC <- neighborhood.connectivity(graph = graph, vertices = vertices, mode = mode)
  BC <- betweenness(graph = graph, v = vertices, directed = directed, weights = weights)
  CI <- collective.influence(graph = graph, vertices = vertices, mode = mode, d = d)

  #Generating temporary measures

  temp.CR <- CR
  temp.NC <- unlist(NC)
  temp.BC <- BC
  temp.CI <- CI

  #Removing the NAN and NA values

  temp.CR[c(which(is.nan(temp.CR)), which(is.na(temp.CR)))] <- 0
  temp.NC[c(which(is.nan(temp.NC)), which(is.na(temp.NC)))] <- 0
  temp.BC[c(which(is.nan(temp.BC)), which(is.na(temp.BC)))] <- 0
  temp.CI[c(which(is.nan(temp.CI)), which(is.na(temp.CI)))] <- 0

  #1-100 normalization of centrality measures

  if(any(temp.CR > 0)) {
    temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
  }

  if(any(temp.NC > 0)) {
    temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
  }

  if(any(temp.BC > 0)) {
    temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
  }

  if(any(temp.CI > 0)) {
    temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
  }

  #Calculation of spreading.score

  temp.spreading.score <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))

  #1-100 normalization of spreading score

  if(scaled == TRUE) {

    temp.spreading.score <- 1+(((temp.spreading.score-min(temp.spreading.score))*(100-1))/(max(temp.spreading.score)-min(temp.spreading.score)))

  }

  return(temp.spreading.score)
}

#=============================================================================
#
#    Code chunk 12: Calculation of Hubness score
#
#=============================================================================

#' Hubness score
#'
#' This function calculates the Hubness score of the desired nodes from a graph.
#' Hubness score reflects the power of each node in its surrounding environment and is
#' one of the major components of the IVI.
#' @param graph A graph (network) of the igraph class.
#' @param vertices A vector of desired vertices, which could be obtained by the V function.
#' @param directed Logical scalar, whether to directed graph is analyzed. This argument
#' is ignored for undirected graphs.
#' @param mode The mode of Hubness score depending on the directedness of the graph.
#' If the graph is undirected, the mode "all" should be specified.
#' Otherwise, for the calculation of Hubness score based on
#' incoming connections select "in" and for the outgoing connections select "out".
#' Also, if all of the connections are desired, specify the "all" mode. Default mode is set to "all".
#' @param loops Logical; whether the loop edges are also counted.
#' @param scaled Logical; whether the end result should be 1-100 range normalized or not (default is TRUE).
#' @return A numeric vector with the Hubness scores.
#' @keywords hubness.score
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{spreading.score}}
#' @export
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' Hubness.score <- hubness.score(graph = My_graph, vertices = GraphVertices,
#'                                directed = FALSE, mode = "all",
#'                                loops = TRUE, scaled = TRUE)
#' }
hubness.score <- function(graph, vertices = V(graph), directed = FALSE,
                          mode = "all", loops = TRUE, scaled = TRUE) {

  #Calculation of required centrality measures

  DC <- igraph::degree(graph = graph, v = vertices, mode = mode, loops = loops)
  LH_index <- lh_index(graph = graph, vertices = vertices, mode = mode)

  #Generating temporary measures

  temp.DC <- DC
  temp.LH_index <- LH_index

  #Removing the NAN and NA values

  temp.DC[c(which(is.nan(temp.DC)), which(is.na(temp.DC)))] <- 0
  temp.LH_index[c(which(is.nan(temp.LH_index)), which(is.na(temp.LH_index)))] <- 0

  #1-100 normalization of centrality measures

  if(any(temp.DC > 0)) {
    temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
  }

  if(any(temp.LH_index > 0)) {
    temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
  }

  #Calculation of Hubness score

  temp.hubness.score <- (temp.DC+temp.LH_index)

  #1-100 normalization of Hubness score

  if(scaled == TRUE) {

    temp.hubness.score <- 1+(((temp.hubness.score-min(temp.hubness.score))*(100-1))/(max(temp.hubness.score)-min(temp.hubness.score)))

  }

  return(temp.hubness.score)
}

#=============================================================================
#
#    Code chunk 13: Calculation of SIRIR
#
#=============================================================================

#' SIR-based Influence Ranking
#'
#' This function is achieved by the integration susceptible-infected-recovered (SIR) model
#' with the leave-one-out cross validation technique and ranks network nodes based on their
#' true universal influence. One of the applications of this function is the assessment of
#' performance of a novel algorithm in identification of network influential nodes.
#' @param graph A graph (network) of the igraph class.
#' @param vertices A vector of desired vertices, which could be obtained by the V function.
#' @param beta Non-negative scalar. The rate of infection of an individual that is susceptible
#' and has a single infected neighbor. The infection rate of a susceptible individual with n
#' infected neighbors is n times beta. Formally this is the rate parameter of an exponential
#' distribution.
#' @param gamma Positive scalar. The rate of recovery of an infected individual.
#' Formally, this is the rate parameter of an exponential distribution.
#' @param no.sim Integer scalar, the number of simulation runs to perform SIR model on for the
#' original network as well perturbed networks generated by leave-one-out technique.
#' You may choose a different no.sim based on the available memory on your system.
#' @param seed A single value, interpreted as an integer to be used for random number generation
#' @return A two-column dataframe; a column containing the difference values of the original and
#' perturbed networks and a column containing node influence rankings
#' @aliases SIRIR
#' @keywords sirir
#' @seealso \code{\link[igraph]{sir}} for a complete description on SIR model.
#' @export
#' @examples
#' set.seed(1234)
#' My_graph <- igraph::sample_gnp(n=50, p=0.05)
#' GraphVertices <- V(My_graph)
#' Influence.Ranks <- sirir(graph = My_graph, vertices = GraphVertices,
#'                          beta = 0.5, gamma = 1, no.sim = 10, seed = 1234)
#' @importFrom igraph vcount as_ids sir
sirir <- function(graph, vertices = V(graph),
                  beta = 0.5, gamma = 1,
                  no.sim = igraph::vcount(graph)*100,  seed = 1234) {

  #Define a data frame
  temp.loocr.table <- data.frame(difference.value = vector("numeric", length = length(vertices)),
                                 rank = vector("integer", length = length(vertices)))

  if(class(vertices) == "character") {
    rownames(temp.loocr.table) <- vertices
  } else if(class(vertices) == "igraph.vs") {
    rownames(temp.loocr.table) <- igraph::as_ids(vertices)
  }


  #Model the spreading based on all nodes
  set.seed(seed)
  all.included.spread <- igraph::sir(graph = graph, beta = beta,
                             gamma = gamma, no.sim = no.sim)

  #Getting the mean of spread in each independent experiment
  all.mean.spread <- vector("numeric", length = length(all.included.spread))

  for (i in 1:length(all.included.spread)) {
    all.mean.spread[i] <- max(all.included.spread[[i]]$NR)
  }
  all.mean.spread <- mean(all.mean.spread)

  #Model the spread based on leave one out cross ranking (LOOCR)

  for(s in 1:length(vertices)) {

    temp.graph <- igraph::delete_vertices(graph, unlist(vertices[s]))

    set.seed(seed)

    loocr.spread <- igraph::sir(graph = temp.graph, beta = beta,
                        gamma = gamma, no.sim = no.sim)

    loocr.mean.spread <- vector("numeric", length = length(loocr.spread))

    for (h in 1:length(loocr.spread)) {
      loocr.mean.spread[h] <- max(loocr.spread[[h]]$NR)
    }
    loocr.mean.spread <- mean(loocr.mean.spread)
    temp.loocr.table$difference.value[s] <- all.mean.spread - loocr.mean.spread
  }

  temp.loocr.table$rank <- rank(-temp.loocr.table$difference.value, ties.method = "min")

  return(temp.loocr.table)
}

#=============================================================================
#
#    Code chunk 14: Some required functions from the igraph package
#
#=============================================================================

#' Creating igraph graphs from data frames
#'
#' This function and all of its descriptions have been obtained from the igraph package.
#' For a complete description if the function and its arguments try this:
#' ?igraph::graph_from_data_frame
#' @param d A data frame containing a symbolic edge list in the first two columns.
#' Additional columns are considered as edge attributes.
#' Since version 0.7 this argument is coerced to a data frame with as.data.frame.
#' @param directed Logical scalar, whether or not to create a directed graph.
#' @param vertices A data frame with vertex metadata, or NULL.
#' Since version 0.7 of igraph this argument is coerced to a data frame with as.data.frame, if not NULL.
#' @return An igraph graph object.
#' @aliases dataframe2graph
#' @keywords graph_from_dataframe
#' @family network_reconstruction functions
#' @seealso \code{\link[igraph]{graph_from_adjacency_matrix}} for a complete description on this function
#' @export
#' @examples
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(d=MyData)
#' @importFrom igraph graph_from_data_frame
  graph_from_data_frame <- igraph::graph_from_data_frame


  #*****************************************************************#

  #' Creating igraph graphs from adjacency matrices
  #'
  #' This function and all of its descriptions have been obtained from the igraph package.
  #' For a complete description if the function and its arguments try this:
  #' ?igraph::graph_from_adjacency_matrix
  #' @param adjmatrix A square adjacency matrix. From igraph version 0.5.1 this
  #' can be a sparse matrix created with the Matrix package.
  #' @param mode Character scalar, specifies how igraph should interpret the supplied matrix.
  #' See also the weighted argument, the interpretation depends on that too.
  #' Possible values are: directed, undirected, upper, lower, max, min, plus.
  #' @param weighted This argument specifies whether to create a weighted graph from an adjacency matrix.
  #' If it is NULL then an unweighted graph is created and the elements of the adjacency matrix gives the
  #' number of edges between the vertices. If it is a character constant then for every non-zero matrix
  #' entry an edge is created and the value of the entry is added as an edge attribute named by the weighted argument.
  #' If it is TRUE then a weighted graph is created and the name of the edge attribute will be weight.
  #' @param diag Logical scalar, whether to include the diagonal of the matrix in the calculation.
  #' If this is FALSE then the diagonal is zerod out first.
  #' @param add.colnames Character scalar, whether to add the column names as vertex attributes.
  #' If it is ‘NULL’ (the default) then, if present, column names are added as vertex attribute ‘name’.
  #' If ‘NA’ then they will not be added. If a character constant, then it gives the name of the vertex attribute to add.
  #' @param add.rownames Character scalar, whether to add the row names as vertex attributes.
  #' Possible values the same as the previous argument. By default row names are not added.
  #' If ‘add.rownames’ and ‘add.colnames’ specify the same vertex attribute, then the former is ignored.
  #' @return An igraph graph object.
  #' @aliases adjmatrix2graph
  #' @keywords graph_from_adjacencymatrices
  #' @family network_reconstruction functions
  #' @seealso \code{\link[igraph]{graph_from_adjacency_matrix}} for a complete description on this function
  #' @export
  #' @examples
  #' MyData <- coexpression.adjacency
  #' My_graph <- graph_from_adjacency_matrix(MyData)
  #' @importFrom igraph graph_from_adjacency_matrix
  graph_from_adjacency_matrix <- igraph::graph_from_adjacency_matrix

  #*****************************************************************#

  #' Vertices of an igraph graph
  #'
  #' This function and all of its descriptions have been obtained from the igraph package.
  #' @param graph The graph (an igraph graph)
  #' @return A vertex sequence containing all vertices, in the order of their numeric vertex ids.
  #' @aliases vertices
  #' @keywords graph_vertices
  #' @seealso \code{\link[igraph]{V}} for a complete description on this function
  #' @export
  #' @examples
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' My_graph_vertices <- V(My_graph)
  #' @importFrom igraph V
  V <- igraph::V

  #*****************************************************************#

  #' Vertex betweenness centrality
  #'
  #' This function and all of its descriptions have been obtained from the igraph package.
  #' @param graph The graph to analyze (an igraph graph).
  #' @param v The vertices for which the vertex betweenness will be calculated.
  #' @param directed Logical, whether directed paths should be considered while determining the shortest paths.
  #' @param weights Optional positive weight vector for calculating weighted betweenness.
  #' If the graph has a weight edge attribute, then this is used by default. Weights are used to calculate weighted shortest paths, so they are interpreted as distances.
  #' @param nobigint Logical scalar, whether to use big integers during the calculation.
  #' This is only required for lattice-like graphs that have very many shortest paths between a pair of vertices.
  #' If TRUE (the default), then big integers are not used.
  #' @param normalized Logical scalar, whether to normalize the betweenness scores. If TRUE, then the results are normalized.
  #' @return A numeric vector with the betweenness score for each vertex in v.
  #' @aliases BC
  #' @keywords betweenness_centrality
  #' @family centrality functions
  #' @seealso \code{\link[igraph]{betweenness}} for a complete description on this function
  #' @export
  #' @examples
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,
  #'                         directed = FALSE, normalized = FALSE)
  #' @importFrom igraph betweenness
  betweenness <- igraph::betweenness

  #*****************************************************************#

  #' Degree of the vertices
  #'
  #' This function and all of its descriptions have been obtained from the igraph package.
  #' @param graph The graph to analyze (an igraph graph).
  #' @param v The ids of vertices of which the degree will be calculated.
  #' @param mode Character string, “out” for out-degree, “in” for in-degree or “total” for the sum of the two.
  #' For undirected graphs this argument is ignored. “all” is a synonym of “total”.
  #' @param loops Logical; whether the loop edges are also counted.
  #' If the graph has a weight edge attribute, then this is used by default. Weights are used to calculate weighted shortest paths, so they are interpreted as distances.
  #' @param normalized Logical scalar, whether to normalize the degree.
  #' If TRUE then the result is divided by n-1, where n is the number of vertices in the graph.
  #' @return A numeric vector of the same length as argument v.
  #' @aliases DC
  #' @keywords degree_centrality
  #' @family centrality functions
  #' @seealso \code{\link[igraph]{degree}} for a complete description on this function
  #' @export
  #' @examples
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' My_graph_degree <- degree(My_graph, v = GraphVertices, normalized = FALSE)
  #' @importFrom igraph degree
  degree <- igraph::degree

#=============================================================================
#
#    Code chunk 15: Funcstions for supporting additional file formats
#
#=============================================================================

  #' SIF to igraph
  #'
  #' This function imports and converts a SIF file from your local hard drive, cloud space,
  #' or internet into a graph with an igraph class, which can then be used for the identification
  #' of most influential nodes via the ivi function.
  #' @param Path A string or character vector indicating the path to the desired SIF file. The SIF file
  #' could be on your local hard drive, cloud space, or on the internet.
  #' @param directed Logical scalar, whether or not to create a directed graph.
  #' @return An igraph graph object.
  #' @keywords SIF.to.igraph
  #' @seealso \code{\link[influential]{graph_from_data_frame}}
  #' @export
  #' @examples
  #' \dontrun{
  #' MyGraph <- sif2igraph(Path = "/Users/User1/Desktop/mygraph.sif", directed=FALSE)
  #' }
  sif2igraph <- function(Path, directed=FALSE) {

    graph_from_data_frame(d = utils::read.delim(Path, header = FALSE)[,c(1,3)],directed = directed)
  }

#=============================================================================
#
#    Code chunk 16: Calculation of ExIR
#
#=============================================================================

  #' Experimental data-based Integrated Ranking
  #'
  #' This function runs the Experimental data-based Integrated Ranking (ExIR)
  #' model for the classification and ranking of top candidate features. The input
  #' data could come from any type of experiment such as transcriptomics and proteomics.
  #' @param Desired_list (Optional) A character vector of your desired features. This vector could be, for
  #' instance, a list of features obtained from cluster analysis, time-course analysis,
  #' or a list of dysregulated features with a specific sign.
  #' @param Diff_data A dataframe of all significant differential/regression data and their
  #' statistical significance values (p-value/adjusted p-value).
  #' You may have selected a proportion of the differential data as the significant ones according
  #' to your desired thresholds. A function, named diff.data.assembly, has also been
  #' provided for the convenient assembling of the Diff_data dataframe.
  #' @param Diff_value A numeric vector containing the column number(s) of the differential
  #' data in the Diff_data dataframe. The differential data could result from any type of
  #' differential data analysis. One example could be the fold changes (FCs) obtained from differential
  #' expression analyses. The user may provide as many differential data as he/she wish.
  #' @param Regr_value (Optional) A numeric vector containing the column number(s) of the regression
  #' data in the Diff_data dataframe. The regression data could result from any type of regression
  #' data analysis or other analyses such as time-course data analyses that are based on regression models.
  #' @param Sig_value A numeric vector containing the column number(s) of the significance values (p-value/adjusted p-value) of
  #' both differential and regression data (if provided). Providing significance values for the regression data is optional.
  #' @param Exptl_data A dataframe containing all of the experimental including a column for specifying the conditions.
  #' The features/variables of the dataframe should be as the columns and the samples should come in the rows.
  #' The condition column should be of the character class. For example, if the study includes several replicates of
  #' cancer and normal samples, the condition column should include "cancer" and "normal" as the conditions of different samples.
  #' Also, the prior normalization of the experimental data is highly recommended. Otherwise,
  #' the user may set the Normalize argument to TRUE for a simple log2 transformation of the data.
  #' The experimental data could come from a variety sources such as transcriptomics and proteomics assays.
  #' @param Condition_colname A string or character vector specifying the name of the condition column of the Exptl_data dataframe.
  #' @param Normalize Logical; whether the experimental data should be normalized or not (default is FALSE). If TRUE, the
  #' experimental data will be log2 transformed.
  #' @param r The threshold of Pearson correlation coefficient for the selection of correlated features (default is 0).
  #' @param alpha The threshold of the statistical significance (p-value) used throughout the entir model (default is 0.05)
  #' @param num_trees Number of trees to be used for the random forest classification (supervised machine learning) Default is set to 10000.
  #' @param num_permutations Number of permutations to be used for computation of the statistical significances (p-values) of
  #' the importance scores resulted from random forest classification (default is 100).
  #' @param seed The seed to be used for all of the random processes throughout the model (default is 1234).
  #' @param verbose Logical; whether the accomplishment of different stages of the model should be printed (default is TRUE).
  #' @return A list of one to four tables including:
  #'
  #' - Driver table: Top candidate drivers
  #'
  #' - DE-mediator table: Top candidate differentially expressed/abundant mediators
  #'
  #' - nonDE-mediators table: Top candidate non-differentially expressed/abundant mediators
  #'
  #' - Biomarker table: Top candidate biomarkers
  #'
  #' The number of returned tables depends on the input data and specified arguments.
  #' @aliases ExIR
  #' @keywords exir
  #' @family integrative ranking functions
  #' @seealso \code{\link[influential]{diff.data.assembly}},
  #' \code{\link[influential]{ivi}},
  #' \code{\link[coop]{pcor}},
  #' \code{\link[stats]{prcomp}},
  #' \code{\link[ranger]{ranger}},
  #' \code{\link[ranger]{importance_pvalues}}
  #' @export
  #' @example
  #' \dontrun{
  #' MyDesired_list <- Desiredlist
  #' MyDiff_data <- Diffdata
  #' Diff_value <- c(1,3,5)
  #' Regr_value <- 7
  #' Sig_value <- c(2,4,6,8)
  #' MyExptl_data <- Exptldata
  #' Condition_colname <- "condition"
  #' My.exir <- exir(Desired_list = MyDesired_list,
  #'                Diff_data = MyDiff_data, Diff_value = Diff_value,
  #'                Regr_value = Regr_value, Sig_value = Sig_value,
  #'                Exptl_data = MyExptl_data, Condition_colname = Condition_colname)
  #' }
  exir <- function(Desired_list = NULL,
                   Diff_data, Diff_value, Regr_value = NULL, Sig_value,
                   Exptl_data, Condition_colname, Normalize = FALSE,
                   r = 0, alpha = 0.05, num_trees = 10000, num_permutations = 100,
                   seed = 1234, verbose = TRUE) {

    # Setup progress bar
    if(verbose) {
      pb <- utils::txtProgressBar(min = 0, max = 100, initial = 0, char = "=",
                           width = NA, style = 3, file = "")
    }

    #ProgressBar: Preparing the input data
    if(verbose) {
      print("Preparing the input data", quote = FALSE)
    }

    #make sure the input data is of data frame class
    Diff_data <- as.data.frame(Diff_data)
    Exptl_data <- as.data.frame(Exptl_data)

    # Get the column number of condition column
    condition.index <- match(Condition_colname, colnames(Exptl_data))

    # Normalize the experimental data (if required)
    if(Normalize) {
      Exptl_data[,-condition.index] <- log2(Exptl_data[,-condition.index]+1)
    }

    #ProgressBar: Preparing the input data
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 5)
    }

    #ProgressBar: Calculating the differential score
    if(verbose) {
      print("Calculating the differential score", quote = FALSE)
    }

    #1 Calculate differential score
    Diff_data$sum.Diff_value <- base::abs(base::apply(Diff_data[,Diff_value, drop = FALSE],1,sum))
    #range normalize the differential score
    Diff_data$sum.Diff_value <- 1+(((Diff_data$sum.Diff_value-min(Diff_data$sum.Diff_value))*(100-1))/
                                     (max(Diff_data$sum.Diff_value)-min(Diff_data$sum.Diff_value)))

    #ProgressBar: Calculating the differential score
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 10)
    }

    #ProgressBar: Calculating the regression/time-course R-squared score
    if(verbose & !is.null(Regr_value)) {
      print("Calculating the regression/time-course R-squared score", quote = FALSE)
    }

    #2 Calculate regression/time-course R-squared score (if provided)
    if (!is.null(Regr_value)) {
      Diff_data$sum.Regr_value <- base::apply(Diff_data[,Regr_value, drop = FALSE],1,sum)
      #range normalize the R-squared score
      Diff_data$sum.Regr_value <- 1+(((Diff_data$sum.Regr_value-min(Diff_data$sum.Regr_value))*(100-1))/
                                       (max(Diff_data$sum.Regr_value)-min(Diff_data$sum.Regr_value)))
    }

    #ProgressBar: Calculating the regression/time-course R-squared score
    if(verbose & !is.null(Regr_value)) {
      utils::setTxtProgressBar(pb = pb, value = 15)
    }

    #ProgressBar: Calculating the collective statistical significance of differential/regression factors
    if(verbose) {
      print("Calculating the collective statistical significance of differential/regression factors", quote = FALSE)
    }

    #3 Calculate statistical significance of differential/regression factors
    if (max(Diff_data[,Sig_value]) > 1 | min(Diff_data[,Sig_value]) < 0) {
      stop("input Sig-values (p-value/padj) must all be in the range 0 to 1!")
    }

    if(min(Diff_data[,Sig_value])==0) {
      Diff_data[,Sig_value, drop = FALSE] <- Diff_data[,Sig_value, drop = FALSE] + sort(as.matrix(Diff_data[,Sig_value, drop = FALSE]))[2]

      for (i in 1:ncol(Diff_data[,Sig_value, drop = FALSE])) {
        Diff_data[,Sig_value, drop = FALSE][which(Diff_data[,Sig_value, drop = FALSE][,i] > 1),i] <- 1
      }
    }

    Diff_data$sum.Sig_value <- base::apply(-log10(Diff_data[,Sig_value, drop = FALSE]),1,sum)
    #range normalize the statistical significance
    Diff_data$sum.Sig_value <- 1+(((Diff_data$sum.Sig_value-min(Diff_data$sum.Sig_value))*(100-1))/
                                    (max(Diff_data$sum.Sig_value)-min(Diff_data$sum.Sig_value)))

    #ProgressBar: Calculating the collective statistical significance of differential/regression factors
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 20)
    }

    #ProgressBar: Performing random forest classification (supervised machine learning)
    if(verbose) {
      print("Performing random forest classification (supervised machine learning)", quote = FALSE)
    }

    #4 Calculation of the Integrated Value of Influence (IVI)

    #a Separate a transcriptomic profile of diff features
    if(!is.null(Desired_list)) {
      sig.diff.index <- stats::na.omit(base::unique(base::match(Desired_list,
                                                                colnames(Exptl_data))))
    } else {
      sig.diff.index <- stats::na.omit(base::unique(base::match(rownames(Diff_data),
                                                                colnames(Exptl_data))))
    }

    exptl.for.super.learn <- Exptl_data[,sig.diff.index]
    exptl.for.super.learn$condition <- Exptl_data[,condition.index]

    #b Perform random forest classification
    rf.diff.exptl <- ranger::ranger(seed = seed,
                                    formula = condition ~ .,
                                    data = exptl.for.super.learn,
                                    num.trees = num_trees, importance = "permutation",
                                    write.forest = FALSE)

    rf.diff.exptl.pvalue <- as.data.frame(ranger::importance_pvalues(seed = seed,
                                                                     x = rf.diff.exptl,
                                                                     formula = condition ~ .,
                                                                     num.permutations = num_permutations,
                                                                     data = exptl.for.super.learn,
                                                                     method = "altmann"))

    rf.diff.exptl.pvalue <- base::subset(rf.diff.exptl.pvalue, rf.diff.exptl.pvalue$pvalue <alpha)

    if(min(rf.diff.exptl.pvalue[,"pvalue"])==0) {
      rf.diff.exptl.pvalue[,"pvalue"] <- rf.diff.exptl.pvalue[,"pvalue"] + sort(as.matrix(rf.diff.exptl.pvalue[,"pvalue"]))[2]
    }

    #ProgressBar: Performing random forest classification (supervised machine learning)
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 35)
    }

    #ProgressBar: Performing the first round association analysis
    if(verbose) {
      print("Performing the first round association analysis", quote = FALSE)
    }

    #c Performing correlation analysis
    temp.corr <- coop::pcor(Exptl_data[,-condition.index])

    #filter corr data for only those corr between diff features and themselves/others
    filter.corr.index <- stats::na.omit(base::unique(base::match(rownames(rf.diff.exptl.pvalue),
                                                                 colnames(temp.corr))))
    temp.corr <- temp.corr[,filter.corr.index]

    temp.corr <- reshape2::melt(data = temp.corr, value.name = "cor")

    #filterring low level correlations
    cor.thresh <- r
    temp.corr <- base::subset(temp.corr, temp.corr[,3]>cor.thresh)

    temp.corr <- temp.corr[-which(as.character(temp.corr$Var1)==
                                    as.character(temp.corr$Var2)),]

    if(nrow(temp.corr)>20000) {

      repeat {

        cor.thresh <- cor.thresh+((1-cor.thresh)/2)
        temp.corr <- base::subset(temp.corr, temp.corr[,3]>cor.thresh)

        if(nrow(temp.corr)<=20000) {
          break
        }
      }
    }

    #getting the list of diff features and their correlated features
    diff.plus.corr.features <- base::unique(c(base::as.character(temp.corr[,1]),
                                              base::as.character(temp.corr[,2])))

    #ProgressBar: Performing first round association analysis
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 45)
    }

    #ProgressBar: Performing the second round association analysis
    if(verbose) {
      print("Performing the second round association analysis", quote = FALSE)
    }

    #redo correlation analysis
    temp.corr <- coop::pcor(Exptl_data[,-condition.index])

    #filter corr data for only those corr between diff.plus.corr.features and themselves/others
    filter.corr.index <- stats::na.omit(match(diff.plus.corr.features,
                                              colnames(temp.corr)))
    temp.corr <- temp.corr[,filter.corr.index]

    temp.corr <- reshape2::melt(data = temp.corr, value.name = "cor")

    #filterring low level correlations
    cor.thresh <- r
    temp.corr <- base::subset(temp.corr, temp.corr[,3]>cor.thresh)

    temp.corr <- temp.corr[-which(as.character(temp.corr$Var1)==
                                    as.character(temp.corr$Var2)),]

    if(nrow(temp.corr)>20000) {

      repeat {

        cor.thresh <- cor.thresh+((1-cor.thresh)/2)
        temp.corr <- base::subset(temp.corr, temp.corr[,3]>cor.thresh)

        if(nrow(temp.corr)<=20000) {
          break
        }
      }
    }

    #ProgressBar: Performing second round association analysis
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 55)
    }

    #ProgressBar: Network reconstruction
    if(verbose) {
      print("Network reconstruction", quote = FALSE)
    }

    #d Graph reconstruction
    temp.corr.graph <- igraph::graph_from_data_frame(temp.corr[,c(1:2)])

    #ProgressBar: Network reconstruction
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 60)
    }

    #ProgressBar: Calculation of the integrated value of influence (IVI)
    if(verbose) {
      print("Calculation of the integrated value of influence (IVI)", quote = FALSE)
    }

    #e Calculation of IVI
    temp.corr.ivi <- ivi(temp.corr.graph)

    #ProgressBar: Calculation of the integrated value of influence (IVI)
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 65)
    }

    #ProgressBar: Performing PCA (unsupervised machine learning)
    if(verbose) {
      print("Performing PCA (unsupervised machine learning)", quote = FALSE)
    }

    #5 Unsupervised machine learning (PCA)

    Exptl_data.for.PCA.index <- stats::na.omit(base::match(base::rownames(rf.diff.exptl.pvalue),
                                                           base::colnames(Exptl_data)))
    temp.Exptl_data.for.PCA <- Exptl_data[,Exptl_data.for.PCA.index]

    temp.PCA <- stats::prcomp(temp.Exptl_data.for.PCA)
    temp.PCA.r <- base::abs(temp.PCA$rotation[,1])

    #range normalize the rotation values
    temp.PCA.r <- 1+(((temp.PCA.r-min(temp.PCA.r))*(100-1))/
                       (max(temp.PCA.r)-min(temp.PCA.r)))

    #ProgressBar: Performing PCA (unsupervised machine learning)
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 70)
    }

    #ProgressBar: Calculation of the primitive driver score
    if(verbose) {
      print("Calculation of the primitive driver score", quote = FALSE)
    }

    ## Driver score and ranking

    #a calculate first level driver score based on #3 and #4

    Diff_data$IVI <- 0
    Diff_data.IVI.index <- stats::na.omit(match(names(temp.corr.ivi),
                                                rownames(Diff_data)))

    temp.corr.ivi.for.Diff_data.IVI.index <- stats::na.omit(match(rownames(Diff_data)[Diff_data.IVI.index],
                                                                  names(temp.corr.ivi)))

    Diff_data$IVI[Diff_data.IVI.index] <- temp.corr.ivi[temp.corr.ivi.for.Diff_data.IVI.index]

    #range normalize the IVI
    Diff_data$IVI <- 1+(((Diff_data$IVI-min(Diff_data$IVI))*(100-1))/
                          (max(Diff_data$IVI)-min(Diff_data$IVI)))

    Diff_data$first.Driver.Rank <- 1
    for (i in 1:nrow(Diff_data)) {
      if(c(any(Diff_data[i,Diff_value, drop = FALSE]<0) & any(Diff_data[i,Diff_value, drop = FALSE]>0))) {
        Diff_data$first.Driver.Rank[i] <- 0
      } else {
        Diff_data$first.Driver.Rank[i] <- Diff_data$sum.Sig_value[i]*Diff_data$IVI[i]
      }
    }
    #range normalize (0,100) the first driver rank
    Diff_data$first.Driver.Rank <- 0+(((Diff_data$first.Driver.Rank-min(Diff_data$first.Driver.Rank))*(100-0))/
                                        (max(Diff_data$first.Driver.Rank)-min(Diff_data$first.Driver.Rank)))

    #ProgressBar: Calculation of the primitive driver score
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 75)
    }

    #ProgressBar: Calculation of the neighborhood driver score
    if(verbose) {
      print("Calculation of the neighborhood driver score", quote = FALSE)
    }

    #b (#6) calculate neighborhood score

    #get the list of network nodes
    network.nodes <- igraph::as_ids(igraph::V(temp.corr.graph))

    neighborehood.score.table <- data.frame(node = network.nodes,
                                            N.score = 0)
    for (n in 1:nrow(neighborehood.score.table)) {
      first.neighbors <- igraph::as_ids(igraph::neighbors(graph = temp.corr.graph,
                                                          v = neighborehood.score.table$node[n],
                                                          mode = "all"))
      first.neighbors.index <- stats::na.omit(match(first.neighbors,
                                                    rownames(Diff_data)))

      neighborehood.score.table$N.score[n] <- sum(Diff_data$first.Driver.Rank[first.neighbors.index])
    }

    #ProgressBar: Calculation of the neighborhood driver score
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 80)
    }

    #ProgressBar: Preparation of the driver table
    if(verbose) {
      print("Preparation of the driver table", quote = FALSE)
    }

    Diff_data$N.score <- 0
    Diff_data.N.score.index <- stats::na.omit(match(neighborehood.score.table$node,
                                                    rownames(Diff_data)))

    neighborehood.score.table.for.Diff_data.N.score.index <- stats::na.omit(match(rownames(Diff_data)[Diff_data.N.score.index],
                                                                                  neighborehood.score.table$node))

    Diff_data$N.score[Diff_data.N.score.index] <- neighborehood.score.table$N.score[neighborehood.score.table.for.Diff_data.N.score.index]

    #range normalize (1,100) the neighborhood score
    Diff_data$N.score <- 1+(((Diff_data$N.score-min(Diff_data$N.score))*(100-1))/
                              (max(Diff_data$N.score)-min(Diff_data$N.score)))

    #c calculate the final driver score

    Diff_data$final.Driver.score <- (Diff_data$first.Driver.Rank)*(Diff_data$N.score)
    Diff_data$final.Driver.score[Diff_data$final.Driver.score==0] <- NA

    # Create the Drivers table

    Driver.table <- Diff_data

    #remove the rows/features with NA in the final driver score
    Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

    #filter the driver table by either the desired list or the list of network nodes
    if(!is.null(Desired_list)) {
      Driver.table.row.index <- stats::na.omit(match(Desired_list,
                                                     rownames(Driver.table)))
    } else {
      Driver.table.row.index <- stats::na.omit(match(network.nodes,
                                                     rownames(Driver.table)))
    }

    Driver.table <- Driver.table[Driver.table.row.index,]
    if(nrow(Driver.table)==0) {Driver.table <- NULL} else {

      #range normalize final driver score
      Driver.table$final.Driver.score <- 1+(((Driver.table$final.Driver.score-min(Driver.table$final.Driver.score))*(100-1))/
                                              (max(Driver.table$final.Driver.score)-min(Driver.table$final.Driver.score)))

      #add driver rank
      Driver.table$rank <- rank(-Driver.table$final.Driver.score, ties.method = "min")

      #add driver type
      Driver.table$driver.type <- ""

      for (d in 1:nrow(Driver.table)) {

        if(sum(Driver.table[d,Diff_value])<0) {
          Driver.table$driver.type[d] <- "Decelerator"

        } else if(sum(Driver.table[d,Diff_value])>0) {
          Driver.table$driver.type[d] <- "Accelerator"
        } else {
          Driver.table$driver.type[d] <- NA
        }
      }

      Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

      #remove redundent columns
      Driver.table <- Driver.table[,c("final.Driver.score",
                                      "rank",
                                      "driver.type")]

      #rename column names
      colnames(Driver.table) <- c("Score", "Rank", "Type")

      #filtering redundant (NaN) results
      Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

      if(nrow(Driver.table)==0) {Driver.table <- NULL}

    }

    #ProgressBar: Preparation of the driver table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 85)
    }

    #ProgressBar: Preparation of the DE-mediator table
    if(verbose) {
      print("Preparation of the DE-mediator table", quote = FALSE)
    }

    # Create the DE mediators table

    DE.mediator.table <- Diff_data

    #include only rows/features with NA in the final driver score (which are mediators)
    DE.mediator.index <- which(is.na(DE.mediator.table$final.Driver.score))

    DE.mediator.table <- DE.mediator.table[DE.mediator.index,]

    DE.mediator.table$DE.mediator.score <- DE.mediator.table$sum.Sig_value*
      DE.mediator.table$IVI*
      DE.mediator.table$N.score

    #filter the DE mediators table by either the desired list or the list of network nodes
    if(!is.null(Desired_list)) {
      DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                    rownames(DE.mediator.table)))
    } else {
      DE.mediator.row.index <- stats::na.omit(match(network.nodes,
                                                    rownames(DE.mediator.table)))
    }

    DE.mediator.table <- DE.mediator.table[DE.mediator.row.index,]
    if(nrow(DE.mediator.table)==0) {DE.mediator.table <- NULL} else {

      #range normalize DE mediators score
      DE.mediator.table$DE.mediator.score <- 1+(((DE.mediator.table$DE.mediator.score-min(DE.mediator.table$DE.mediator.score))*(100-1))/
                                                  (max(DE.mediator.table$DE.mediator.score)-min(DE.mediator.table$DE.mediator.score)))

      #add DE mediators rank
      DE.mediator.table$rank <- rank(-DE.mediator.table$DE.mediator.score, ties.method = "min")

      #remove redundent columns
      DE.mediator.table <- DE.mediator.table[,c("DE.mediator.score", "rank")]

      #rename column names
      colnames(DE.mediator.table) <- c("Score", "Rank")

      #filtering redundant (NaN) results
      DE.mediator.table <- DE.mediator.table[stats::complete.cases(DE.mediator.table),]

      if(nrow(DE.mediator.table)==0) {DE.mediator.table <- NULL}

    }

    #ProgressBar: Preparation of the DE-mediator table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 90)
    }

    #ProgressBar: Preparation of the nonDE-mediator table
    if(verbose) {
      print("Preparation of the nonDE-mediator table", quote = FALSE)
    }

    # Create the non-DE mediators table
    non.DE.mediators.index <- stats::na.omit(unique(match(rownames(Diff_data),
                                                          neighborehood.score.table$node)))

    non.DE.mediators.table <- neighborehood.score.table[-c(non.DE.mediators.index),]
    if(nrow(non.DE.mediators.table)==0) {non.DE.mediators.table <- NULL} else {

      #filter the non-DE mediators table by either the desired list or the list of network nodes
      if(!is.null(Desired_list)) {
        non.DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                          non.DE.mediators.table$node))
        non.DE.mediators.table <- non.DE.mediators.table[non.DE.mediator.row.index,]
      }
    }

    if(nrow(non.DE.mediators.table)==0) {non.DE.mediators.table <- NULL} else {

      rownames(non.DE.mediators.table) <- non.DE.mediators.table$node

      non.DE.mediators.ivi.index <- stats::na.omit(match(rownames(non.DE.mediators.table),
                                                         names(temp.corr.ivi)))

      non.DE.mediators.table$ivi <- temp.corr.ivi[non.DE.mediators.ivi.index]

      non.DE.mediators.table$non.DE.mediator.score <- non.DE.mediators.table$N.score*non.DE.mediators.table$ivi

      #range normalize DE mediators score
      non.DE.mediators.table$non.DE.mediator.score <- 1+(((non.DE.mediators.table$non.DE.mediator.score-min(non.DE.mediators.table$non.DE.mediator.score))*(100-1))/
                                                           (max(non.DE.mediators.table$non.DE.mediator.score)-min(non.DE.mediators.table$non.DE.mediator.score)))

      #add non-DE mediators rank
      non.DE.mediators.table$rank <- rank(-non.DE.mediators.table$non.DE.mediator.score, ties.method = "min")

      #remove redundent columns
      non.DE.mediators.table <- non.DE.mediators.table[,c("non.DE.mediator.score", "rank")]

      #rename column names
      colnames(non.DE.mediators.table) <- c("Score", "Rank")

      #filtering redundant (NaN) results
      non.DE.mediators.table <- non.DE.mediators.table[stats::complete.cases(non.DE.mediators.table),]

      if(nrow(non.DE.mediators.table)==0) {non.DE.mediators.table <- NULL}

    }

    #ProgressBar: Preparation of the nonDE-mediator table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 95)
    }

    #ProgressBar: Preparation of the biomarker table
    if(verbose) {
      print("Preparation of the biomarker table", quote = FALSE)
    }

    # Create the Biomarkers table

    Biomarker.table <- Diff_data

    #remove the rows/features with NA in the final driver score
    Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]

    #filter by random forest table
    Biomarker.rf.index <- stats::na.omit(match(rownames(rf.diff.exptl.pvalue),
                                               rownames(Biomarker.table)))

    Biomarker.table <- Biomarker.table[Biomarker.rf.index,]

    if(nrow(Biomarker.table)==0) {Biomarker.table <- NULL} else {

      #add RF importance score and p-value
      rf.for.Biomarker.table <- stats::na.omit(match(rownames(Biomarker.table),
                                                     rownames(rf.diff.exptl.pvalue)))

      Biomarker.table$rf.importance <- rf.diff.exptl.pvalue$importance[rf.for.Biomarker.table]
      Biomarker.table$rf.pvalue <- rf.diff.exptl.pvalue$pvalue[rf.for.Biomarker.table]

      #range normalize rf.importance and rf.pvalue
      Biomarker.table$rf.importance <- 1+(((Biomarker.table$rf.importance-min(Biomarker.table$rf.importance))*(100-1))/
                                            (max(Biomarker.table$rf.importance)-min(Biomarker.table$rf.importance)))

      Biomarker.table$rf.pvalue <- -log10(Biomarker.table$rf.pvalue)
      Biomarker.table$rf.pvalue <- 1+(((Biomarker.table$rf.pvalue-min(Biomarker.table$rf.pvalue))*(100-1))/
                                        (max(Biomarker.table$rf.pvalue)-min(Biomarker.table$rf.pvalue)))

      #add rotation values
      Biomarker.table.rotation.index <- stats::na.omit(match(rownames(Biomarker.table),
                                                             names(temp.PCA.r)))

      Biomarker.table$rotation <- temp.PCA.r[Biomarker.table.rotation.index]

      #calculate biomarker score
      if(!is.null(Regr_value)) {
        Biomarker.table$final.biomarker.score <- (Biomarker.table$sum.Diff_value)*
          (Biomarker.table$sum.Regr_value)*(Biomarker.table$sum.Sig_value)*
          (Biomarker.table$rf.pvalue)*(Biomarker.table$rf.importance)*
          (Biomarker.table$rotation)
      } else {
        Biomarker.table$final.biomarker.score <- (Biomarker.table$sum.Diff_value)*
          (Biomarker.table$sum.Sig_value)*
          (Biomarker.table$rf.pvalue)*(Biomarker.table$rf.importance)*
          (Biomarker.table$rotation)
      }

      #range normalize biomarker score
      Biomarker.table$final.biomarker.score <- 1+(((Biomarker.table$final.biomarker.score-min(Biomarker.table$final.biomarker.score))*(100-1))/
                                                    (max(Biomarker.table$final.biomarker.score)-min(Biomarker.table$final.biomarker.score)))

      #add biomarker rank
      Biomarker.table$rank <- rank(-Biomarker.table$final.biomarker.score, ties.method = "min")

      #remove redundent columns
      Biomarker.table <- Biomarker.table[,c("final.biomarker.score", "rank")]

      #rename column names
      colnames(Biomarker.table) <- c("Score", "Rank")

      #filtering redundant (NaN) results
      Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]

      if(nrow(Biomarker.table)==0) {Biomarker.table <- NULL}

    }

    #ProgressBar: Preparation of the biomarker table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 100)
    }

    Results <- list("Driver table" = Driver.table,
                    "DE-mediator table" = DE.mediator.table,
                    "nonDE-mediators table" = non.DE.mediators.table,
                    "Biomarker table" = Biomarker.table)

    return(Results)
  }

  #=============================================================================
  #
  #    Code chunk 17: Assembling the differential/regression data in a dataframe
  #
  #=============================================================================

  #' Assembling the differential/regression data
  #'
  #' This function assembles a dataframe required for running the `ExIR` model. You may provide
  #' as many differential/regression data as you wish. Also, the datasets should be filtered
  #' beforehand according to your desired thresholds and, consequently, should only include the significant data.
  #' Each dataset provided should be a dataframe with one or two columns.
  #' The first column should always include differential/regression values
  #' and the second one (if provided) the significance values.
  #' @param ... Desired datasets/dataframes.
  #' @return A dataframe including the collective list of features in rows and all of the
  #' differential/regression data and their statistical significance in columns with the same
  #' order provided by the user.
  #' @aliases DDA
  #' @keywords diff.data.assembly
  #' @seealso \code{\link[influential]{exir}}
  #' @export
  #' @examples
  #' \dontrun{
  #' my.Diff_data <- diff.data.assembly(Differential_data1,
  #'                                    Differential_data2,
  #'                                    Regression_data1)
  #' }
  diff.data.assembly <- function(...) {

    #Getting the list of all datasets provided
    datasets <- lapply(list(...), as.data.frame)

    #Getting the feature names
    feature.names <- unique(unlist(lapply(X = datasets, FUN = rownames)))

    #Creating the Diff_data dataset
    Diff_data <- data.frame(Diff_value1 = rep(0,length(feature.names)),
                            row.names = feature.names)

    for (i in 1:length(datasets)) {

      feature.names.index <- match(rownames(datasets[[i]]),
                                   rownames(Diff_data))

      if(ncol(datasets[[i]]) == 2) {

        Diff_data[,paste("Diff_value", i, sep = "")] <- 0
        Diff_data[,paste("Diff_value", i, sep = "")][feature.names.index] <- datasets[[i]][,1]

        Diff_data[,paste("Sig_value", i, sep = "")] <- 1
        Diff_data[,paste("Sig_value", i, sep = "")][feature.names.index] <- datasets[[i]][,2]

      } else if(ncol(datasets[[i]]) == 1) {

        Diff_data[,paste("Diff_value", i, sep = "")] <- 0
        Diff_data[,paste("Diff_value", i, sep = "")][feature.names.index] <- datasets[[i]][,1]

      }
    }

    return(Diff_data)
  }
