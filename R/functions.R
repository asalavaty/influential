#=============================================================================
#
#    Code chunk 0: Influential package description
#
#=============================================================================

#' @keywords internal
#' @title Influential package
#' @description
#' The goal of \emph{\strong{\code{influential}}} is to help identification of the most influential nodes in a network.
#' This package contains functions for reconstruction of networks from adjacency matrices and
#' data frames, analysis of the topology of the network and calculation of centrality measures
#' as well as a novel and powerful influential node ranking. The first integrative method,
#' namely the \strong{Integrated Vector of Influence (IVI)}, that captures all topological dimensions
#' of the network for the identification of network most influential nodes is also provided as
#' a function. Also, neighborhood connectivity, H-index, local H-index, and collective
#' influence (CI), all of which required centrality measures for the calculation of IVI,
#' are for the first time provided in an R package. Furthermore, some functions have been
#' provided for the assessment of dependence and correlation of two network centrality
#' measures as well as the conditional probability of deviation from their corresponding
#' means in opposite directions.
#'
#' You may check the latest developmental version of the \emph{influential} package on its
#' \href{https://github.com/asalavaty/influential}{GitHub repository}
#'
#' @details
#' \itemize{
#'   \item Package: influential
#'   \item Type: Package
#'   \item Version: 1.0.0
#'   \item Date: 23-04-2020
#'   \item License: GPL-3
#' }
#'
#' @author
#' Author: Adrian (Abbas) Salavaty
#' Maintainer: Adrian (Abbas) Salavaty \email{abbas.salavaty@@gmail.com}
#'
#' You may find more information on my personal website at \href{https://www.abbassalavaty.com/}{www.AbbasSalavaty.com}
#'
#' @references
#' \itemize{
#'   \item Fred Viole and David Nawrocki (2013, ISBN:1490523995).
#'   \item Csardi G, Nepusz T (2006). “The igraph software package for complex network research.”
#' InterJournal, Complex Systems, 1695. http://igraph.org.
#' }
#' \strong{Note:} Adapted algorithms and sources are referenced in function document.
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

  # Getting the first neighbors of each node
  first.neighbors <- igraph::neighborhood(graph,
                                    nodes = vertices, mode = mode)

  # Getting the names of vertices with the order in use
  node.names <- sapply(first.neighbors, function(n) rownames(as.matrix(n[[]][1])))

  # Getting the neibors of each vertex
  node.neighbors <- sapply(first.neighbors, function(n) rownames(as.matrix(n[[]][-1])))

  # Getting the neighborhood size of each node
    first.neighbors.size <- sapply(node.neighbors,
                                   function(s) igraph::neighborhood.size(graph, s,
                                                                 mode = mode, order = 1) - 1)

  first.neighbors.size.sum <- sapply(first.neighbors.size, sum)

  # Calculation of neighborhood connectivity
  nc <- vector(mode = "numeric", length = length(first.neighbors.size))

  for (co in 1:length(first.neighbors.size.sum)) {
    nc[co] <- first.neighbors.size.sum[co]/length(unlist(node.neighbors[co]))
  }

  nc.table <- data.frame(Neighborhood_connectivity = nc)

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
  #' This function and all of its descriptions have been obtained from the centiserve package.
  #' For a complete description if the function and its arguments try this:
  #' ?centiserve::clusterrank
  #' @param graph The input graph as igraph object
  #' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices.
  #' @param directed Logical scalar, whether to directed graph is analyzed. This argument is ignored for undirected graphs.
  #' @param loops Logical; whether the loop edges are also counted.
  #' @return A numeric vector contaning the ClusterRank centrality scores for the selected vertices.
  #' @aliases CR
  #' @keywords clusterrank
  #' @family centrality functions
  #' @seealso \code{\link[centiserve]{clusterrank}} for a complete description on this function
  #' @export
  #' @examples
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' cr <- clusterrank(graph = My_graph, vids = GraphVertices, directed = FALSE, loops = TRUE)
  #' @importFrom centiserve clusterrank
  clusterrank <- centiserve::clusterrank

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
#' MyData <- centrality.measures
#' My.metrics.assessment <- double.cent.assess(data = MyData,
#'                                             nodes.colname = rownames(MyData),
#'                                             dependent.colname = "BC",
#'                                             independent.colname = "NC")
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
#' MyData <- centrality.measures
#' My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,
#'                                             nodes.colname = rownames(MyData),
#'                                             centrality1.colname = "BC",
#'                                             centrality2.colname = "NC")
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

#' Integrated Vector of Influence (IVI)
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
#' @seealso \code{\link[influential]{ivi}}
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

#' Integrated Vector of Influence (IVI)
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
#' @keywords IVI integrated_vector_of_influence
#' @seealso \code{\link[influential]{ivi.from.indices}}
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
  CR <- clusterrank(graph = graph, vids = vertices, directed = directed, loops = loops)
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

  CR <- clusterrank(graph = graph, vids = vertices, directed = directed, loops = loops)
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
