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
#' Also, a web-based \href{https://influential.erc.monash.edu/}{Influential Software Package} with a convenient
#' user-interface (UI) has been developed for the comfort of all users including those without a coding background.
#'
#' @details
#' \itemize{
#'   \item Package: influential
#'   \item Type: Package
#'   \item Version: 2.2.9.9000
#'   \item Date: 01-12-2023
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
#' You may find more information on my personal website at \href{https://asalavaty.com/}{www.ASalavaty.com}
#'
#' @references
#' \itemize{
#'   \item Fred Viole and David Nawrocki (2013, ISBN:1490523995).
#'   \item Csardi G, Nepusz T (2006). “The igraph software package for complex network research.”
#' InterJournal, Complex Systems, 1695. \url{https://igraph.org/}.
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
#    Code chunk 0.5: Run shiny apps (examples)
#
#=============================================================================

#' @title Run shiny app
#' @description Run shiny apps included in the influential R package.
#' Also, a web-based \href{https://influential.erc.monash.edu/}{Influential Software Package} with a convenient
#' user-interface (UI) has been developed for the comfort of all users including those without a coding background.
#' @param shinyApp The name of the shiny app you want to run. You can get the exact name of the available
#' shiny apps via the following command.
#' \emph{list.files(system.file("ShinyApps", package = "influential"))}. Please also note this function is
#' case-sensitive.
#' @return A shiny app.
#' @keywords runShinyApp
#' @export runShinyApp
#' @examples
#' \dontrun{
#' runShinyApp(shinyApp = "IVI")
#' }
runShinyApp <- function(shinyApp) {
  
  cat("Checking the requirements for running the shiny app.\n")
  
  # Check if the BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager", dependencies = TRUE)
  }
  
  # Loop through each package
  for (pkg in c("shiny", "shinythemes", "shinyWidgets", "shinyjs",
                "shinycssloaders", "colourpicker", "DT", "magrittr", "janitor",
                "ranger", "coop", "influential", "ggplot2", "igraph")) {
    # Check if the package namespace is available
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # Install the package if it's not available
      BiocManager::install(pkg, update = TRUE, ask = FALSE, quiet = TRUE)
        cat(paste0("The package '", pkg, "' required for running the shiny app is installed.\n"))
    }
  }
  
    # locate all the shiny app examples that exist
    validExamples <- list.files(system.file("ShinyApps", package = "influential"))

    validExamplesMsg <-
      paste0(
        "Valid shiny apps are: '",
        paste(validExamples, collapse = "', '"),
        "'")

    # if an invalid shiny app is given, throw an error
    if (missing(shinyApp) || !nzchar(shinyApp) ||
        !shinyApp %in% validExamples) {
      stop(
        'Please run `influential::runShinyApp()` with a valid shiny app name as an argument.\n',
        validExamplesMsg,
        call. = FALSE)
    }
    
    cat("Loading the shiny app ...")

    # find and launch the app
    appDir <- system.file("ShinyApps", shinyApp, package = "influential")
    shiny::runApp(appDir, display.mode = "normal")
}

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
#' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
#' @return A vector including the neighborhood connectivity score of each vertex inputted.
#' @aliases NC
#' @keywords neighborhood_connectivity
#' @family centrality functions
#' @seealso \code{\link[influential]{ivi}},
#' \code{\link[influential]{cent_network.vis}}
#' @export neighborhood.connectivity
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' neighrhood.co <- neighborhood.connectivity(graph = My_graph,
#'                                            vertices = GraphVertices,
#'                                            mode = "all")
#'                                            }
neighborhood.connectivity <- function(graph, vertices = V(graph), mode = "all", verbose = FALSE) {

  # Getting the names of vertices
  if(inherits(vertices, "igraph.vs")) {
    node.names <- as.character(igraph::as_ids(vertices))
  } else {
    node.names <- as.character(vertices)
  }

  if(verbose) {
    cat("Getting the first neighbors of each node\n")
  }
  # Getting the first neighbors of each node
  node.neighbors <- sapply(as.list(node.names),
                           FUN = function(i) as.character(igraph::as_ids(igraph::neighbors(graph = graph,
                                                                                           v = i,
                                                                                           mode = mode))))
  if(verbose) {
    cat("Getting the neighborhood size of each node\n")
  }
  # Getting the neighborhood size of each node
    first.neighbors.size <- sapply(node.neighbors,
                                   function(s) igraph::neighborhood.size(graph = graph,
                                                                         nodes = s,
                                                                         mode = mode,
                                                                         order = 1) - 1)

    if(verbose) {
      cat("Calculating the neighborhood connectivity of nodes\n")
    }
  # Calculation of neighborhood connectivity

  if(length(vertices) == 1) {
    first.neighbors.size.sum <- sum(first.neighbors.size)
    temp.nc <- first.neighbors.size.sum/nrow(node.neighbors)
  } else {

    first.neighbors.size.sum <- sapply(first.neighbors.size, sum)
    temp.nc <- vector(mode = "numeric", length = length(vertices))

    for (i in 1:length(vertices)) {
      temp.nc[i] <- first.neighbors.size.sum[i]/length(node.neighbors[[i]])
    }
  }

  temp.nc[c(which(is.nan(temp.nc)), which(is.na(temp.nc)))] <- 0

  names(temp.nc) <- node.names

  return(temp.nc)

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
#' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
#' @return A vector including the H-index of each vertex inputted.
#' @aliases h.index
#' @keywords h_index
#' @family centrality functions
#' @seealso \code{\link[influential]{ivi}},
#' \code{\link[influential]{cent_network.vis}}
#' @export h_index
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' h.index <- h_index(graph = My_graph, vertices = GraphVertices, mode = "all")
#' }
#' @importFrom utils tail
  h_index <- function(graph, vertices = V(graph), mode = "all", verbose = FALSE) {

    if(verbose) {
      cat("Getting the first neighbors of each node\n")
    }
  # Getting the first neighbors of each node
  first.neighbors <- igraph::neighborhood(graph, nodes = vertices, mode = mode)

  if(verbose) {
    cat("Getting the neighbors of each node\n")
  }
  # Getting the neighbors of each node
  node.neighbors <- sapply(first.neighbors, function(n) rownames(as.matrix(n[[]][-1])), simplify = F)

  if(verbose) {
    cat("Getting the neighborhood size of each node\n")
  }
  # Getting the neighborhood size of each node
  first.neighbors.size <- lapply(node.neighbors, function(s) igraph::neighborhood.size(graph, s,
                                                                               mode = mode, order = 1) - 1)
  if(verbose) {
    cat("Calculating the H-index\n")
  }
  # Calculating the H-index
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

  # Getting the names of vertices
  if(inherits(vertices, "igraph.vs")) {
    node.names <- as.character(igraph::as_ids(vertices))
  } else {
    node.names <- as.character(vertices)
  }
  names(hindex) <- node.names
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
  #' @param ncores Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of 
  #' cores to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").
  #' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
  #' @return A vector including the local H-index of each vertex inputted.
  #' @aliases lh.index
  #' @keywords lh_index
  #' @family centrality functions
  #' @seealso \code{\link[influential]{ivi}},
  #' \code{\link[influential]{cent_network.vis}}
  #' @export lh_index
  #' @examples
  #' \dontrun{
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' lh.index <- lh_index(graph = My_graph, vertices = GraphVertices, mode = "all", ncores = 1)
  #' }
  #' @importFrom foreach %dopar%
  lh_index <- function(graph, vertices = V(graph), mode = "all", ncores = "default", verbose = FALSE) {
    
    # Make clusters for parallel processing
    cl <- parallel::makeCluster(ifelse(ncores == "default", parallel::detectCores() - 1, ncores))
    doParallel::registerDoParallel(cl)

    if(verbose) {
      cat("Getting the first neighbors of each node\n")
    }
    # Getting the first neighbors of each node
    first.neighbors <- igraph::neighborhood(graph, nodes = vertices, mode = mode)

    if(verbose) {
      cat("Calculating H-index\n")
    }
    # Calculation of local H-index (LH-index)
    lhindex <- foreach::foreach(i = 1:length(vertices), .combine = "c", .multicombine = TRUE, 
                                .packages = c("igraph", "influential")) %dopar% {
                                  sum(h_index(graph = graph, vertices = unlist(first.neighbors[i]), mode = mode, verbose = FALSE))
                                }
    parallel::stopCluster(cl)
    
    if(verbose) {
      cat("Preparing the LH-index\n")
    }
    
    # Getting the names of vertices
    if(inherits(vertices, "igraph.vs")) {
      node.names <- as.character(igraph::as_ids(vertices))
    } else {
      node.names <- as.character(vertices)
    }
    names(lhindex) <- node.names
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
  #' it is the product of the reduced degree (degree - 1) of a node and the total (sum of) reduced
  #' degrees of all nodes at a distance d from the node.
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
  #' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
  #' @return A vector of collective influence for each vertex of the graph corresponding to
  #' the order of vertices output by V(graph).
  #' @aliases CI
  #' @keywords collective.influence
  #' @family centrality functions
  #' @seealso \code{\link[influential]{ivi}},
  #' \code{\link[influential]{cent_network.vis}}
  #' @export collective.influence
  #' @examples
  #' \dontrun{
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' ci <- collective.influence(graph = My_graph, vertices = GraphVertices, mode = "all", d=3)
  #' }
  collective.influence <- function(graph, vertices = V(graph), mode = "all", d=3, verbose = FALSE) {

    ci <- vector(mode="numeric", length=length(vertices))  # collective influence output

    if(verbose) {
      cat("Calculating the reduced degrees of nodes\n")
    }
    
    # Calculate the reduced degree of nodes
    reduced.degrees <- igraph::degree(graph = graph,
                              v = vertices,
                              mode = mode) - 1

    if(verbose) {
      cat("Identifing only neighbors at distance d\n")
    }
    
    # Identify only neighbors at distance d
    nodes.at.distance <- igraph::neighborhood(graph = graph, nodes = vertices,
                                              mode = mode, order=d, mindist=d)
    
    # Get the non-duplicated vector of node names of neighbors at distance d
    nodes.at.distance_names <- base::unique(base::unlist(base::sapply(nodes.at.distance, igraph::as_ids)))
    
    if(verbose) {
      cat("Calculating the reduced degrees of neighbors at distance d\n")
    }
    
    # Calculate the reduced degree of neighbors at distance d
    nodes.at.distance_reduced.degrees <- igraph::degree(graph = graph,
                                                v = nodes.at.distance_names,
                                                mode = mode) - 1
    
    if(verbose) {
      cat("Calculating the the Collective Influence\n")
    }
    
    # Calculate the collective influence
    for (i in 1:length(nodes.at.distance)) {
      rd <- reduced.degrees[i]  # i is the index of the node
      rd.neighbours <- sum(nodes.at.distance_reduced.degrees[igraph::as_ids(nodes.at.distance[[i]])]) # calculate the cumulative reduced degree of neighbors
      ci[i] <- rd * rd.neighbours
    }

    # Getting the names of vertices
    if(inherits(vertices, "igraph.vs")) {
      node.names <- as.character(igraph::as_ids(vertices))
    } else {
      node.names <- as.character(vertices)
    }
    names(ci) <- node.names
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
  #' @param ncores Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of cores 
  #' to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").
  #' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
  #' @return A numeric vector contaning the ClusterRank centrality scores for the selected vertices.
  #' @aliases CR
  #' @keywords clusterRank
  #' @family centrality functions
  #' @seealso \code{\link[influential]{ivi}},
  #' \code{\link[influential]{cent_network.vis}}
  #' @export clusterRank
  #' @examples
  #' \dontrun{
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' cr <- clusterRank(graph = My_graph, vids = GraphVertices, 
  #' directed = FALSE, loops = TRUE, ncores = 1)
  #' }
  #' @importFrom foreach %dopar%
  clusterRank <- function(graph, vids = V(graph),
                        directed = FALSE, loops = TRUE, ncores = "default", verbose = FALSE) {
    
    # Make clusters for parallel processing
    cl <- parallel::makeCluster(ifelse(ncores == "default", parallel::detectCores() - 1, ncores))
    doParallel::registerDoParallel(cl)

  vertex.transitivity <- vector(mode = "numeric")
  
  if(verbose) {
    cat("Calculating the transitivity of nodes\n")
  }

  if (directed) {
    cl.Rank.mode <- "out"
    foreach::foreach(i = V(graph), .combine = c) %dopar% {
      vertex.neighborhood <- igraph::neighborhood(graph = graph,
                                                  order = 1, nodes = i,
                                                  mode = cl.Rank.mode)[[1]][-1]
      if (length(vertex.neighborhood) < 2) {
        NaN
      } else {
        indc.subgraph <- igraph::induced.subgraph(graph = graph, vertex.neighborhood)
        igraph::ecount(indc.subgraph)/(igraph::vcount(indc.subgraph)*(igraph::vcount(indc.subgraph)-1))
      }
    } -> vertex.transitivity
  } else {
    cl.Rank.mode <- "all"
    vertex.transitivity <- igraph::transitivity(graph = graph, type = "local")
  }
  parallel::stopCluster(cl)

  if(inherits(vids, "igraph.vs")) {
    vertices.index <- stats::na.omit(match(vids, V(graph)))
  } else {
    vertices.index <- stats::na.omit(match(vids, igraph::as_ids(V(graph))))
  }
  
  cl <- parallel::makeCluster(ifelse(ncores == "default", parallel::detectCores() - 1, ncores))
  doParallel::registerDoParallel(cl)

  if(verbose) {
    cat("Getting the neighborhood of selected nodes and calculating the ClusterRank\n")
  }

  cl.Rank <- foreach::foreach(i = V(graph)[vertices.index], .combine = c, .multicombine = TRUE,
                              .packages = "igraph") %dopar% {
                                if (is.nan(vertex.transitivity[i])) {
                                  NaN
                                } else {
                                  selected.v.neighborhood <- igraph::neighborhood(graph = graph,
                                                                                  order = 1, nodes = i,
                                                                                  mode = cl.Rank.mode)[[1]][-1]
                                  temp.cl.Rank <- 0
                                  for (j in selected.v.neighborhood) {
                                    temp.cl.Rank <- temp.cl.Rank + igraph::degree(graph = graph,
                                                                                  v = j, mode = cl.Rank.mode,
                                                                                  loops = loops) + 1
                                  }
                                  temp.cl.Rank * vertex.transitivity[i]
                                }
                              }
  
  parallel::stopCluster(cl)

  if (igraph::is.named(graph)) {
    names(cl.Rank) <- igraph::V(graph)$name[vertices.index]
  }

  return(cl.Rank)
  }

#=============================================================================
#
#    Code chunk 6: Calculation of the conditional probability of deviation from
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
#' @export cond.prob.analysis
#' @examples
#' \dontrun{
#' MyData <- centrality.measures
#' My.conditional.prob <- cond.prob.analysis(data = MyData,
#'                                           nodes.colname = rownames(MyData),
#'                                           Desired.colname = "BC",
#'                                           Condition.colname = "NC")
#'                                           }
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
#' @export double.cent.assess
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
#' @export double.cent.assess.noRegression
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
#' A shiny app has also been developed for the calculation of IVI as well as IVI-based network
#' visualization, which is accessible using the `influential::runShinyApp("IVI")` command.
#' You can also access the shiny app online at https://influential.erc.monash.edu/.
#' @param DC A vector containing the values of degree centrality of the desired vertices.
#' @param CR A vector containing the values of ClusterRank of the desired vertices.
#' @param LH_index A vector containing the values of local H-index of the desired vertices.
#' @param NC A vector containing the values of neighborhood connectivity of the desired vertices.
#' @param BC A vector containing the values of betweenness centrality of the desired vertices.
#' @param CI A vector containing the values of Collective Influence of the desired vertices.
#' @param scale Character string; the method used for scaling/normalizing the results. Options include 'range' (normalization within a 1-100 range), 
#' 'z-scale' (standardization using the z-score), and 'none' (no data scaling). The default selection is 'range'. Opting for the 'range' method is 
#' suitable when exploring a single network, allowing you to observe the complete spectrum and distribution of node influences. In this case, there is 
#' no intention to establish a specific threshold for the outcomes. However, it is possible to identify and present the top influential nodes 
#' based on their rankings. Conversely, the 'z-scale' option proves advantageous if the aim is to compare node influences across multiple networks or 
#' if there is a desire to establish a threshold (usually z-score > 1.645) for generating a list of the most influential nodes without manual intervention.
#' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
#' @return A numeric vector with the IVI values based on the provided centrality measures.
#' @aliases IVI.FI
#' @keywords ivi.from.indices
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{cent_network.vis}}
#' @export ivi.from.indices
#' @examples
#' \dontrun{
#' MyData <- centrality.measures
#' My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,
#'                                     CR = centrality.measures$CR,
#'                                     NC = centrality.measures$NC,
#'                                     LH_index = centrality.measures$LH_index,
#'                                     BC = centrality.measures$BC,
#'                                     CI = centrality.measures$CI)
#'                                     }
ivi.from.indices <- function(DC, CR, LH_index, NC, BC, CI, scale = "range", verbose = FALSE) {

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
  
  if(verbose) {
    cat("1-100 normalization of centrality measures\n")
  }

  #1-100 normalization of centrality measures

  if(length(temp.DC) > 1 & any(temp.DC > 0)) {
    temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
  }
  
  if(length(temp.CR) > 1 & any(temp.CR > 0)) {
    temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
  }
  
  if(length(temp.LH_index) > 1 & any(temp.LH_index > 0)) {
    temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
  }
  
  if(length(temp.NC) > 1 & any(temp.NC > 0)) {
    temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
  }
  
  if(length(temp.BC) > 1 & any(temp.BC > 0)) {
    temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
  }
  
  if(length(temp.CI) > 1 & any(temp.CI > 0)) {
    temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
  }
  
  #Calculation of IVI
  
  if(verbose) {
    cat("Calculating the Spreading Rank\n")
  }

  spreading.rank <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))
  
  suppressWarnings(
    if(any(stats::na.omit(spreading.rank) == 0 | is.na(spreading.rank))) {
      spreading.rank[which(spreading.rank == 0 | is.na(spreading.rank))] <- 1
    }
  )
  
  if(verbose) {
    cat("Calculating the Hubness Rank\n")
  }

  hubness.rank <- (temp.DC+temp.LH_index)
  
  suppressWarnings(
    if(any(stats::na.omit(hubness.rank) == 0 | is.na(hubness.rank))) {
      hubness.rank[which(hubness.rank == 0 | is.na(hubness.rank))] <- 1
    }
  )
  
  
  if(verbose) {
    cat("Calculating the IVI\n")
  }
  
  temp.ivi <- (hubness.rank)*(spreading.rank)

  #1-100 normalization of IVI

  if(scale == "range") {
    
    if(verbose) {
      cat("1-100 normalization of IVI\n")
    }

    if(length(unique(temp.ivi)) > 1) {
      temp.ivi <- 1+(((temp.ivi-min(temp.ivi))*(100-1))/(max(temp.ivi)-min(temp.ivi)))
    }
  } else if(scale == 'z-scale') {
    if(verbose) {
      cat("Z-score standardization of IVI\n")
    }
    
    temp.ivi <- base::scale(temp.ivi)
    temp.ivi.names <- rownames(temp.ivi)
    temp.ivi <- c(temp.ivi)
    names(temp.ivi) <- temp.ivi.names
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
#' #' A shiny app has also been developed for the calculation of IVI as well as IVI-based network
#' visualization, which is accessible using the `influential::runShinyApp("IVI")` command.
#' You can also access the shiny app online at https://influential.erc.monash.edu/.
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
#' @param scale Character string; the method used for scaling/normalizing the results. Options include 'range' (normalization within a 1-100 range), 
#' 'z-scale' (standardization using the z-score), and 'none' (no data scaling). The default selection is 'range'. Opting for the 'range' method is 
#' suitable when exploring a single network, allowing you to observe the complete spectrum and distribution of node influences. In this case, there is 
#' no intention to establish a specific threshold for the outcomes. However, it is possible to identify and present the top influential nodes 
#' based on their rankings. Conversely, the 'z-scale' option proves advantageous if the aim is to compare node influences across multiple networks or 
#' if there is a desire to establish a threshold (usually z-score > 1.645) for generating a list of the most influential nodes without manual intervention.
#' @param ncores Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of 
#' cores to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").
#' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
#' @return A numeric vector with the IVI values based on the provided centrality measures.
#' @aliases IVI
#' @keywords IVI integrated_value_of_influence
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{cent_network.vis}}
#' @export ivi
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' My.vertices.IVI <- ivi(graph = My_graph, vertices = GraphVertices,
#'                        weights = NULL, directed = FALSE, mode = "all",
#'                        loops = TRUE, d = 3, scale = "range")
#' }
ivi <- function(graph, vertices = V(graph), weights = NULL, directed = FALSE,
                mode = "all", loops = TRUE, d = 3, scale = "range", ncores = "default", verbose = FALSE) {

  #Calculation of required centrality measures
  
  if(verbose) {
    message("Calculating the Degree Centrality of Nodes\n")
  }

  DC <- igraph::degree(graph = graph, v = vertices, mode = mode, loops = loops)
  
  if(verbose) {
    message("Calculating the ClusterRank of Nodes\n")
  }
  
  CR <- clusterRank(graph = graph, vids = vertices, directed = directed, loops = loops, ncores = ncores, verbose = verbose)
  
  if(verbose) {
    message("Calculating the Local H-index of Nodes\n")
  }
  
  LH_index <- lh_index(graph = graph, vertices = vertices, mode = mode, ncores = ncores, verbose = verbose)
  
  if(verbose) {
    message("Calculating the Neighborhood Connectivity of Nodes\n")
  }
  
  NC <- neighborhood.connectivity(graph = graph, vertices = vertices, mode = mode, verbose = verbose)
  
  if(verbose) {
    message("Calculating the Betweenness Centrality of Nodes\n")
  }
  
  BC <- betweenness(graph = graph, v = vertices, directed = directed, weights = weights)
  
  if(verbose) {
    message("Calculating the Collective Influence of Nodes\n")
  }
  
  CI <- collective.influence(graph = graph, vertices = vertices, mode = mode, d = d, verbose = verbose)

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
  
  if(verbose) {
    cat("1-100 normalization of centrality measures\n")
  }

  #1-100 normalization of centrality measures

  if(length(temp.DC) > 1 & any(temp.DC > 0)) {
    temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
  }
  
  if(length(temp.CR) > 1 & any(temp.CR > 0)) {
    temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
  }
  
  if(length(temp.LH_index) > 1 & any(temp.LH_index > 0)) {
    temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
  }
  
  if(length(temp.NC) > 1 & any(temp.NC > 0)) {
    temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
  }
  
  if(length(temp.BC) > 1 & any(temp.BC > 0)) {
    temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
  }
  
  if(length(temp.CI) > 1 & any(temp.CI > 0)) {
    temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
  }

  #Calculation of IVI
  
  if(verbose) {
    message("Calculation of IVI\n")
  }
  
  if(verbose) {
    cat("Calculating the Spreading Rank\n")
  }

  spreading.rank <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))

  suppressWarnings(
    if(any(stats::na.omit(spreading.rank) == 0 | is.na(spreading.rank))) {
      spreading.rank[which(spreading.rank == 0 | is.na(spreading.rank))] <- 1
    }
  )
  
  if(verbose) {
    cat("Calculating the Hubness Rank\n")
  }

  hubness.rank <- (temp.DC+temp.LH_index)

  suppressWarnings(
    if(any(stats::na.omit(hubness.rank) == 0 | is.na(hubness.rank))) {
      hubness.rank[which(hubness.rank == 0 | is.na(hubness.rank))] <- 1
    }
  )
  

  temp.ivi <- (hubness.rank)*(spreading.rank)

  #1-100 normalization of IVI

  if(scale == "range") {
    
    if(verbose) {
      cat("1-100 normalization of IVI\n")
    }

    if(length(unique(temp.ivi)) > 1) {
      temp.ivi <- 1+(((temp.ivi-min(temp.ivi))*(100-1))/(max(temp.ivi)-min(temp.ivi)))
    }
  } else if(scale == 'z-scale') {
    if(verbose) {
      cat("Z-score standardization of IVI\n")
    }
    
    temp.ivi <- base::scale(temp.ivi)
    temp.ivi.names <- rownames(temp.ivi)
    temp.ivi <- c(temp.ivi)
    names(temp.ivi) <- temp.ivi.names
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
#' of nodes as a requirement for calculation of spreading score. If the graph has a weight edge attribute,
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
#' @param scale Character string; the method used for scaling/normalizing the results. Options include 'range' (normalization within a 1-100 range), 
#' 'z-scale' (standardization using the z-score), and 'none' (no data scaling). The default selection is 'range'. Opting for the 'range' method is 
#' suitable when exploring a single network, allowing you to observe the complete spectrum and distribution of node influences. In this case, there is 
#' no intention to establish a specific threshold for the outcomes. However, it is possible to identify and present the top spreading nodes 
#' based on their rankings. Conversely, the 'z-scale' option proves advantageous if the aim is to compare node influences across multiple networks or 
#' if there is a desire to establish a threshold (usually z-score > 1.645) for generating a list of the most spreading nodes without manual intervention.
#' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
#' @return A numeric vector with Spreading scores.
#' @keywords spreading.score
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{cent_network.vis}}
#' @export spreading.score
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' Spreading.score <- spreading.score(graph = My_graph, vertices = GraphVertices,
#'                                    weights = NULL, directed = FALSE, mode = "all",
#'                                    loops = TRUE, d = 3, scale = "range")
#' }
spreading.score <- function(graph, vertices = V(graph), weights = NULL, directed = FALSE,
                            mode = "all", loops = TRUE, d = 3, scale = "range", verbose = FALSE) {

  #Calculation of required centrality measures
  
  if(verbose) {
    message("Calculating the ClusterRank of Nodes\n")
  }

  CR <- clusterRank(graph = graph, vids = vertices, directed = directed, loops = loops, verbose = verbose)
  CR[which(is.nan(CR))] <- 0
  
  if(verbose) {
    message("Calculating the Neighborhood Connectivity of Nodes\n")
  }
  
  NC <- neighborhood.connectivity(graph = graph, vertices = vertices, mode = mode, verbose = verbose)
  
  if(verbose) {
    message("Calculating the Betweenness Centrality of Nodes\n")
  }
  
  BC <- betweenness(graph = graph, v = vertices, directed = directed, weights = weights)
  
  if(verbose) {
    message("Calculating the Collective Influence of Nodes\n")
  }
  
  CI <- collective.influence(graph = graph, vertices = vertices, mode = mode, d = d, verbose = verbose)

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
  
  if(verbose) {
    cat("1-100 normalization of centrality measures\n")
  }

  #1-100 normalization of centrality measures

  if(length(temp.CR) > 1 & any(temp.CR > 0)) {
    temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
  }
  
  if(length(temp.NC) > 1 & any(temp.NC > 0)) {
    temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
  }
  
  if(length(temp.BC) > 1 & any(temp.BC > 0)) {
    temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
  }
  
  if(length(temp.CI) > 1 & any(temp.CI > 0)) {
    temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
  }

  #Calculation of spreading.score

  temp.spreading.score <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))

  #1-100 normalization of spreading score

  if(scale == 'range') {
    
    if(verbose) {
      cat("1-100 normalization of Spreading Score\n")
    }

    if(length(unique(temp.spreading.score)) > 1) {
    temp.spreading.score <- 1+(((temp.spreading.score-min(temp.spreading.score))*(100-1))/(max(temp.spreading.score)-min(temp.spreading.score)))
    }
  } else if(scale == 'z-scale') {
    if(verbose) {
      cat("Z-score standardization of Spreading Score\n")
    }
    
    temp.spreading.score <- base::scale(temp.spreading.score)
    temp.spreading.score.names <- rownames(temp.spreading.score)
    temp.spreading.score <- c(temp.spreading.score)
    names(temp.spreading.score) <- temp.spreading.score.names
    
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
#' @param scale Character string; the method used for scaling/normalizing the results. Options include 'range' (normalization within a 1-100 range), 
#' 'z-scale' (standardization using the z-score), and 'none' (no data scaling). The default selection is 'range'. Opting for the 'range' method is 
#' suitable when exploring a single network, allowing you to observe the complete spectrum and distribution of node influences. In this case, there is 
#' no intention to establish a specific threshold for the outcomes. However, it is possible to identify and present the top hub nodes 
#' based on their rankings. Conversely, the 'z-scale' option proves advantageous if the aim is to compare node influences across multiple networks or 
#' if there is a desire to establish a threshold (usually z-score > 1.645) for generating a list of the most hub nodes without manual intervention.
#' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
#' @return A numeric vector with the Hubness scores.
#' @keywords hubness.score
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{cent_network.vis}}
#' @export hubness.score
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' GraphVertices <- V(My_graph)
#' Hubness.score <- hubness.score(graph = My_graph, vertices = GraphVertices,
#'                                directed = FALSE, mode = "all",
#'                                loops = TRUE, scale = "range")
#' }
hubness.score <- function(graph, vertices = V(graph), directed = FALSE,
                          mode = "all", loops = TRUE, scale = "range", verbose = FALSE) {

  #Calculation of required centrality measures
  
  if(verbose) {
    message("Calculating the Degree Centrality of Nodes\n")
  }

  DC <- igraph::degree(graph = graph, v = vertices, mode = mode, loops = loops)
  
  if(verbose) {
    message("Calculating the Local H-index of Nodes\n")
  }
  
  LH_index <- lh_index(graph = graph, vertices = vertices, mode = mode, verbose = verbose)

  #Generating temporary measures

  temp.DC <- DC
  temp.LH_index <- LH_index

  #Removing the NAN and NA values

  temp.DC[c(which(is.nan(temp.DC)), which(is.na(temp.DC)))] <- 0
  temp.LH_index[c(which(is.nan(temp.LH_index)), which(is.na(temp.LH_index)))] <- 0
  
  if(verbose) {
    cat("1-100 normalization of centrality measures\n")
  }

  #1-100 normalization of centrality measures

  if(length(temp.DC) > 1 & any(temp.DC > 0)) {
    temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
  }
  
  if(length(temp.LH_index) > 1 & any(temp.LH_index > 0)) {
    temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
  }

  #Calculation of Hubness score

  temp.hubness.score <- (temp.DC+temp.LH_index)

  #1-100 normalization of Hubness score

  if(scale == "range") {
    
    if(verbose) {
      cat("1-100 normalization of Hubness Score\n")
    }

    if(length(unique(temp.hubness.score)) > 1) {
    temp.hubness.score <- 1+(((temp.hubness.score-min(temp.hubness.score))*(100-1))/(max(temp.hubness.score)-min(temp.hubness.score)))
    
    }
  } else if(scale == 'z-scale') {
    if(verbose) {
      cat("Z-score standardization of Hubness Score\n")
    }
    
    temp.hubness.score <- base::scale(temp.hubness.score)
    temp.hubness.score.names <- rownames(temp.hubness.score)
    temp.hubness.score <- c(temp.hubness.score)
    names(temp.hubness.score) <- temp.hubness.score.names
    
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
#' performance of a novel algorithm in identification of network influential nodes by considering
#' the SIRIR ranks as the ground truth (gold standard).
#' @param graph A graph (network) of the igraph class.
#' @param vertices A vector of desired vertices, which could be obtained by the V function.
#' @param beta Non-negative scalar. The rate of infection of an individual that is susceptible
#' and has a single infected neighbor. The infection rate of a susceptible individual with n
#' infected neighbors is n times beta. Formally this is the rate parameter of an exponential
#' distribution.
#' @param gamma Positive scalar. The rate of recovery of an infected individual.
#' Formally, this is the rate parameter of an exponential distribution.
#' @param no.sim Integer scalar, the number of simulation runs to perform SIR model on the
#' original network as well as perturbed networks generated by leave-one-out technique.
#' You may choose a different no.sim based on the available memory on your system.
#' @param ncores Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of 
#' cores to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").
#' @param seed A single value, interpreted as an integer to be used for random number generation.
#' @param node_verbose Logical; whether the process of Parallel Socket Cluster creation should be printed (default is FALSE).
#' @param loop_verbose Logical; whether the accomplishment of the evaluation of network nodes in each loop should be printed (default is TRUE).
#' @return A two-column dataframe; a column containing the difference values of the original and
#' perturbed networks and a column containing node influence rankings
#' @aliases SIRIR
#' @keywords sirir
#' @family centrality functions
#' @seealso \code{\link[influential]{cent_network.vis}},
#' and \code{\link[igraph]{sir}} for a complete description on SIR model
#' @export sirir
#' @examples
#' \dontrun{
#' set.seed(1234)
#' My_graph <- igraph::sample_gnp(n=50, p=0.05)
#' GraphVertices <- V(My_graph)
#' Influence.Ranks <- sirir(graph = My_graph, vertices = GraphVertices, 
#'                          beta = 0.5, gamma = 1, ncores = "default", no.sim = 10, seed = 1234)
#' }
#' @importFrom igraph vcount as_ids sir
#' @importFrom foreach %dopar%
sirir <- function(graph, vertices = V(graph),
                  beta = 0.5, gamma = 1, no.sim = 100,  
                  ncores = "default", seed = 1234, loop_verbose = TRUE, node_verbose = FALSE) {
  
  suppressWarnings({
    
    # Make clusters for parallel processing
    cl <- parallel::makeCluster(ifelse(ncores == "default", parallel::detectCores() - 1, ncores), 
                                outfile=ifelse(node_verbose, "", "NULL"))
    doParallel::registerDoParallel(cl)
    
    #Define a data frame
    temp.loocr.table <- data.frame(difference.value = vector("numeric", length = length(vertices)),
                                   rank = vector("integer", length = length(vertices)))
    
    if(inherits(vertices, "character")) {
      rownames(temp.loocr.table) <- vertices
    } else if(inherits(vertices, "igraph.vs")) {
      rownames(temp.loocr.table) <- igraph::as_ids(vertices)
    }
    
    #Model the spreading based on all nodes
    set.seed(seed)
    all.included.spread <- igraph::sir(graph = graph, beta = beta,
                                       gamma = gamma, no.sim = no.sim)
    
    #Getting the mean of spread in each independent experiment
    all.mean.spread <- sapply(1:length(all.included.spread), function(i) {
      max(all.included.spread[[i]]$NR)
    })
    
    all.mean.spread <- mean(all.mean.spread)
    
    #Model the spread based on leave one out cross ranking (LOOCR)
    loocr.mean.spread_vec <- 
      foreach::foreach(s = vertices,
                       .combine = "c", 
                       .verbose = loop_verbose,
                       .multicombine = TRUE,
                       .packages = c("igraph", "influential"),
                       .export = c("graph", "beta", "gamma", "no.sim", "all.mean.spread")
      ) %dopar% {
        
        temp.graph <- igraph::delete_vertices(graph, s)
        set.seed(seed)
        loocr.spread <- igraph::sir(graph = temp.graph, beta = beta,
                                    gamma = gamma, no.sim = no.sim)
        
        loocr.mean.spread <- sapply(1:length(loocr.spread), function(h) {
          max(loocr.spread[[h]]$NR)
        })
        
        cat(paste("\nProcessing vertex ", as_ids(s), " is done!", "\n", sep = ""))  # Print message for each iteration
        
        # return mean of loocr spreads
        mean(loocr.mean.spread)
      }
    
    # Stop the parallel backend
    parallel::stopCluster(cl)
    
    temp.loocr.table$difference.value <- all.mean.spread - loocr.mean.spread_vec
    temp.loocr.table$rank <- rank(-temp.loocr.table$difference.value, ties.method = "min")
    
  })
  
  return(temp.loocr.table)
}

#=============================================================================
#
#    Code chunk 14: Some required functions from the igraph package
#
#=============================================================================

  #' @importFrom igraph graph_from_data_frame
  #' @export
  igraph::graph_from_data_frame


  #*****************************************************************#

  #' @importFrom igraph graph_from_adjacency_matrix
  #' @export
  igraph::graph_from_adjacency_matrix

  #*****************************************************************#

  #' @importFrom igraph graph_from_incidence_matrix
  #' @export
  igraph::graph_from_incidence_matrix

  #*****************************************************************#

  #' @importFrom igraph V
  #' @export
  igraph::V

  #*****************************************************************#

  #' Vertex betweenness centrality
  #'
  #' This function and all of its descriptions have been obtained from the igraph package.
  #' @param graph The graph to analyze (an igraph graph).
  #' @param v The vertices for which the vertex betweenness will be calculated.
  #' @param directed Logical, whether directed paths should be considered while determining the shortest paths.
  #' @param weights Optional positive weight vector for calculating weighted betweenness.
  #' If the graph has a weight edge attribute, then this is used by default. Weights are used to calculate weighted shortest paths, so they are interpreted as distances.
  #' @param normalized Logical scalar, whether to normalize the betweenness scores. If TRUE, then the results are normalized.
  #' @param ... Additional arguments according to the original \code{\link[igraph]{betweenness}} function in the package igraph.
  #' @return A numeric vector with the betweenness score for each vertex in v.
  #' @aliases BC
  #' @keywords betweenness_centrality
  #' @family centrality functions
  #' @seealso \code{\link[influential]{ivi}},
  #' \code{\link[influential]{cent_network.vis}},
  #' and \code{\link[igraph]{betweenness}} for a complete description on this function
  #' @export betweenness
  #' @examples
  #' \dontrun{
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' GraphVertices <- V(My_graph)
  #' My_graph_betweenness <- betweenness(My_graph, v = GraphVertices,
  #'                         directed = FALSE, normalized = FALSE)
  #'                         }
  betweenness <- function(graph,
                          v = V(graph),
                          directed = TRUE,
                          weights = NULL,
                          normalized = FALSE, ...) {
    igraph::betweenness(graph = graph,
                        v = v,
                        directed=directed,
                        weights = weights,
                        normalized = normalized, ...)
  }

  #*****************************************************************#

  #' @importFrom igraph degree
  #' @export
  igraph::degree

#=============================================================================
#
#    Code chunk 15: Funcstions for supporting additional file formats
#
#=============================================================================

  #' SIF to igraph
  #'
  #' This function imports and converts a SIF file from your local hard drive, cloud space,
  #' or internet into a graph with an igraph class, which can then be used for the identification
  #' of most influential nodes via the ivi function, for instance.
  #' @param Path A string or character vector indicating the path to the desired SIF file. The SIF file
  #' could be on your local hard drive, cloud space, or on the internet.
  #' @param directed Logical scalar, whether or not to create a directed graph.
  #' @return An igraph graph object.
  #' @keywords SIF.to.igraph
  #' @family network_reconstruction functions
  #' @export sif2igraph
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
  #' A shiny app has also been developed for Running the ExIR model, visualization of its results as well as computational
  #' simulation of knockout and/or up-regulation of its top candidate outputs, which is accessible using
  #' the `influential::runShinyApp("ExIR")` command.
  #' You can also access the shiny app online at https://influential.erc.monash.edu/.
  #' @param Desired_list (Optional) A character vector of your desired features. This vector could be, for
  #' instance, a list of features obtained from cluster analysis, time-course analysis,
  #' or a list of dysregulated features with a specific sign.
  #' @param Diff_data A dataframe of all significant differential/regression data and their
  #' statistical significance values (p-value/adjusted p-value). Note that the differential data
  #' should be in the log fold-change (log2FC) format.
  #' You may have selected a proportion of the differential data as the significant ones according
  #' to your desired thresholds. A function, named \code{\link[influential]{diff_data.assembly}}, has also been
  #' provided for the convenient assembling of the Diff_data dataframe.
  #' @param Diff_value An integer vector containing the column number(s) of the differential
  #' data in the Diff_data dataframe. The differential data could result from any type of
  #' differential data analysis. One example could be the fold changes (FCs) obtained from differential
  #' expression analyses. The user may provide as many differential data as he/she wish.
  #' @param Regr_value (Optional) An integer vector containing the column number(s) of the regression
  #' data in the Diff_data dataframe. The regression data could result from any type of regression
  #' data analysis or other analyses such as time-course data analyses that are based on regression models.
  #' @param Sig_value An integer vector containing the column number(s) of the significance values (p-value/adjusted p-value) of
  #' both differential and regression data (if provided). Providing significance values for the regression data is optional.
  #' @param Exptl_data A dataframe containing all of the experimental data including a column for specifying the conditions.
  #' The features/variables of the dataframe should be as the columns and the samples should come in the rows.
  #' The condition column should be of the character class. For example, if the study includes several replicates of
  #' cancer and normal samples, the condition column should include "cancer" and "normal" as the conditions of different samples.
  #' Also, the prior normalization of the experimental data is highly recommended. Otherwise,
  #' the user may set the Normalize argument to TRUE for a simple log2 transformation of the data.
  #' The experimental data could come from a variety sources such as transcriptomics and proteomics assays.
  #' @param Condition_colname A string or character vector specifying the name of the column "condition" of the Exptl_data dataframe.
  #' @param Normalize Logical; whether the experimental data should be normalized or not (default is FALSE). If TRUE, the
  #' experimental data will be log2 transformed.
  #' @param cor_thresh_method A character string indicating the method for filtering the correlation results, either
  #' "mr" (default; Mutual Rank) or "cor.coefficient".
  #' @param mr An integer determining the threshold of mutual rank for the selection of correlated features (default is 20). Note that
  #' higher mr values considerably increase the computation time.
  #' @param r The threshold of Spearman correlation coefficient for the selection of correlated features (default is 0.5).
  #' @param max.connections The maximum number of connections to be included in the association network.
  #' Higher max.connections might increase the computation time, cost, and accuracy of the results (default is 50,000).
  #' @param alpha The threshold of the statistical significance (p-value) used throughout the entire model (default is 0.05)
  #' @param num_trees Number of trees to be used for the random forests classification (supervised machine learning). Default is set to 10000.
  #' @param mtry Number of features to possibly split at in each node. Default is the (rounded down) square root of the
  #' number of variables. Alternatively, a single argument function returning an integer, given the number of independent variables.
  #' @param num_permutations Number of permutations to be used for computation of the statistical significance (p-values) of
  #' the importance scores resulted from random forests classification (default is 100).
  #' @param inf_const The constant value to be multiplied by the maximum absolute value of differential (logFC)
  #' values for the substitution with infinite differential values. This results in noticeably high biomarker values for features
  #' with infinite differential values compared with other features. Having said that, the user can still use the
  #' biomarker rank to compare all of the features. This parameter is ignored if no infinite value
  #' is present within Diff_data. However, this is used in the case of sc-seq experiments where some genes are uniquely
  #' expressed in a specific cell-type and consequently get infinite differential values. Note that the sign of differential
  #' value is preserved (default is 10^10).
  #' @param ncores Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of 
  #' cores to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").
  #' @param seed The seed to be used for all of the random processes throughout the model (default is 1234).
  #' @param verbose Logical; whether the accomplishment of different stages of the model should be printed (default is TRUE).
  #' @return A list of one graph and one to four tables including:
  #'
  #' - Driver table: Top candidate drivers
  #'
  #' - DE-mediator table: Top candidate differentially expressed/abundant mediators
  #'
  #' - nonDE-mediator table: Top candidate non-differentially expressed/abundant mediators
  #'
  #' - Biomarker table: Top candidate biomarkers
  #'
  #' The number of returned tables depends on the input data and specified arguments.
  #' @aliases ExIR
  #' @keywords exir
  #' @family integrative ranking functions
  #' @seealso \code{\link[influential]{exir.vis}},
  #' \code{\link[influential]{diff_data.assembly}},
  #' \code{\link[coop]{pcor}},
  #' \code{\link[stats]{prcomp}},
  #' \code{\link[ranger]{ranger}},
  #' \code{\link[ranger]{importance_pvalues}}
  #' @export exir
  #' @examples
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
                   cor_thresh_method = "mr", r = 0.5, mr = 20,
                   max.connections = 50000, alpha = 0.05,
                   num_trees = 10000, mtry = NULL, num_permutations = 100,
                   inf_const = 10^10, ncores = "default", seed = 1234, verbose = TRUE) {

    # Setup progress bar
    if(verbose) {
      pb <- utils::txtProgressBar(min = 0, max = 100, initial = 0, char = "=",
                           width = NA, style = 3, file = "")
    }

    #ProgressBar: Preparing the input data
    if(verbose) {
      print(unname(as.data.frame("Preparing the input data")),quote = FALSE, row.names = FALSE)
    }

    #make sure the input data is of data frame class
    Diff_data <- as.data.frame(Diff_data)
    Exptl_data <- as.data.frame(Exptl_data)
    Exptl_data[is.na(Exptl_data)] <- 0

    #change the colnames of Diff_data
    base::colnames(Diff_data) <- base::paste("source",
                                             base::colnames(Diff_data),
                                             sep = ".")

    #change the Inf/-Inf diff values (applicable to sc-Data)
    for(i in 1:base::length(Diff_value)) {

      if(any(base::is.infinite(Diff_data[,Diff_value[i]]))) {

        temp.max.abs.diff.value <-
      base::max(base::abs(Diff_data[,Diff_value[i]][!base::is.infinite(Diff_data[,Diff_value[i]])]))

        temp.inf.index <- base::which(base::is.infinite(Diff_data[,Diff_value[i]]))

        Diff_data[temp.inf.index, Diff_value[i]] <-
          base::ifelse(base::unlist(Diff_data[temp.inf.index,Diff_value[i]]) > 0,
                 temp.max.abs.diff.value*inf_const,
                 -1*temp.max.abs.diff.value*inf_const)
      }
    }

    # Get the column number of condition column
    condition.index <- match(Condition_colname, colnames(Exptl_data))
    
    # transform the data to numeric
    if(any(sapply(Exptl_data[,-condition.index], mode)  == "character")) {
      Exptl_data[,-condition.index] <- data.frame(sapply(Exptl_data[,-condition.index], as.numeric))
    }

    # Transform the condition column to a factor
    Exptl_data[,condition.index] <- base::as.factor(Exptl_data[,condition.index])

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
      print(unname(as.data.frame("Calculating the differential score")),quote = FALSE, row.names = FALSE)
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
      print(unname(as.data.frame("Calculating the regression/time-course R-squared score")),quote = FALSE, row.names = FALSE)
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
      print(unname(as.data.frame("Calculating the collective statistical significance of differential/regression factors")),quote = FALSE, row.names = FALSE)
    }

    #3 Calculate statistical significance of differential/regression factors
    if (max(Diff_data[,Sig_value]) > 1 | min(Diff_data[,Sig_value]) < 0) {
      stop("input Sig-values (p-value/padj) must all be in the range 0 to 1!")
    }

    for(m in 1:length(Sig_value)) {

      if(min(Diff_data[,Sig_value[m]])==0) {

        #range normalize the primitive Sig_value
        temp.min_Sig_value <- base::sort(base::unique(Diff_data[,Sig_value[m]]))[2]

        Diff_data[,Sig_value[m]] <- temp.min_Sig_value+
          (((Diff_data[,Sig_value[m]]-min(Diff_data[,Sig_value[m]]))*(max(Diff_data[,Sig_value[m]])-temp.min_Sig_value))/
             (max(Diff_data[,Sig_value[m]])-min(Diff_data[,Sig_value[m]])))
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

    #ProgressBar: Performing the random forests classification (supervised machine learning)
    if(verbose) {
      print(unname(as.data.frame("Performing the random forests classification (supervised machine learning)")),quote = FALSE, row.names = FALSE)
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

    #correct the names of features
      #first preserve a copy of original names
      features.exptl.for.super.learn <- colnames(exptl.for.super.learn)[-ncol(exptl.for.super.learn)]

    colnames(exptl.for.super.learn) <- janitor::make_clean_names(colnames(exptl.for.super.learn))
    
    # Defining the num.threads
    if(ncores == "default") {
      num.threads <- parallel::detectCores() - 1
    } else {
      num.threads <- ncores
    }

    #b Perform random forests classification
    base::set.seed(seed = seed)
    rf.diff.exptl <- ranger::ranger(formula = condition ~ .,
                                    data = exptl.for.super.learn,
                                    num.trees = num_trees,
                                    mtry = mtry,
                                    importance = "impurity_corrected",
                                    write.forest = FALSE, 
                                    num.threads = num.threads,
                                    seed = seed)

    base::set.seed(seed = seed)
    rf.diff.exptl.pvalue <- as.data.frame(ranger::importance_pvalues(x = rf.diff.exptl,
                                                                     formula = condition ~ .,
                                                                     num.permutations = num_permutations,
                                                                     data = exptl.for.super.learn,
                                                                     method = "altmann",
                                                                     num.threads = num.threads,
                                                                     seed = seed))

    #replace feature names (rownames) with their original names
    rownames(rf.diff.exptl.pvalue) <- features.exptl.for.super.learn

    if(any(is.na(rf.diff.exptl.pvalue[,"pvalue"])) |
       any(is.nan(rf.diff.exptl.pvalue[,"pvalue"]))) {
      rf.diff.exptl.pvalue[c(which(is.na(rf.diff.exptl.pvalue[,"pvalue"])),
                             which(is.nan(rf.diff.exptl.pvalue[,"pvalue"]))),
                           "pvalue"] <- 1
    }

    # filtering the RF output data
    select.number <- ifelse(!is.null(Desired_list),
                            round(length(Desired_list)/2),
                            100)

      if(length(which(rf.diff.exptl.pvalue[,"pvalue"] <alpha)) >= select.number) {
        rf.diff.exptl.pvalue <- base::subset(rf.diff.exptl.pvalue, rf.diff.exptl.pvalue$pvalue < alpha)

      } else {
        if(length(which(rf.diff.exptl.pvalue[, "pvalue"] < alpha)) >= 10) {
          rf.pval.select <- which(rf.diff.exptl.pvalue[, "pvalue"] < alpha)
        } else {
          rf.pval.select <- which(order(rf.diff.exptl.pvalue[, "pvalue"]) <= 10)
        }
        rf.nonSig <- seq(nrow(rf.diff.exptl.pvalue))[-rf.pval.select]
        required.pos.importance <- select.number - length(rf.pval.select)

        temp.rf.diff.exptl.pvalue <- rf.diff.exptl.pvalue[rf.nonSig,]

        rf.importance.select <- utils::tail(order(temp.rf.diff.exptl.pvalue[,"importance"]),
                                           n = required.pos.importance)

        temp.rf.diff.exptl.pvalue <- temp.rf.diff.exptl.pvalue[rf.importance.select,]

        #combine pvalue-based and importance-based tables
        rf.diff.exptl.pvalue <- rbind(rf.diff.exptl.pvalue[rf.pval.select,],
                                      temp.rf.diff.exptl.pvalue)
      }

    # negative importance values could be considered as 0
    if(any(rf.diff.exptl.pvalue[,"importance"] < 0)) {
      rf.diff.exptl.pvalue[which(rf.diff.exptl.pvalue[,"importance"] < 0),
                           "importance"] <- 0
    }

    # taking care of zero p-values
    if(min(rf.diff.exptl.pvalue[,"pvalue"])==0) {

      #range normalize the primitive pvalue
      temp.min_pvalue <- base::sort(base::unique(rf.diff.exptl.pvalue[,"pvalue"]))[2]

      rf.diff.exptl.pvalue[,"pvalue"] <- temp.min_pvalue+
        (((rf.diff.exptl.pvalue[,"pvalue"]-min(rf.diff.exptl.pvalue[,"pvalue"]))*(max(rf.diff.exptl.pvalue[,"pvalue"])-temp.min_pvalue))/
           (max(rf.diff.exptl.pvalue[,"pvalue"])-min(rf.diff.exptl.pvalue[,"pvalue"])))
    }

    #ProgressBar: Performing the random forests classification (supervised machine learning)
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 35)
    }


    #ProgressBar: Performing PCA (unsupervised machine learning)
    if(verbose) {
      print(unname(as.data.frame("Performing PCA (unsupervised machine learning)")),quote = FALSE, row.names = FALSE)
    }

    #5 Unsupervised machine learning (PCA)

    Exptl_data.for.PCA.index <- stats::na.omit(base::match(base::rownames(rf.diff.exptl.pvalue),
                                                           base::colnames(Exptl_data)))
    temp.Exptl_data.for.PCA <- Exptl_data[,Exptl_data.for.PCA.index]

    set.seed(seed)
    temp.PCA <- stats::prcomp(temp.Exptl_data.for.PCA)
    temp.PCA.r <- base::abs(temp.PCA$rotation[,1])

    #range normalize the rotation values
    temp.PCA.r <- 1+(((temp.PCA.r-min(temp.PCA.r))*(100-1))/
                       (max(temp.PCA.r)-min(temp.PCA.r)))

    #ProgressBar: Performing PCA (unsupervised machine learning)
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 40)
    }

    #ProgressBar: Performing the first round of association analysis
    if(verbose) {
      print(unname(as.data.frame("Performing the first round of association analysis")),quote = FALSE, row.names = FALSE)
    }

    #c Performing correlation analysis
    temp.corr <- fcor(data = Exptl_data[,-condition.index],
                      method = "spearman", mutualRank = ifelse(cor_thresh_method == "mr", TRUE, FALSE))

    #save a second copy of all cor data
    temp.corr.for.sec.round <- temp.corr

    #filter corr data for only those corr between diff features and themselves/others
    filter.corr.index <- stats::na.omit(base::unique(c(base::which(temp.corr$row %in% rownames(rf.diff.exptl.pvalue)),
                                                       base::which(temp.corr$column %in% rownames(rf.diff.exptl.pvalue)))))
    temp.corr <- temp.corr[filter.corr.index,]

    #filtering low level correlations
    cor.thresh <- r
    mr.thresh <- mr

    if(cor_thresh_method == "mr") {

      temp.corr <- base::subset(temp.corr, temp.corr[,4] < mr.thresh)

      if(nrow(temp.corr)> (max.connections*0.95)) {

        temp.corr.select.index <- utils::head(order(temp.corr$mr),
                                              n = round(max.connections*0.95))

        temp.corr <- temp.corr[temp.corr.select.index,]

      }

    } else if(cor_thresh_method == "cor.coefficient") {

      temp.corr <- base::subset(temp.corr, base::abs(temp.corr[,3])>cor.thresh)

      if(nrow(temp.corr)> (max.connections*0.95)) {

        temp.corr.select.index <- utils::tail(order(temp.corr$cor),
                                              n = round(max.connections*0.95))

        temp.corr <- temp.corr[temp.corr.select.index,]

      }
    }

    diff.only.temp.corr <- temp.corr

    #getting the list of diff features and their correlated features
    diff.plus.corr.features <- base::unique(c(base::as.character(temp.corr[,1]),
                                              base::as.character(temp.corr[,2])))

    #find the diff features amongst diff.plus.corr.features
    diff.only.features.index <- stats::na.omit(base::unique(base::match(rownames(rf.diff.exptl.pvalue),
                                                                        diff.plus.corr.features)))
    non.diff.only.features <- diff.plus.corr.features[-diff.only.features.index]

    #ProgressBar: Performing first round of association analysis
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 50)
    }

    #ProgressBar: Performing the second round of association analysis
    if(verbose) {
      print(unname(as.data.frame("Performing the second round of association analysis")),quote = FALSE, row.names = FALSE)
    }

    if(base::length(non.diff.only.features) > 0) {

    #redo correlation analysis
    temp.corr <- temp.corr.for.sec.round
    rm(temp.corr.for.sec.round)

    #filter corr data for only those corr between non.diff.only.features and themselves/others
    filter.corr.index <- stats::na.omit(base::unique(c(base::which(temp.corr$row %in% non.diff.only.features),
                                                       base::which(temp.corr$column %in% non.diff.only.features))))
    temp.corr <- temp.corr[filter.corr.index,]

    #filtering low level correlations
    cor.thresh <- r
    mr.thresh <- mr

    if(cor_thresh_method == "mr") {
      temp.corr <- base::subset(temp.corr, temp.corr[,4] < mr.thresh)

    } else if(cor_thresh_method == "cor.coefficient") {
      temp.corr <- base::subset(temp.corr, base::abs(temp.corr[,3])>cor.thresh)
    }

    # separate non.diff.only features
    temp.corr.diff.only.index <- stats::na.omit(base::unique(c(base::which(temp.corr$row %in% rownames(rf.diff.exptl.pvalue)),
                                                               base::which(temp.corr$column %in% rownames(rf.diff.exptl.pvalue)))))

    if(base::length(temp.corr.diff.only.index)>0) {
      temp.corr <- temp.corr[-temp.corr.diff.only.index,]
    }

    if(nrow(temp.corr)>(max.connections-nrow(diff.only.temp.corr))) {

      if(cor_thresh_method == "mr") {
        temp.corr.select.index <- utils::head(order(temp.corr$mr),
                                              n = (max.connections-nrow(diff.only.temp.corr)))

      } else if(cor_thresh_method == "cor.coefficient") {
        temp.corr.select.index <- utils::tail(order(temp.corr$cor),
                                              n = (max.connections-nrow(diff.only.temp.corr)))
      }

      temp.corr <- temp.corr[temp.corr.select.index,]
    }

    # recombine the diff.only.temp.corr data and temp.corr
    temp.corr <- base::rbind(temp.corr, diff.only.temp.corr)

    } else {
      temp.corr <- diff.only.temp.corr
      rm(temp.corr.for.sec.round, diff.only.temp.corr)
    }

    #ProgressBar: Performing second round of association analysis
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 60)
    }

    #ProgressBar: Network reconstruction
    if(verbose) {
      print(unname(as.data.frame("Network reconstruction")),quote = FALSE, row.names = FALSE)
    }

    #d Graph reconstruction
    temp.corr.graph <- igraph::graph_from_data_frame(temp.corr[,c(1:2)])

    #ProgressBar: Network reconstruction
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 65)
    }

    #ProgressBar: Calculation of the integrated value of influence (IVI)
    if(verbose) {
      print(unname(as.data.frame("Calculation of the integrated value of influence (IVI)")),quote = FALSE, row.names = FALSE)
    }

    #e Calculation of IVI
    temp.corr.ivi <- ivi(temp.corr.graph, ncores = ncores)

    #ProgressBar: Calculation of the integrated value of influence (IVI)
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 70)
    }

    #ProgressBar: Calculation of the primitive driver score
    if(verbose) {
      print(unname(as.data.frame("Calculation of the primitive driver score")),quote = FALSE, row.names = FALSE)
    }

    ## Driver score and ranking

    #a calculate first level driver score based on #3 and #4

    Diff_data.IVI.index <- stats::na.omit(match(names(temp.corr.ivi),
                                                rownames(Diff_data)))

    if(length(Diff_data.IVI.index) > 0) {

      Diff_data$IVI <- 0

      temp.corr.ivi.for.Diff_data.IVI.index <- stats::na.omit(match(rownames(Diff_data)[Diff_data.IVI.index],
                                                                    names(temp.corr.ivi)))

      Diff_data$IVI[Diff_data.IVI.index] <- temp.corr.ivi[temp.corr.ivi.for.Diff_data.IVI.index]

      #range normalize the IVI
      Diff_data$IVI <- 1+(((Diff_data$IVI-min(Diff_data$IVI))*(100-1))/
                            (max(Diff_data$IVI)-min(Diff_data$IVI)))

    } else {
      Diff_data$IVI <- 1
    }

    Diff_data$first.Driver.Rank <- 1

    for (i in 1:nrow(Diff_data)) {
      if(c(any(Diff_data[i,Diff_value, drop = FALSE]<0) & any(Diff_data[i,Diff_value, drop = FALSE]>0))) {
        Diff_data$first.Driver.Rank[i] <- 0
      } else {
        Diff_data$first.Driver.Rank[i] <- Diff_data$sum.Sig_value[i]*Diff_data$IVI[i]
      }
    }
    #range normalize the first driver rank
    if(any(Diff_data$first.Driver.Rank == 0)) {
      Diff_data$first.Driver.Rank <- 0+(((Diff_data$first.Driver.Rank-min(Diff_data$first.Driver.Rank))*(100-0))/
                                          (max(Diff_data$first.Driver.Rank)-min(Diff_data$first.Driver.Rank)))
    } else {
      Diff_data$first.Driver.Rank <- 1+(((Diff_data$first.Driver.Rank-min(Diff_data$first.Driver.Rank))*(100-1))/
                                          (max(Diff_data$first.Driver.Rank)-min(Diff_data$first.Driver.Rank)))
    }

    #ProgressBar: Calculation of the primitive driver score
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 75)
    }

    #ProgressBar: Calculation of the neighborhood driver score
    if(verbose) {
      print(unname(as.data.frame("Calculation of the neighborhood driver score")),quote = FALSE, row.names = FALSE)
    }

    #b (#6) calculate neighborhood score

    #get the list of network nodes
    network.nodes <- base::as.character(igraph::as_ids(V(temp.corr.graph)))

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
      print(unname(as.data.frame("Preparation of the driver table")),quote = FALSE, row.names = FALSE)
    }

    Diff_data$N.score <- 0
    Diff_data.N.score.index <- stats::na.omit(match(neighborehood.score.table$node,
                                                    rownames(Diff_data)))

    neighborehood.score.table.for.Diff_data.N.score.index <- stats::na.omit(match(rownames(Diff_data)[Diff_data.N.score.index],
                                                                                  neighborehood.score.table$node))

    Diff_data$N.score[Diff_data.N.score.index] <- neighborehood.score.table$N.score[neighborehood.score.table.for.Diff_data.N.score.index]

    #range normalize (1,100) the neighborhood score
    Diff_data$N.score <- ifelse(sum(Diff_data$N.score) == 0, 1,
                                1+(((Diff_data$N.score-min(Diff_data$N.score))*(100-1))/
                                     (max(Diff_data$N.score)-min(Diff_data$N.score))))

    #c calculate the final driver score

    Diff_data$final.Driver.score <- (Diff_data$first.Driver.Rank)*(Diff_data$N.score)
    Diff_data$final.Driver.score[Diff_data$final.Driver.score==0] <- NA

    # Create the Drivers table

    Driver.table <- Diff_data

    #remove the rows/features with NA in the final driver score
    Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

    #filter the driver table by the desired list (if provided)
    if(!is.null(Desired_list)) {
      Driver.table.row.index <- stats::na.omit(match(Desired_list,
                                                     rownames(Driver.table)))
      Driver.table <- Driver.table[Driver.table.row.index,]
    }

    if(nrow(as.data.frame(Driver.table))==0) {Driver.table <- NULL} else {
      
      #add Z.score
      Driver.table$Z.score <- base::scale(Driver.table$final.Driver.score)

      #range normalize final driver score
      ifelse(length(unique(Driver.table$final.Driver.score)) > 1,
             Driver.table$final.Driver.score <- 1+(((Driver.table$final.Driver.score-min(Driver.table$final.Driver.score))*(100-1))/
                                                     (max(Driver.table$final.Driver.score)-min(Driver.table$final.Driver.score))),
             Driver.table$final.Driver.score <- 1)

      #add driver rank
      Driver.table$rank <- rank(-Driver.table$final.Driver.score,
                                ties.method = "min")

      #add P-value
      Driver.table$p.value <- stats::pnorm(Driver.table$Z.score,
                                           lower.tail = FALSE)

      #add adjusted pvalue
      Driver.table$padj <- stats::p.adjust(p = Driver.table$p.value,
                                           method = "BH")

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
                                      "Z.score",
                                      "rank",
                                      "p.value",
                                      "padj",
                                      "driver.type")]

      #rename column names
      colnames(Driver.table) <- c("Score", "Z.score",
                                  "Rank", "P.value",
                                  "P.adj", "Type")

      #filtering redundant (NaN) results
      Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

      if(nrow(as.data.frame(Driver.table))==0) {Driver.table <- NULL}

    }

    #ProgressBar: Preparation of the driver table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 85)
    }


    #ProgressBar: Preparation of the biomarker table
    if(verbose) {
      print(unname(as.data.frame("Preparation of the biomarker table")),quote = FALSE, row.names = FALSE)
    }

    # Create the Biomarker table

    Biomarker.table <- Diff_data

    #remove the rows/features with NA in the final driver score
    Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]

    #filter the biomarker table by the desired list (if provided)
    if(!is.null(Desired_list)) {
      Biomarker.table.row.index <- stats::na.omit(match(Desired_list,
                                                     rownames(Biomarker.table)))
      Biomarker.table <- Biomarker.table[Biomarker.table.row.index,]
    }

    if(nrow(as.data.frame(Biomarker.table))==0) {Biomarker.table <- NULL} else {

      #add RF importance score and p-value
      Biomarker.table$rf.importance <- 0
      Biomarker.table$rf.pvalue <- 1

      Biomarker.table.rf.index <- stats::na.omit(match(rownames(rf.diff.exptl.pvalue),
                                                       rownames(Biomarker.table)))

      rf.for.Biomarker.table <- stats::na.omit(match(rownames(Biomarker.table)[Biomarker.table.rf.index],
                                                     rownames(rf.diff.exptl.pvalue)))

      Biomarker.table$rf.importance[Biomarker.table.rf.index] <- rf.diff.exptl.pvalue$importance[rf.for.Biomarker.table]
      Biomarker.table$rf.pvalue[Biomarker.table.rf.index] <- rf.diff.exptl.pvalue$pvalue[rf.for.Biomarker.table]

      #range normalize rf.importance and rf.pvalue
      Biomarker.table$rf.importance <- 1+(((Biomarker.table$rf.importance-min(Biomarker.table$rf.importance))*(100-1))/
                                            (max(Biomarker.table$rf.importance)-min(Biomarker.table$rf.importance)))

      Biomarker.table$rf.pvalue <- -log10(Biomarker.table$rf.pvalue)
      Biomarker.table$rf.pvalue <- 1+(((Biomarker.table$rf.pvalue-min(Biomarker.table$rf.pvalue))*(100-1))/
                                        (max(Biomarker.table$rf.pvalue)-min(Biomarker.table$rf.pvalue)))

      #add rotation values
      Biomarker.table$rotation <- 0
      Biomarker.table.for.rotation <- stats::na.omit(match(names(temp.PCA.r),
                                                           rownames(Biomarker.table)))

      Biomarker.table.rotation.index <- stats::na.omit(match(rownames(Biomarker.table)[Biomarker.table.for.rotation],
                                                             names(temp.PCA.r)))

      Biomarker.table$rotation[Biomarker.table.for.rotation] <- temp.PCA.r[Biomarker.table.rotation.index]

      #range normalize rotation values
      Biomarker.table$rotation <- 1+(((Biomarker.table$rotation-min(Biomarker.table$rotation))*(100-1))/
                                            (max(Biomarker.table$rotation)-min(Biomarker.table$rotation)))

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
      
      #add biomarker Z.score
      Biomarker.table$Z.score <- base::scale(Biomarker.table$final.biomarker.score)

      #range normalize biomarker score
      ifelse(length(unique(Biomarker.table$final.biomarker.score)) > 1,
             Biomarker.table$final.biomarker.score <- 1+(((Biomarker.table$final.biomarker.score-min(Biomarker.table$final.biomarker.score))*(100-1))/
                                                           (max(Biomarker.table$final.biomarker.score)-min(Biomarker.table$final.biomarker.score))),
             Biomarker.table$final.biomarker.score <- 1)

      #add biomarker rank
      Biomarker.table$rank <- rank(-Biomarker.table$final.biomarker.score, ties.method = "min")

      #add biomarker P-value
      Biomarker.table$P.value <- stats::pnorm(Biomarker.table$Z.score,
                                              lower.tail = FALSE)

      #add biomarker adjusted p-value
      Biomarker.table$padj <- stats::p.adjust(p = Biomarker.table$P.value,
                                              method = "BH")

      #add biomarker type
      Biomarker.table$type <- ""

      for (b in 1:nrow(Biomarker.table)) {

        if(sum(Biomarker.table[b,Diff_value])<0) {
          Biomarker.table$type[b] <- "Down-regulated"

        } else if(sum(Biomarker.table[b,Diff_value])>0) {
          Biomarker.table$type[b] <- "Up-regulated"
        } else {
          Biomarker.table$type[b] <- NA
        }
      }

      #remove redundent columns
      Biomarker.table <- Biomarker.table[,c("final.biomarker.score",
                                            "Z.score", "rank",
                                            "P.value", "padj", "type")]

      #rename column names
      colnames(Biomarker.table) <- c("Score", "Z.score",
                                     "Rank", "P.value",
                                     "P.adj", "Type")

      #filtering redundant (NaN) results
      Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]

      if(nrow(as.data.frame(Biomarker.table))==0) {Biomarker.table <- NULL}

    }

    #ProgressBar: Preparation of the biomarker table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 90)
    }

    #ProgressBar: Preparation of the DE-mediator table
    if(verbose) {
      print(unname(as.data.frame("Preparation of the DE-mediator table")),quote = FALSE, row.names = FALSE)
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

      DE.mediator.row.index <- stats::na.omit(match(network.nodes,
                                                    rownames(DE.mediator.table)))

      if(!is.null(Desired_list)) {
        desired.DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                           rownames(DE.mediator.table)[DE.mediator.row.index]))
        DE.mediator.row.index <- DE.mediator.row.index[desired.DE.mediator.row.index]
      }

    DE.mediator.table <- DE.mediator.table[DE.mediator.row.index,]
    if(nrow(as.data.frame(DE.mediator.table))==0) {DE.mediator.table <- NULL} else {
      
      #add DE mediators Z score
      DE.mediator.table$Z.score <- base::scale(DE.mediator.table$DE.mediator.score)

      #range normalize DE mediators score
      ifelse(length(unique(DE.mediator.table$DE.mediator.score)) > 1,
             DE.mediator.table$DE.mediator.score <- 1+(((DE.mediator.table$DE.mediator.score-min(DE.mediator.table$DE.mediator.score))*(100-1))/
                                                         (max(DE.mediator.table$DE.mediator.score)-min(DE.mediator.table$DE.mediator.score))),
             DE.mediator.table$DE.mediator.score <- 1)

      #add DE mediators rank
      DE.mediator.table$rank <- rank(-DE.mediator.table$DE.mediator.score, ties.method = "min")

      #add DE mediators P-value
      DE.mediator.table$P.value <- stats::pnorm(DE.mediator.table$Z.score,
                                                lower.tail = FALSE)

      #add DE mediators adjusted P-value
      DE.mediator.table$padj <- stats::p.adjust(p = DE.mediator.table$P.value,
                                                method = "BH")

      #remove redundent columns
      DE.mediator.table <- DE.mediator.table[,c("DE.mediator.score", "Z.score",
                                                "rank", "P.value", "padj")]

      #rename column names
      colnames(DE.mediator.table) <- c("Score", "Z.score",
                                       "Rank", "P.value",
                                       "P.adj")

      #filtering redundant (NaN) results
      DE.mediator.table <- DE.mediator.table[stats::complete.cases(DE.mediator.table),]

      if(nrow(as.data.frame(DE.mediator.table))==0) {
        DE.mediator.table <- NULL
        } else {
          # Adding the associated drivers to the table
          
          ## First order drivers
          first_order_drivers <- 
            lapply(1:nrow(DE.mediator.table), function(i) {
              
              first_order_assoc_drivers <-
                Driver.table[grep(paste0(paste("^",
                                               igraph::as_ids(igraph::neighborhood(temp.corr.graph, 
                                                                                   nodes = rownames(DE.mediator.table)[i], 
                                                                                   order = 1)[[1]]),
                                               "$", sep = ""), 
                                         collapse = "|"), 
                                  rownames(Driver.table), value = TRUE),]
              
              first_order_assoc_drivers <- rownames(first_order_assoc_drivers[order(first_order_assoc_drivers$Rank),])
              first_order_assoc_drivers
            })
          
          ## Second order drivers
          second_order_drivers <- 
            lapply(1:nrow(DE.mediator.table), function(i) {
              
              second_order_assoc_drivers <-
                Driver.table[grep(paste0(paste("^",
                                               igraph::as_ids(igraph::neighborhood(temp.corr.graph, 
                                                                                   nodes = rownames(DE.mediator.table)[i], 
                                                                                   order = 2)[[1]]),
                                               "$", sep = ""), 
                                         collapse = "|"), 
                                  rownames(Driver.table), value = TRUE),]
              
              second_order_assoc_drivers <- second_order_assoc_drivers[!(rownames(second_order_assoc_drivers) %in% first_order_drivers[[i]]),]
              
              second_order_assoc_drivers <- rownames(second_order_assoc_drivers[order(second_order_assoc_drivers$Rank),])
              second_order_assoc_drivers
            })
          
          ## Collapsing the associated drivers
          first_order_drivers <- lapply(first_order_drivers, function(i) {
            paste0(i, collapse = ", ")
          })
          
          second_order_drivers <- lapply(second_order_drivers, function(i) {
            paste0(i, collapse = ", ")
          })
          
          ## Adding to the table
          DE.mediator.table$First.order.Drivers <- unlist(first_order_drivers)
          DE.mediator.table$Second.order.Drivers <- unlist(second_order_drivers)
        }
    }

    #ProgressBar: Preparation of the DE-mediator table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 95)
    }

    #ProgressBar: Preparation of the nonDE-mediator table
    if(verbose) {
      print(unname(as.data.frame("Preparation of the nonDE-mediator table")),quote = FALSE, row.names = FALSE)
    }

    # Create the non-DE mediators table
    non.DE.mediators.index <- stats::na.omit(unique(match(rownames(Diff_data),
                                                          neighborehood.score.table$node)))

    non.DE.mediators.table <- neighborehood.score.table[-c(non.DE.mediators.index),]
    if(nrow(as.data.frame(non.DE.mediators.table))==0) {non.DE.mediators.table <- NULL} else {

      #filter the non-DE mediators table by either the desired list
      if(!is.null(Desired_list)) {
        non.DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                          non.DE.mediators.table$node))
        non.DE.mediators.table <- non.DE.mediators.table[non.DE.mediator.row.index,]
      }
    }

    if(nrow(as.data.frame(non.DE.mediators.table))==0) {non.DE.mediators.table <- NULL} else {

      rownames(non.DE.mediators.table) <- non.DE.mediators.table$node

      non.DE.mediators.ivi.index <- stats::na.omit(match(rownames(non.DE.mediators.table),
                                                         names(temp.corr.ivi)))

      non.DE.mediators.table$ivi <- temp.corr.ivi[non.DE.mediators.ivi.index]

      non.DE.mediators.table$non.DE.mediator.score <- non.DE.mediators.table$N.score*non.DE.mediators.table$ivi
      
      #add non-DE mediators Z.score
      non.DE.mediators.table$Z.score <- base::scale(non.DE.mediators.table$non.DE.mediator.score)

      #range normalize nonDE mediators score
      ifelse(length(unique(non.DE.mediators.table$non.DE.mediator.score)) > 1,
             non.DE.mediators.table$non.DE.mediator.score <- 1+(((non.DE.mediators.table$non.DE.mediator.score-min(non.DE.mediators.table$non.DE.mediator.score))*(100-1))/
                                                                  (max(non.DE.mediators.table$non.DE.mediator.score)-min(non.DE.mediators.table$non.DE.mediator.score))),
             non.DE.mediators.table$non.DE.mediator.score <- 1)

      #add non-DE mediators P-value
      non.DE.mediators.table$P.value <- stats::pnorm(non.DE.mediators.table$Z.score,
                                                     lower.tail = FALSE)

      #add non-DE mediators adjusted p-value
      non.DE.mediators.table$padj <- stats::p.adjust(p = non.DE.mediators.table$P.value,
                                                     method = "BH")

      #add non-DE mediators rank
      non.DE.mediators.table$rank <- rank(-non.DE.mediators.table$non.DE.mediator.score, ties.method = "min")

      #remove redundent columns
      non.DE.mediators.table <- non.DE.mediators.table[,c("non.DE.mediator.score",
                                                          "Z.score", "rank",
                                                          "P.value", "padj")]

      #rename column names
      colnames(non.DE.mediators.table) <- c("Score", "Z.score",
                                            "Rank", "P.value",
                                            "P.adj")

      #filtering redundant (NaN) results
      non.DE.mediators.table <- non.DE.mediators.table[stats::complete.cases(non.DE.mediators.table),]

      if(nrow(as.data.frame(non.DE.mediators.table))==0) {
        non.DE.mediators.table <- NULL
      } else {
        # Adding the associated drivers to the table
        
        ## First order drivers
        first_order_drivers <- 
          lapply(1:nrow(non.DE.mediators.table), function(i) {
            
            first_order_assoc_drivers <-
              Driver.table[grep(paste0(paste("^",
                                             igraph::as_ids(igraph::neighborhood(temp.corr.graph, 
                                                                                 nodes = rownames(non.DE.mediators.table)[i], 
                                                                                 order = 1)[[1]]),
                                             "$", sep = ""), 
                                       collapse = "|"), 
                                rownames(Driver.table), value = TRUE),]
            
            first_order_assoc_drivers <- rownames(first_order_assoc_drivers[order(first_order_assoc_drivers$Rank),])
            first_order_assoc_drivers
          })
        
        ## Second order drivers
        second_order_drivers <- 
          lapply(1:nrow(non.DE.mediators.table), function(i) {
            
            second_order_assoc_drivers <-
              Driver.table[grep(paste0(paste("^",
                                             igraph::as_ids(igraph::neighborhood(temp.corr.graph, 
                                                                                 nodes = rownames(non.DE.mediators.table)[i], 
                                                                                 order = 2)[[1]]),
                                             "$", sep = ""), 
                                       collapse = "|"), 
                                rownames(Driver.table), value = TRUE),]
            
            second_order_assoc_drivers <- second_order_assoc_drivers[!(rownames(second_order_assoc_drivers) %in% first_order_drivers[[i]]),]
            
            second_order_assoc_drivers <- rownames(second_order_assoc_drivers[order(second_order_assoc_drivers$Rank),])
            second_order_assoc_drivers
          })
        
        ## Collapsing the associated drivers
        first_order_drivers <- lapply(first_order_drivers, function(i) {
          paste0(i, collapse = ", ")
        })
        
        second_order_drivers <- lapply(second_order_drivers, function(i) {
          paste0(i, collapse = ", ")
        })
        
        ## Adding to the table
        non.DE.mediators.table$First.order.Drivers <- unlist(first_order_drivers)
        non.DE.mediators.table$Second.order.Drivers <- unlist(second_order_drivers)
        }
    }

    #ProgressBar: Preparation of the nonDE-mediator table
    if(verbose) {
      utils::setTxtProgressBar(pb = pb, value = 100)
    }


    Results <- list("Driver table" = Driver.table,
                    "DE-mediator table" = DE.mediator.table,
                    "nonDE-mediator table" = non.DE.mediators.table,
                    "Biomarker table" = Biomarker.table,
                    "Graph" = temp.corr.graph)

    Results.Non.Null.vector <- vector()

    for (i in 1:length(Results)) {
      Results.Non.Null.vector[i] <- !is.null(Results[[i]])
    }

    Results <- Results[Results.Non.Null.vector]

    # set the class of Results
    base::class(Results) <- "ExIR_Result"

    return(Results)
  }

  #=============================================================================
  #
  #    Code chunk 17: Assembling the differential/regression data in a dataframe
  #
  #=============================================================================

  #' Assembling the differential/regression data
  #'
  #' This function assembles a dataframe required for running the \emph{\strong{\code{ExIR}}} model. You may provide
  #' as many differential/regression data as you wish. Also, the datasets should be filtered
  #' beforehand according to your desired thresholds and, consequently, should only include the significant data.
  #' Each dataset provided should be a dataframe with one or two columns.
  #' The first column should always include differential/regression values
  #' and the second one (if provided) the significance values. Please also note that the significance (adjusted P-value)
  #' column is mandatory for differential datasets.
  #' @param ... Desired datasets/dataframes.
  #' @return A dataframe including the collective list of features in rows and all of the
  #' differential/regression data and their statistical significance in columns with the same
  #' order provided by the user.
  #' @aliases DDA
  #' @keywords diff_data.assembly
  #' @seealso \code{\link[influential]{exir}}
  #' @export diff_data.assembly
  #' @examples
  #' \dontrun{
  #' my.Diff_data <- diff_data.assembly(Differential_data1,
  #'                                    Differential_data2,
  #'                                    Regression_data1)
  #' }
  diff_data.assembly <- function(...) {

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

  #=============================================================================
  #
  #    Code chunk 18: Visualization of a graph based on centrality measures
  #
  #=============================================================================

  #' Centrality-based network visualization
  #'
  #' This function has been developed for the visualization of a network based on
  #' applying a centrality measure to the size and color of network nodes. You are
  #' also able to adjust the directedness and weight of connections. Some of the documentations
  #' of the arguments of this function have been adapted from ggplot2 and igraph packages.
  #' A shiny app has also been developed for the calculation of IVI as well as IVI-based network
  #' visualization, which is accessible using the `influential::runShinyApp("IVI")` command.
  #' You can also access the shiny app online at https://influential.erc.monash.edu/.
  #' @param graph A graph (network) of the igraph class.
  #' @param cent.metric A numeric vector of the desired centrality measure previously
  #' calculated by any means. For example, you may use the function \code{\link[influential]{ivi}}
  #' for the calculation of the Integrated Value of Influence (IVI) of network nodes. Please note that
  #' if the centrality measure has been calculated by any means other than the \code{influential} package, make
  #' sure that the order of the values in the \code{cent.metric} vector is consistent with the order of vertices
  #' in the network \code{(V(graph))}.
  #' @param layout The layout to be used for organizing network nodes. Current available layouts include
  #' \code{"kk", "star", "tree", "components", "circle", "automatic", "grid",
  #' "sphere", "random", "dh", "drl", "fr", "gem", "graphopt", "lgl", "mds", and "sugiyama"}
  #' (default is set to "kk"). For a complete description of different layouts and their
  #' underlying algorithms please refer to the function \code{\link[igraph]{layout_}}.
  #' @param node.group A vector of the same length as the number of network nodes defining the group each node of the network belongs to.
  #' @param node.color A character string indicating the colormap option to use.
  #' Five options are available: "magma" (or "A"), "inferno" (or "B"), "plasma"
  #' (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E").
  #' @param node.size.min The size of nodes with the lowest value of the centrality measure (default is set to 3).
  #' @param node.size.max The size of nodes with the highest value of the centrality measure (default is set to 15).
  #' @param dist.power The power to be used to visualize more distinction between nodes with high and low
  #' centrality measure values. The higher the power, the smaller the nodes with lower values of the centrality
  #' measure will become. Default is set to 1, meaning the relative sizes of nodes are reflective of their
  #' actual centrality measure values.
  #' @param node.shape The shape of nodes. Current available shapes include \code{"circle",
  #' "square", "diamond", "triangle", and "inverted triangle"} (default is set to "circle"). You can also
  #' set different shapes to different groups of nodes by providing a character vector of shapes of nodes with
  #' the same length and order of network vertices. This is useful when plotting a network that include different
  #' type of node (for example, up- and down-regulated features).
  #' @param stroke.size The size of stroke (border) around the nodes (default is set to 1.5).
  #' @param stroke.color The color of stroke (border) around the nodes (default is set to "identical" meaning that the
  #' stroke color of a node will be identical to its corresponding node color). You can also
  #' set different colors to different groups of nodes by providing a character vector of colors of nodes with
  #' the same length and order of network vertices. This is useful when plotting a network that include different
  #' type of node (for example, up- and down-regulated features).
  #' @param stroke.alpha The transparency of the stroke (border) around the nodes which should
  #' be a number between 0 and 1 (default is set to 0.6).
  #' @param show.labels Logical scalar, whether to show node labels or not (default is set to TRUE).
  #' @param label.cex The amount by which node labels should be scaled relative to the node sizes (default is set to 0.4).
  #' @param label.color The color of node labels (default is set to "black").
  #' @param directed Logical scalar, whether to draw the network as directed or not (default is set to FALSE).
  #' @param arrow.width The width of arrows in the case the network is directed (default is set to 25).
  #' @param arrow.length The length of arrows in inch in the case the network is directed (default is set to 0.07).
  #' @param edge.width The constant width of edges if the network is unweighted (default is set to 0.5).
  #' @param weighted Logical scalar, whether the network is a weighted network or not (default is set to FALSE).
  #' @param edge.width.min The width of edges with the lowest weight (default is set to 0.2).
  #' This parameter is ignored for unweighted networks.
  #' @param edge.width.max The width of edges with the highest weight (default is set to 1).
  #' This parameter is ignored for unweighted networks.
  #' @param edge.color The color of edges (default is set to "grey75").
  #' @param edge.linetype The line type of edges. Current available linetypes include
  #' \code{"twodash", "longdash", "dotdash", "dotted", "dashed", and "solid"} (default is set to "solid").
  #' @param legend.position The position of legends ("none", "left", "right",
  #' "bottom", "top", or two-element numeric vector). The default is set to "right".
  #' @param legend.direction layout of items in legends ("horizontal" or "vertical").
  #' The default is set to "vertical".
  #' @param legend.title The legend title in the string format (default is set to "Centrality measure").
  #' @param boxed.legend Logical scalar, whether to draw a box around the legend or not (default is set to TRUE).
  #' @param show.plot.title Logical scalar, whether to show the plot title or not (default is set to TRUE).
  #' @param plot.title The plot title in the string format (default is set to "Centrality Measure-based Network").
  #' @param title.position The position of title ("left", "center", or "right"). The default is set to "center".
  #' @param show.bottom.border Logical scalar, whether to draw the bottom border line (default is set to TRUE).
  #' @param show.left.border Logical scalar, whether to draw the left border line (default is set to TRUE).
  #' @param seed A single value, interpreted as an integer to be used for random number generation for preparing
  #' the network layout (default is set to 1234).
  #' @return A plot with the class ggplot.
  #' @keywords cent_network.vis
  #' @family visualization functions
  #' @seealso \code{\link[influential]{ivi}}
  #' @export cent_network.vis
  #' @examples
  #' \dontrun{
  #' MyData <- coexpression.data
  #' My_graph <- graph_from_data_frame(MyData)
  #' Graph_IVI <- ivi(graph = My_graph, mode = "all")
  #' Graph_IVI_plot <- cent_network.vis(graph = My_graph, cent.metric = Graph_IVI,
  #'                                    legend.title = "IVI",
  #'                                    plot.title = "IVI-based Network")
  #' }
  cent_network.vis <- function(graph,
                               cent.metric,
                               layout = "kk",
                               node.group = NULL,
                               node.color = "viridis",
                               node.size.min = 3,
                               node.size.max = 15,
                               dist.power = 1,
                               node.shape = "circle",
                               stroke.size = 1.5,
                               stroke.color = "identical",
                               stroke.alpha = 0.6,
                               show.labels = TRUE,
                               label.cex = 0.4,
                               label.color = "black",
                               directed = FALSE,
                               arrow.width = 25,
                               arrow.length = 0.07,
                               edge.width = 0.5,
                               weighted = FALSE,
                               edge.width.min = 0.2,
                               edge.width.max = 1,
                               edge.color = "grey75",
                               edge.linetype = "solid",
                               legend.position = "right",
                               legend.direction = "vertical",
                               legend.title = "Centrality\nmeasure",
                               boxed.legend = TRUE,
                               show.plot.title = TRUE,
                               plot.title = "Centrality Measure-based Network",
                               title.position = "center",
                               show.bottom.border = TRUE,
                               show.left.border = TRUE,
                               seed = 1234) {

  # preparing the layout
  if(layout == "kk") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_kk(graph))
  } else if(layout == "star") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_as_star(graph))
  } else if(layout == "tree") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_as_tree(graph))
  } else if(layout == "components") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_components(graph))
  } else if(layout == "circle") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_in_circle(graph))
  } else if(layout == "automatic") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_nicely(graph))
  } else if(layout == "grid") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_on_grid(graph))
  } else if(layout == "dh") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_dh(graph))
  } else if(layout == "sphere") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_on_sphere(graph))
  } else if(layout == "random") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_randomly(graph))
  } else if(layout == "drl") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_drl(graph))
  } else if(layout == "fr") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_fr(graph))
  } else if(layout == "gem") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_gem(graph))
  } else if(layout == "graphopt") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_graphopt(graph))
  } else if(layout == "lgl") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_lgl(graph))
  } else if(layout == "mds") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_mds(graph))
  } else if(layout == "sugiyama") {
    base::set.seed(seed = seed)
    plotcord <- base::data.frame(igraph::layout_with_sugiyama(graph))
  }

  ####*******************************####

  # preparing the plotcord table
  base::rownames(plotcord) <- igraph::as_ids(V(graph))
  base::colnames(plotcord) = c("X","Y")

  # add the centrality measure
  plotcord$cent.metric <- cent.metric

  # range normalize the Node size based on the centrality measure
  plotcord$Node.size <- node.size.min+((((cent.metric)^dist.power-min((cent.metric)^dist.power))*(node.size.max-node.size.min))/
                       (max((cent.metric)^dist.power)-min((cent.metric)^dist.power)))

  # add the Node name
  plotcord$Node.name <- base::as.character(igraph::as_ids(V(graph)))
  
  # add node.group
  if(!is.null(node.group)) {
    plotcord$Group <- node.group
  } else {
    plotcord$Group <- plotcord$Node.name
  }

  ####*******************************####

  # get edges (pairs of node IDs)
  edgelist <- base::data.frame(igraph::get.edgelist(graph))

  # prepare a four column edge data frame with source and destination coordinates
  edges.cord <- base::data.frame(base::matrix(nrow = nrow(edgelist),
                                 ncol = 4), stringsAsFactors = FALSE)
  base::colnames(edges.cord) <- c("X1", "Y1", "X2", "Y2")
  for(i in 1:nrow(edges.cord)) {
    edges.cord$X1[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,1]),1]
    edges.cord$Y1[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,1]),2]
    edges.cord$X2[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,2]),1]
    edges.cord$Y2[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,2]),2]
  }

  # refine end positions of edges for directerd networks
  if(directed) {
    for (i in 1:nrow(edges.cord)) {

      # correct the x coordinate of arrow
      if(edges.cord$X1[i]>edges.cord$X2[i]) {
        edges.cord$X2[i] <- edges.cord$X2[i]+0.115
      } else if(edges.cord$X1[i]<edges.cord$X2[i]) {
        edges.cord$X2[i] <- edges.cord$X2[i]-0.115
      }

      # correct the y coordinate of arrow
      if(edges.cord$Y1[i]>edges.cord$Y2[i]) {
        edges.cord$Y2[i] <- edges.cord$Y2[i]+0.115
      } else if(edges.cord$Y1[i]<edges.cord$Y2[i]) {
        edges.cord$Y2[i] <- edges.cord$Y2[i]-0.115
      }
    }
  }

  # set the edge width
  if(weighted) {
    edges.cord$Weight <- igraph::E(graph)$weight
    #range normalize the weight
    edges.cord$Weight <- edge.width.min+(((edges.cord$Weight-min(edges.cord$Weight))*(edge.width.max-edge.width.min))/
                                (max(edges.cord$Weight)-min(edges.cord$Weight)))
  } else {
    edges.cord$Weight <- edge.width
  }

  ####*******************************####

  # draw the plot
  temp.plot <- ggplot2::ggplot(data = plotcord, ggplot2::aes(x = X, y = Y))

  ##***********##

    # add the edges
  if(directed) {
    temp.plot <- temp.plot +
      ggplot2::geom_segment(data=edges.cord, ggplot2::aes(x=X1, y=Y1, xend = X2, yend = Y2),
                            size = edges.cord$Weight,
                            arrow = ggplot2::arrow(angle = arrow.width,
                                          length = ggplot2::unit(arrow.length, "in"),
                                          type = "closed"),
                            colour = edge.color,
                            linetype = edge.linetype)
  } else {
    temp.plot <- temp.plot +
      ggplot2::geom_segment(data=edges.cord, ggplot2::aes(x=X1, y=Y1, xend = X2, yend = Y2),
                            size = edges.cord$Weight,
                            colour = edge.color,
                            linetype = edge.linetype)
  }

  ##***********##

  # add nodes

  # define node shapes
  node.shape <- base::as.data.frame(node.shape, stringsAsFactors = FALSE)
  for (i in 1:nrow(node.shape)) {
      if(node.shape[i,1] == "circle") {
        node.shape[i,1] <- 21
      } else if(node.shape[i,1] == "square") {
        node.shape[i,1] <- 22
      } else if(node.shape[i,1] == "diamond") {
        node.shape[i,1] <- 23
      } else if(node.shape[i,1] == "triangle") {
        node.shape[i,1] <- 24
      } else if(node.shape[i,1] == "inverted triangle") {
        node.shape[i,1] <- 25
      }
  }

  node.shape <- base::as.numeric(node.shape[,1])

  # add stroke color
  base::suppressWarnings(
  if(length(stroke.color) == 1 && stroke.color == "identical") {
    temp.plot <- temp.plot +
      ggplot2::geom_point(data = plotcord, ggplot2::aes(x = X, y = Y, colour = cent.metric),
                          shape = node.shape,
                          size = plotcord$Node.size,
                          stroke = stroke.size,
                          alpha = stroke.alpha,
                          show.legend = FALSE) +
      ggplot2::scale_color_viridis_c(option = node.color,
                                     begin = 0.15)
  } else {
    temp.plot <- temp.plot +
      ggplot2::geom_point(data = plotcord, ggplot2::aes(x = X, y = Y, color = Group),
                          shape = node.shape,
                          size = plotcord$Node.size,
                          stroke = stroke.size,
                          alpha = stroke.alpha,
                          show.legend = ifelse(length(stroke.color) == 1, FALSE, TRUE)) +
      ggplot2::scale_color_manual(values = stroke.color)
  }
  )

  # add node objects
  temp.plot <- temp.plot +
    ggplot2::geom_point(data = plotcord, ggplot2::aes(x = X, y = Y, fill = cent.metric),
                        shape = node.shape,
                        stroke = 0,
                        size = plotcord$Node.size) +

    ##***********##

    # add node color
    ggplot2::scale_fill_viridis_c(option = node.color,
                                  begin = 0.15)

    ##***********##

    # add node labels
    if(show.labels) {
      temp.plot <- temp.plot +
        ggplot2::geom_text(data = plotcord,
                           ggplot2::aes(x = X, y = Y, label=Node.name),
                           size = plotcord$Node.size*label.cex,
                           color = label.color)
    }

      ##***********##
    # expand the x and y limits
    temp.plot <- temp.plot +
    ggplot2::scale_x_continuous(expand=c(0,1)) +
    ggplot2::scale_y_continuous(expand=c(0,1))

      ##***********##

    # add title
    if(show.plot.title) {
      temp.plot <- temp.plot +
        ggplot2::ggtitle(label = plot.title)
    }

    # define title position
    if(title.position == "left") {
      title.position <- 0
    } else if(title.position == "center") {
      title.position <- 0.5
    } else if(title.position == "right") {
      title.position <- 1
    }

    title.position <- base::as.numeric(title.position)

    ##***********##

    # add theme elements

    # add main plot theme
    temp.plot <- temp.plot +
      ggplot2::theme_void() +

      # add legend specifications
      ggplot2::theme(legend.position = legend.position,
                     legend.direction = legend.direction,
                     legend.spacing.y = ggplot2::unit(0.12, "in"),
                     plot.margin = ggplot2::unit(c(0.2,0.2,0.2,0.2),
                                                 units = "cm"),
                     panel.border = ggplot2::element_blank()) +
      ggplot2::labs(fill = legend.title)

    if(boxed.legend) {
      temp.plot <- temp.plot +
        ggplot2::theme(legend.box.background = ggplot2::element_rect(color="black", size=0.5),
                       legend.margin = ggplot2::margin(c(3,3,3,3)),
                       legend.box.margin = ggplot2::margin(c(3,1,3,3)),
                       legend.box.spacing = ggplot2::unit(0, "cm"))
    }

    # add border lines
    if(show.bottom.border) {
      temp.plot <- temp.plot +
        ggplot2::theme(axis.line.x.bottom = ggplot2::element_line(color = 'black'))
    }

    if(show.left.border) {
      temp.plot <- temp.plot +
        ggplot2::theme(axis.line.y.left = ggplot2::element_line(color = 'black'))
    }

    # set title position
    if(show.plot.title) {
        temp.plot <- temp.plot +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = title.position))
      }
    return(temp.plot)
  }

  #=============================================================================
  #
  #    Code chunk 19: Visualization of ExIR results
  #
  #=============================================================================

  #' Visualization of ExIR results
  #'
  #' This function has been developed for the visualization of ExIR results. Some of the documentations
  #' of the arguments of this function have been adapted from ggplot2 package.
  #' A shiny app has also been developed for Running the ExIR model, visualization of its results as well as computational
  #' simulation of knockout and/or up-regulation of its top candidate outputs, which is accessible using
  #' the `influential::runShinyApp("ExIR")` command.
  #' You can also access the shiny app online at https://influential.erc.monash.edu/.
  #' @param exir.results An object of class \code{"ExIR_Result"} which is the output of the function \code{"exir"}.
  #' @param synonyms.table (Optional) A data frame or matrix with two columns including a column for the used feature
  #' names in the input data of the \code{"exir"} model and the other column their synonyms. Note, the original feature names should
  #' always come as the first column and the synonyms as the second one. For example, if
  #' the original feature names used for running the \code{"exir"} model are Ensembl gene
  #' symbols, you can use their HGNC synonyms in the second column to be used for the visualization of the ExIR results
  #' @param n An integer specifying the number of top candidates to be selected from each category of ExIR results (default is set to 10).
  #' @param driver.type A string specifying the type of drivers to be used for the selection of top N candidates. The possible types
  #' include \code{"combined"} (meaning both driver types), \code{"accelerator"} and \code{"decelerator"} (default is set to "combined").
  #' @param biomarker.type A string specifying the type of biomarkers to be used for the selection of top N candidates. Possible types
  #' include \code{"combined"} (meaning both biomarker types), \code{"up-regulated"} and \code{"down-regulated"} (default is set to "combined").
  #' @param show.drivers Logical scalar, whether to show Drivers or not (default is set to TRUE).
  #' @param show.biomarkers Logical scalar, whether to show Biomarkers or not (default is set to TRUE).
  #' @param show.de.mediators Logical scalar, whether to show DE-mediators or not (default is set to TRUE).
  #' @param show.nonDE.mediators Logical scalar, whether to show nonDE-mediators or not (default is set to TRUE).
  #' @param basis A string specifying the basis for the selection of top N candidates from each category of the results. Possible options include
  #' \code{"Rank"} and \code{"Adjusted p-value"} (default is set to "Rank").
  #' @param label.position By default, the labels are displayed on the top of the plot. Using label.position it is possible
  #' to place the labels on either of the four sides by setting label.position = c("top", "bottom", "left", "right").
  #' @param nrow Number of rows of the plot (default is set to 1).
  #' @param dot.size.min The size of dots with the lowest statistical significance (default is set to 2).
  #' @param dot.size.max The size of dots with the highest statistical significance (default is set to 5).
  #' @param type.color A character string or function indicating the color palette to be used for the visualization of
  #' different types of candidates. You may choose one of the Viridis palettes including "magma" (or "A"),
  #' "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E"), use a function specifying
  #' your desired palette, or manually specify the vector of colors for different types.
  #' @param stroke.size The size of stroke (border) around the dots (default is set to 1.5).
  #' @param stroke.alpha The transparency of the stroke (border) around the dots which should
  #' be a number between 0 and 1 (default is set to 1).
  #' @param dot.color.low The color to be used for the visualization of dots (features) with the lowest Z-score values (default is set to "blue").
  #' @param dot.color.high The color to be used for the visualization of dots (features) with the highest Z-score values (default is set to "red").
  #' @param legend.position The position of legends ("none", "left", "right",
  #' "bottom", "top", or two-element numeric vector). The default is set to "bottom".
  #' @param legend.direction Layout of items in legends ("horizontal" or "vertical").
  #' The default is set to "vertical".
  #' @param legends.layout Layout of different legends of the plot ("horizontal" or "vertical").
  #' The default is set to "horizontal".
  #' @param boxed.legend Logical scalar, whether to draw a box around the legend or not (default is set to TRUE).
  #' @param show.plot.title Logical scalar, whether to show the plot title or not (default is set to TRUE).
  #' @param plot.title The plot title in the string format (default is set to "auto" which automatically generates a title for the plot).
  #' @param title.position The position of title ("left", "center", or "right"). The default is set to "left".
  #' @param plot.title.size The font size of the plot title (default is set to 12).
  #' @param show.plot.subtitle Logical scalar, whether to show the plot subtitle or not (default is set to TRUE).
  #' @param plot.subtitle The plot subtitle in the string format (default is set to "auto" which automatically generates a subtitle for the plot).
  #' @param subtitle.position The position of subtitle ("left", "center", or "right"). The default is set to "left".
  #' @param y.axis.title The title of the y axis (features title). Default is set to "Features".
  #' @param show.y.axis.grid Logical scalar, whether to draw y axis grid lines (default is set to TRUE).
  #' @return A plot with the class ggplot.
  #' @keywords exir.vis
  #' @family visualization functions
  #' @seealso \code{\link[influential]{exir}}
  #' @export exir.vis
  #' @examples
  #' \dontrun{
  #' MyResults <- exir.results
  #' ExIR.plot <- exir.vis(exir.results = MyResults, n = 5)
  #' }
  exir.vis <- function(exir.results,
                       synonyms.table = NULL,
                       n = 10,
                       driver.type = "combined",
                       biomarker.type = "combined",
                       show.drivers = TRUE,
                       show.biomarkers = TRUE,
                       show.de.mediators = TRUE,
                       show.nonDE.mediators = TRUE,
                       basis = "Rank",
                       label.position = "top",
                       nrow = 1,
                       dot.size.min = 2,
                       dot.size.max = 5,
                       type.color = "viridis",
                       stroke.size = 1.5,
                       stroke.alpha = 1,
                       dot.color.low = "blue",
                       dot.color.high = "red",
                       legend.position = "bottom",
                       legend.direction = "vertical",
                       legends.layout = "horizontal",
                       boxed.legend = TRUE,
                       show.plot.title = TRUE,
                       plot.title = "auto",
                       title.position = "left",
                       plot.title.size = 12,
                       show.plot.subtitle = TRUE,
                       plot.subtitle = "auto",
                       subtitle.position = "left",
                       y.axis.title = "Feature",
                       show.y.axis.grid = TRUE) {

    if(base::identical(base::inherits(exir.results, "ExIR_Result"), FALSE)) {
      stop("The provided ExIR model-result is wrong. The exir.results should be
           the output of the exir model and have a 'ExIR_Result' class.",
           call. = FALSE)
    }

    # prepare a list for storing the results
    exir.for.plot <- base::list()

    # select top N features

    # for driver table
    if(any(base::names(exir.results) == "Driver table")) {
      top.N.driver.table <- exir.results[[which(names(exir.results) %in% "Driver table")]]
      top.N.driver.table$Feature <- base::rownames(top.N.driver.table)
      top.N.driver.table$Class <- "Driver"

      # top N combined
      if(driver.type == "combined" & nrow(top.N.driver.table)> n) {

        if(basis == "Rank") {
          top.drivers.index <- utils::head(order(top.N.driver.table$Rank),
                                           n = n)

        } else if(basis == "Adjusted p-value") {
          top.drivers.index <- utils::head(order(top.N.driver.table$P.adj),
                                           n = n)
        }
        top.N.driver.table <- top.N.driver.table[top.drivers.index,]
        top.N.driver.table$Rank <- base::rank(top.N.driver.table$Rank,
                                              ties.method = "min")

        # top N accelerator
      } else if(driver.type == "accelerator") {
        top.N.driver.table <- base::subset(top.N.driver.table,
                                           Type == "Accelerator")

          if(basis == "Rank") {
            top.drivers.index <- utils::head(order(top.N.driver.table$Rank),
                                             n = n)

          } else if(basis == "Adjusted p-value") {
            top.drivers.index <- utils::head(order(top.N.driver.table$P.adj),
                                             n = n)
          }
          top.N.driver.table <- top.N.driver.table[top.drivers.index,]
          top.N.driver.table$Rank <- base::rank(top.N.driver.table$Rank,
                                                ties.method = "min")
          
          # top N decelerator
      } else if(driver.type == "decelerator") {
        top.N.driver.table <- base::subset(top.N.driver.table,
                                           Type == "Decelerator")
          if(basis == "Rank") {
            top.drivers.index <- utils::head(order(top.N.driver.table$Rank),
                                             n = n)

          } else if(basis == "Adjusted p-value") {
            top.drivers.index <- utils::head(order(top.N.driver.table$P.adj),
                                             n = n)
          }
          top.N.driver.table <- top.N.driver.table[top.drivers.index,]
          top.N.driver.table$Rank <- base::rank(top.N.driver.table$Rank,
                                                ties.method = "min")
      }
      base::rownames(top.N.driver.table) <- NULL
      top.N.driver.table$Type[top.N.driver.table$Type == "Accelerator"] <- "Accelerator\ndriver"
      top.N.driver.table$Type[top.N.driver.table$Type == "Decelerator"] <- "Decelerator\ndriver"
      exir.for.plot <- base::append(x = exir.for.plot,
                                    values = base::list(top.N.driver.table))
    }

    ####********************####

    # for biomarker table
    if(any(base::names(exir.results) == "Biomarker table")) {
      top.N.biomarker.table <- exir.results[[which(names(exir.results) %in% "Biomarker table")]]
      top.N.biomarker.table$Feature <- base::rownames(top.N.biomarker.table)
      top.N.biomarker.table$Class <- "Biomarker"

      # top N combined
      if(biomarker.type == "combined" & nrow(top.N.biomarker.table)> n) {

        if(basis == "Rank") {
          top.biomarkers.index <- utils::head(order(top.N.biomarker.table$Rank),
                                              n = n)

        } else if(basis == "Adjusted p-value") {
          top.biomarkers.index <- utils::head(order(top.N.biomarker.table$P.adj),
                                              n = n)
        }
        top.N.biomarker.table <- top.N.biomarker.table[top.biomarkers.index,]
        top.N.biomarker.table$Rank <- base::rank(top.N.biomarker.table$Rank,
                                                 ties.method = "min")

        # top N up-regulated
      } else if(biomarker.type == "up-regulated") {
        top.N.biomarker.table <- base::subset(top.N.biomarker.table,
                                              Type == "Up-regulated")
          if(basis == "Rank") {
            top.biomarkers.index <- utils::head(order(top.N.biomarker.table$Rank),
                                                n = n)

          } else if(basis == "Adjusted p-value") {
            top.biomarkers.index <- utils::head(order(top.N.biomarker.table$P.adj),
                                                n = n)
          }
          top.N.biomarker.table <- top.N.biomarker.table[top.biomarkers.index,]
          top.N.biomarker.table$Rank <- base::rank(top.N.biomarker.table$Rank,
                                                ties.method = "min")

      } else if(biomarker.type == "down-regulated") {
        top.N.biomarker.table <- base::subset(top.N.biomarker.table,
                                              Type == "Down-regulated")
          if(basis == "Rank") {
            top.biomarkers.index <- utils::head(order(top.N.biomarker.table$Rank),
                                                n = n)

          } else if(basis == "Adjusted p-value") {
            top.biomarkers.index <- utils::head(order(top.N.biomarker.table$P.adj),
                                                n = n)
          }
          top.N.biomarker.table <- top.N.biomarker.table[top.biomarkers.index,]
          top.N.biomarker.table$Rank <- base::rank(top.N.biomarker.table$Rank,
                                                   ties.method = "min")

      }
      base::rownames(top.N.biomarker.table) <- NULL
      top.N.biomarker.table$Type[top.N.biomarker.table$Type == "Up-regulated"] <- "Up-regulated\nbiomarker"
      top.N.biomarker.table$Type[top.N.biomarker.table$Type == "Down-regulated"] <- "Down-regulated\nbiomarker"

      exir.for.plot <- base::append(x = exir.for.plot,
                                    values = base::list(top.N.biomarker.table))
    }

    ####********************####

    # for nonDE-mediator table
    if(any(base::names(exir.results) == "nonDE-mediator table")) {
      top.N.nonDE.mediator.table <- exir.results[[which(names(exir.results) %in% "nonDE-mediator table")]]
      top.N.nonDE.mediator.table$Type <- "nonDE-mediator"
      top.N.nonDE.mediator.table$Feature <- base::rownames(top.N.nonDE.mediator.table)
      top.N.nonDE.mediator.table$Class <- "nonDE-mediator"

      # top N combined

        if(basis == "Rank") {
          top.nonDE.mediators.index <- utils::head(order(top.N.nonDE.mediator.table$Rank),
                                                   n = n)

        } else if(basis == "Adjusted p-value") {
          top.nonDE.mediators.index <- utils::head(order(top.N.nonDE.mediator.table$P.adj),
                                                   n = n)
        }
        top.N.nonDE.mediator.table <- top.N.nonDE.mediator.table[top.nonDE.mediators.index,]
        top.N.nonDE.mediator.table$Rank <- base::rank(top.N.nonDE.mediator.table$Rank,
                                                 ties.method = "min")
      base::rownames(top.N.nonDE.mediator.table) <- NULL
      exir.for.plot <- base::append(x = exir.for.plot,
                                    values = base::list(top.N.nonDE.mediator.table))
    }

    ####********************####

    # for DE-mediator table
    if(any(base::names(exir.results) == "DE-mediator table")) {
      top.N.DE.mediator.table <- exir.results[[which(names(exir.results) %in% "DE-mediator table")]]
      top.N.DE.mediator.table$Type <- "DE-mediator"
      top.N.DE.mediator.table$Feature <- base::rownames(top.N.DE.mediator.table)
      top.N.DE.mediator.table$Class <- "DE-mediator"

      # top N combined

        if(basis == "Rank") {
          top.DE.mediators.index <- utils::head(order(top.N.DE.mediator.table$Rank),
                                                n = n)

        } else if(basis == "Adjusted p-value") {
          top.DE.mediators.index <- utils::head(order(top.N.DE.mediator.table$P.adj),
                                                n = n)
        }
        top.N.DE.mediator.table <- top.N.DE.mediator.table[top.DE.mediators.index,]
        top.N.DE.mediator.table$Rank <- base::rank(top.N.DE.mediator.table$Rank,
                                                      ties.method = "min")
      base::rownames(top.N.DE.mediator.table) <- NULL
      exir.for.plot <- base::append(x = exir.for.plot,
                                    values = base::list(top.N.DE.mediator.table))
    }

    # combine the results for plotting
    exir.for.plot <- base::Reduce(function(dtf1, dtf2) base::merge(dtf1, dtf2,
                                                                   all = TRUE,
                                                                   all.x = TRUE),
                                  exir.for.plot)

    # correct the features names
    if(!is.null(synonyms.table)) {
      synonyms.table <- base::as.data.frame(synonyms.table, stringsAsFactors = FALSE)
      synonyms.index <- base::match(exir.for.plot$Feature,
                                    synonyms.table[,1])
      exir.for.plot$Feature <- synonyms.table[synonyms.index,2]
    }

    # correct the Type levels
    driver.levels <- base::unique(exir.for.plot$Type[base::grep("driver", exir.for.plot$Type)])
    biomarker.levels <- base::unique(exir.for.plot$Type[base::grep("biomarker", exir.for.plot$Type)])
    mediators.levels <- base::unique(exir.for.plot$Type[base::seq_along(rownames(exir.for.plot))[-c(base::grep("driver", exir.for.plot$Type),
                                                                                                    base::grep("biomarker", exir.for.plot$Type))]])
    exir.for.plot$Type <- base::factor(exir.for.plot$Type,
                                       levels = c(driver.levels,
                                                  biomarker.levels,
                                                  mediators.levels))

    # correct the Class levels
    mediators.class.levels <- base::unique(exir.for.plot$Class[base::seq_along(rownames(exir.for.plot))[-c(base::grep("Driver", exir.for.plot$Class),
                                                                                                           base::grep("Biomarker", exir.for.plot$Class))]])
    exir.for.plot$Class <- base::factor(exir.for.plot$Class,
                                        levels = c("Driver",
                                                   "Biomarker",
                                                   mediators.class.levels))

    # remove undesired classes
    if(isFALSE(show.drivers) & any(exir.for.plot$Class == "Driver")) {
      exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "Driver")),]
    }

    if(isFALSE(show.biomarkers) & any(exir.for.plot$Class == "Biomarker")) {
      exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "Biomarker")),]
    }

    if(isFALSE(show.de.mediators) & any(exir.for.plot$Class == "DE-mediator")) {
      exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "DE-mediator")),]
    }

    if(isFALSE(show.nonDE.mediators) & any(exir.for.plot$Class == "nonDE-mediator")) {
      exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "nonDE-mediator")),]
    }

    # correct the P.adj to be used as the dot size

    exir.for.plot$P.value[is.nan(exir.for.plot$P.value)] <- 1
    exir.for.plot$P.adj[is.nan(exir.for.plot$P.adj)] <- 1
    exir.for.plot$P.value[is.na(exir.for.plot$P.value)] <- 1
    exir.for.plot$P.adj[is.na(exir.for.plot$P.adj)] <- 1

    if(min(exir.for.plot$P.adj)==0) {

      #range normalize the primitive P.adj
      temp.min_P.adj <- base::sort(base::unique(exir.for.plot$P.adj))[2]

      exir.for.plot$P.adj <- temp.min_P.adj+
        (((exir.for.plot$P.adj-min(exir.for.plot$P.adj))*(max(exir.for.plot$P.adj)-temp.min_P.adj))/
           (max(exir.for.plot$P.adj)-min(exir.for.plot$P.adj)))
    }

    # Set the P.adj based on min and max arguments
    if(length(unique(exir.for.plot$P.adj)) == 1) {
      exir.for.plot$P.adj <- mean(c(dot.size.min, dot.size.max))

    } else {
      exir.for.plot$P.adj <- dot.size.min+(((-log10(exir.for.plot$P.adj)-min(-log10(exir.for.plot$P.adj)))*(dot.size.max-dot.size.min))/
                                             (max(-log10(exir.for.plot$P.adj))-min(-log10(exir.for.plot$P.adj))))
    }

    # Correct the levels of Features based on each class
    visClassLength <- base::length(base::unique(exir.for.plot$Class))

    visFeatureLevels <- base::unique(base::as.character(base::unlist(base::sapply(X = 1:visClassLength,
                                                               FUN = function(i) {
                                                                 exir.for.plot$Feature[which(exir.for.plot$Class %in% base::unique(exir.for.plot$Class)[i])]
                                                               }))))

    exir.for.plot$Feature <- base::factor(exir.for.plot$Feature,
                                          levels = visFeatureLevels)

    ####*******************************####

    # draw the plot
    temp.exir.plot <- ggplot2::ggplot(data = exir.for.plot,
                                      ggplot2::aes(x = Rank, y = Feature)) +

      # add node objects
      ggplot2::geom_point(ggplot2::aes(fill = Z.score,
                              colour = Type,
                              size = P.adj),
                          shape = 21,
                          stroke = 0,
                          alpha = 1) +

      ##***********##

      # add color of Type
      ggplot2::geom_point(ggplot2::aes(colour = Type,
                              size = P.adj),
                          shape = 21,
                          stroke = stroke.size,
                          alpha = stroke.alpha)


    ##***********##

    # add node and stroke colors
    if(base::inherits(type.color, "character") & base::length(type.color) == 1) {
    if(base::length(base::grep(type.color,
                               c("magma", "inferno",
                                 "plasma", "viridis", "cividis",
                                 "A", "B", "C", "D", "E"))) == 1) {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::scale_colour_viridis_d(option = type.color)
    } else {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::scale_colour_manual(values = type.color)
    }
    } else {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::scale_colour_manual(values = type.color)
    }

    temp.exir.plot <- temp.exir.plot +
      ggplot2::scale_fill_gradient(name = "Z-score",
                                   low = dot.color.low,
                                   high = dot.color.high) +

      ##***********##

      # correct size identity inside of aes
      ggplot2::scale_size_identity(guide = "legend") +


      ##***********##

      # add y axis title
      ggplot2::ylab(y.axis.title) +


      ##***********##

      # add facets
      ggplot2::facet_wrap(. ~ Class,
                          scales = "free_x",
                          strip.position = label.position,
                          nrow = nrow) +


      ##***********##

      # add theme elements

      ggplot2::theme_bw()

    # set the x axis breaks
    rank.uniqueness <- base::vector(mode = "integer")
    for(i in base::unique(exir.for.plot$Class)) {
      rank.uniqueness <- base::append(x = rank.uniqueness, values =
                                        length(unique(exir.for.plot$Rank[exir.for.plot$Class == i]))
      )
    }
    by.x_continuous <- base::ifelse(any(rank.uniqueness == 1), 1, base::round(n/5))

    # set the x axis numbers and start
    x.axis.numbers <- function(x) {
      if(by.x_continuous %% 2 == 0 & n > 7) {
        base::seq(2, max(x), by = by.x_continuous)
      } else if(n > 7) {
        base::seq(1, max(x), by = by.x_continuous)
      } else {
        base::seq(1, max(x), by = 1)
      }
    }

    # set the x axis limits
    x.axis.limits <- function(x) {
      c((min(x)-(n/50)), (max(x)+(n/50)))
    }

    temp.exir.plot <- temp.exir.plot +
      ggplot2::scale_x_continuous(breaks = x.axis.numbers,
                                  limits = x.axis.limits) +

      ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                     legend.box = legends.layout,
                     legend.position = legend.position,
                     legend.direction = legend.direction) +

      ggplot2::guides(colour = ggplot2::guide_legend(keyheight = 1.5),
                      fill = ggplot2::guide_colorbar(frame.colour = "black",
                                            barwidth = 1.5),
                      size = ggplot2::guide_legend(title = "Statistical\nsignificance"))

    if(boxed.legend) {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::theme(legend.box.background = ggplot2::element_rect(color="black", size=0.5),
                       legend.margin = ggplot2::margin(c(3,3,3,3)),
                       legend.box.margin = ggplot2::margin(c(3,1,3,3)))
    }

    ##***********##

    # add title

    if(plot.title == "auto") {
      plot.title <- "ExIR model-based prioritized features"
    }

    if(show.plot.title) {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::labs(title = plot.title)
    }

    # add subtitle

    if(plot.subtitle == "auto") {
      plot.subtitle <- base::paste(basis,
                                   "-based selection of top ",
                                   n,
                                   " candidates", sep = "")
    }

    if(show.plot.subtitle) {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::labs(subtitle = plot.subtitle)
    }

    # define title position
    if(title.position == "left") {
      title.position <- 0
    } else if(title.position == "center") {
      title.position <- 0.5
    } else if(title.position == "right") {
      title.position <- 1
    }

    title.position <- base::as.numeric(title.position)

    # define subtitle position
    if(subtitle.position == "left") {
      subtitle.position <- 0
    } else if(subtitle.position == "center") {
      subtitle.position <- 0.5
    } else if(subtitle.position == "right") {
      subtitle.position <- 1
    }

    subtitle.position <- base::as.numeric(subtitle.position)

    title.size <- plot.title.size - 2

    # set title position
    if(show.plot.title) {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::theme(plot.title = ggplot2::element_text(size = title.size,
                                                          hjust = title.position))
    }

    # set subtitle position
    if(show.plot.subtitle) {
      subtitle.size <- title.size - 2
      temp.exir.plot <- temp.exir.plot +
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = subtitle.size,
                                                             hjust = subtitle.position))
    }

    ##***********##

    # Set the order of legends

    temp.exir.plot <- temp.exir.plot + ggplot2::guides(color = ggplot2::guide_legend(order = 1),
                                                       size = ggplot2::guide_legend(order = 2))

    ##***********##

    if(base::identical(show.y.axis.grid, FALSE)) {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
    }

    return(temp.exir.plot)
  }

#=============================================================================
#
#    Code chunk 20: Computational manipulation of cells
#
#=============================================================================

  #' Computational manipulation of cells
  #'
  #' This function works based on the SIRIR (SIR-based Influence Ranking) model and could be applied on the
  #' output of the ExIR model or any other independent association network. For feature (gene/protein/etc.)
  #' knockout the SIRIR model is used to remove the feature from the network and assess its impact on the
  #' flow of information (signaling) within the network. On the other hand, in case of up-regulation a node similar
  #' to the desired node is added to the network with exactly the same connections (edges) as of the original node.
  #' Next, the SIRIR model is used to evaluate the difference in the flow of information/signaling after adding (up-regulating)
  #' the desired feature/node compared with the original network. In case you are applying this function on the output of
  #' ExIR model, you may note that as the gene/protein knockout would impact on the integrity of the under-investigation network
  #' as well as the networks of other overlapping biological processes/pathways, it is recommended to select those features that
  #' simultaneously have the highest (most significant) ExIR-based rank and lowest knockout rank. In contrast, as the up-regulation
  #' would not affect the integrity of the network, you may select the features with highest (most significant) ExIR-based and
  #' up-regulation-based ranks.
  #' A shiny app has also been developed for Running the ExIR model, visualization of its results as well as computational
  #' simulation of knockout and/or up-regulation of its top candidate outputs, which is accessible using
  #' the `influential::runShinyApp("ExIR")` command.
  #' You can also access the shiny app online at https://influential.erc.monash.edu/.
  #' @param exir_output The output of the ExIR model (optional).
  #' @param graph A graph (network) of the igraph class (not required if the exir_output is inputted).
  #' @param ko_vertices A vector of desired vertices/features to knockout. Default is set to V(graph) meaning to assess
  #' the knockout of all vertices/features.
  #' @param upregulate_vertices A vector of desired vertices/features to up-regulate. Default is set to V(graph) meaning to assess
  #' the up-regulation of all vertices/features.
  #' @param beta Non-negative scalar corresponding to the SIRIR model. The rate of infection of an individual that is susceptible
  #' and has a single infected neighbor. The infection rate of a susceptible individual with n
  #' infected neighbors is n times beta. Formally this is the rate parameter of an exponential
  #' distribution.
  #' @param gamma Positive scalar corresponding to the SIRIR model. The rate of recovery of an infected individual.
  #' Formally, this is the rate parameter of an exponential distribution.
  #' @param no.sim Integer scalar corresponding to the SIRIR model. The number of simulation runs to perform SIR model on for the
  #' original network as well perturbed networks generated by leave-one-out technique.
  #' You may choose a different no.sim based on the available memory on your system.
  #' @param node_verbose Logical; whether the process of Parallel Socket Cluster creation should be printed (default is FALSE).
  #' @param loop_verbose Logical; whether the accomplishment of the evaluation of network nodes in each loop should be printed (default is TRUE).
  #' @param ncores Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of 
  #' cores to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").
  #' @param seed A single value, interpreted as an integer to be used for random number generation.
  #' @return Depending on the input data, a list including one to three data frames of knockout/up-regulation rankings.
  #' @keywords comp_manipulate
  #' @family integrative ranking functions
  #' @seealso \code{\link[influential]{exir}}, \code{\link[influential]{sirir}},
  #' and \code{\link[igraph]{sir}} for a complete description on SIR model
  #' @export comp_manipulate
  #' @examples
  #' \dontrun{
  #' set.seed(1234)
  #' My_graph <- igraph::sample_gnp(n=50, p=0.05)
  #' GraphVertices <- V(My_graph)
  #' Computational_manipulation <- comp_manipulate(graph = My_graph, beta = 0.5,
  #'                                               gamma = 1, no.sim = 10, seed = 1234)
  #'                                               }
  #' @importFrom igraph vcount as_ids sir
  #' @importFrom foreach %dopar% foreach

  comp_manipulate <- function(exir_output = NULL,
                              graph = NULL,
                              ko_vertices = igraph::V(graph),
                              upregulate_vertices = igraph::V(graph),
                              beta = 0.5,
                              gamma = 1,
                              no.sim = 100,
                              node_verbose = FALSE,
                              loop_verbose = TRUE,
                              ncores = "default", 
                              seed = 1234) {

    ##**************************##
    # Take care of input graph
    if(!is.null(exir_output) && inherits(exir_output, "ExIR_Result")) {
      graph <- exir_output$Graph
    }

    ##**************************##
    # Over-expression function

    overexpr <- function(graph, vertices = upregulate_vertices, beta = beta, gamma = gamma,
                          no.sim = no.sim, loop_verbose = loop_verbose, 
                         ncores = ncores,
                         node_verbose = node_verbose, seed = seed) {
      
      suppressWarnings({
        
        # Make clusters for parallel processing
        cl <- parallel::makeCluster(ifelse(ncores == "default", parallel::detectCores() - 1, ncores), 
                                    outfile=ifelse(node_verbose, "", "NULL"))
        doParallel::registerDoParallel(cl)
      
        temp.loocr.table <- data.frame(difference.value = rep(NA, length(vertices)), rank = rep(NA, length(vertices)))
        if (inherits(vertices, "character")) {
          rownames(temp.loocr.table) <- vertices
        } else if (inherits(vertices, "igraph.vs")) {
          rownames(temp.loocr.table) <- igraph::as_ids(vertices)
        }
        
        set.seed(seed)
        all.included.spread <- igraph::sir(graph = graph, beta = beta,
                                           gamma = gamma, no.sim = no.sim)
        all.mean_spread <- foreach(i = 1:length(all.included.spread), .combine = "c", .verbose = loop_verbose) %dopar% {
          max(all.included.spread[[i]]$NR)
        }
        all.mean_spread <- mean(all.mean_spread)
        
        results <- foreach(s = 1:length(vertices), .combine = "c", .verbose = loop_verbose) %dopar% {
          all.edges <- as.data.frame(igraph::as_edgelist(graph))
          all.desired.edges.index <- unlist(apply(X = all.edges, MARGIN = 2,
                                                  FUN = function(i) which(i %in% rownames(temp.loocr.table)[s])))
          
          all.edges <- all.edges[c(all.desired.edges.index),]
          
          all.edges.from.indices <- unlist(lapply(X = all.edges[,1],
                                                  FUN = function(i) which(igraph::as_ids(V(graph)) %in% i)))
          
          all.edges.to.indices <- unlist(lapply(X = all.edges[,2],
                                                FUN = function(i) which(igraph::as_ids(V(graph)) %in% i)))
          
          all.edges.desired.indices <- data.frame(from = all.edges.from.indices, to = all.edges.to.indices)
          
          temp.graph <- igraph::add_vertices(graph = graph, nv = 1,
                                             name = paste0(rownames(temp.loocr.table)[s], "Fold1"))
          
          desired.vertex.index <- match(rownames(temp.loocr.table)[s], igraph::as_ids(V(temp.graph)))
          
          desired.vertexFold1.index <- match(paste0(rownames(temp.loocr.table)[s], "Fold1"),
                                             igraph::as_ids(V(temp.graph))
          )
          
          all.edges.desired.indices[which(all.edges.desired.indices[,1] == desired.vertex.index),1] <- desired.vertexFold1.index
          all.edges.desired.indices[which(all.edges.desired.indices[,2] == desired.vertex.index),2] <- desired.vertexFold1.index
          all.edges.desired.indices <- apply(X = all.edges.desired.indices, MARGIN = 1,
                                             FUN = function(i) paste(i[1], i[2], sep = ","))
          
          all.edges.desired.indices <- paste(all.edges.desired.indices, collapse = ",")
          
          all.edges.desired.indices <- as.numeric(unlist(strsplit(x = all.edges.desired.indices, split = ",")))
          
          temp.graph <- igraph::add_edges(graph = temp.graph, edges = all.edges.desired.indices)
          
          set.seed(seed)
          loocr.spread <- igraph::sir(graph = temp.graph, beta = beta,
                                      gamma = gamma, no.sim = no.sim)
          loocr.mean_spread <- foreach(h = 1:length(loocr.spread), .combine = "c", .verbose = loop_verbose) %dopar% {
            max(loocr.spread[[h]]$NR)
          }
          loocr.mean_spread <- mean(loocr.mean_spread)
          
          loocr.mean_spread - all.mean_spread
        }
        
        temp.loocr.table$difference.value <- results
        
        # Stop the parallel backend
        parallel::stopCluster(cl)
        
        temp.loocr.table$rank <- rank(-temp.loocr.table$difference.value, ties.method = "min")
      
      })
      
      return(temp.loocr.table)
    }

    ##**************************##
    # Knockout results
    if(!is.null(ko_vertices)) {
      base::suppressWarnings(
      ko_results <- influential::sirir(
        graph = graph,
        vertices = ko_vertices,
        beta = beta,
        gamma = gamma,
        no.sim = no.sim,
        loop_verbose = loop_verbose, 
        node_verbose = node_verbose,
        ncores = ncores,
        seed = seed
      )
      )

      ko_results <- cbind(Feature_name = rownames(ko_results),
                          ko_results,
                          Manipulation_type = "Knockout")

      rownames(ko_results) <- NULL
      colnames(ko_results)[3] <- "Rank"

      # correct the orders
      ko_results <- ko_results[order(ko_results[,3]),]

    } else {ko_results <- NULL}

    ##**************************##

    # Over-expression results
    if(!is.null(upregulate_vertices)) {
      base::suppressWarnings(
      overexpr_results <- overexpr(
        graph = graph,
        vertices = upregulate_vertices,
        beta = beta,
        gamma = gamma,
        no.sim = no.sim,
        loop_verbose = loop_verbose, 
        node_verbose = node_verbose,
        ncores = ncores,
        seed = seed
      )
      )

      overexpr_results <- cbind(Feature_name = rownames(overexpr_results),
                                overexpr_results,
                                Manipulation_type = "Up-regulation")

      rownames(overexpr_results) <- NULL
      colnames(overexpr_results)[3] <- "Rank"

      # correct the orders
      overexpr_results <- overexpr_results[order(overexpr_results[,3]),]

    } else {overexpr_results <- NULL}

    ##**************************##

    # Combined results
    if(!is.null(ko_results) & !is.null(overexpr_results)) {
      combined_results <- rbind(ko_results,
                                overexpr_results)

      # Correct the combined ranks
      combined_results[,3] <- rank(-combined_results[,2],
                                   ties.method = "min")

      # correct the orders
      combined_results <- combined_results[order(combined_results[,3]),]

    } else {
      combined_results <- NULL
    }

    ##**************************##

    # Correct results data frames

    ko_results <- ko_results[,-2]
    overexpr_results <- overexpr_results[,-2]
    combined_results <- combined_results[,-2]

    ##**************************##

    # Final results
    final.results <- list(Knockout = ko_results,
                          Up_regulation = overexpr_results,
                          Combined = combined_results
    )
    if(sum(sapply(final.results, is.null)) == 3) {
      cat("You should input the name of at least a vertix (feature/gene/etc)
        in the 'ko_vertices' or 'upregulate_vertices' argument.")
    } else {
      non_null_index <- which(sapply(final.results, function(i) (!is.null(i))))
      final.results <- final.results[non_null_index]
      names(final.results) <- names(non_null_index)
      return(final.results)
    }
  }

#=============================================================================
#
#    Code chunk 21: Fast correlation and mutual rank analysis
#
#=============================================================================

  #' Fast correlation and mutual rank analysis
  #'
  #' This function calculates Pearson/Spearman correlations between all pairs of features in a matrix/dataframe much faster than the base R cor function.
  #' It is also possible to simultaneously calculate mutual rank (MR) of correlations as well as their p-values and adjusted p-values.
  #' Additionally, this function can automatically combine and flatten the result matrices. Selecting correlated features using an MR-based threshold
  #' rather than based on their correlation coefficients or an arbitrary p-value is more efficient and accurate in inferring
  #' functional associations in systems, for example in gene regulatory networks.
  #' @param data a numeric dataframe/matrix (features on columns and samples on rows).
  #' @param use The NA handler, as in R's cov() and cor() functions. Options are "everything", "all.obs", and "complete.obs".
  #' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson" or "spearman" (default).
  #' @param mutualRank logical, whether to calculate mutual ranks of correlations or not.
  #' @param mutualRank_mode a character string indicating whether to rank based on "signed" or "unsigned" (default) correlation values. 
  #' In the "unsigned" mode, only the level of a correlation value is important and not its sign (the function ranks the absolutes of correlations). 
  #' Options are "unsigned", and "signed".
  #' @param pvalue logical, whether to calculate p-values of correlations or not.
  #' @param adjust p-value correction method (when pvalue = TRUE), a character string including any of "BH" (default),
  #' "bonferroni", "holm", "hochberg", "hommel", or "none".
  #' @param flat logical, whether to combine and flatten the result matrices or not.
  #' @return Depending on the input data, a dataframe or list including cor (correlation coefficients),
  #' mr (mutual ranks of correlation coefficients), p (p-values of correlation coefficients), and p.adj (adjusted p-values).
  #' @keywords fcor
  #' @seealso \code{\link[coop]{pcor}}, \code{\link[stats]{p.adjust}},
  #' and \code{\link[influential]{graph_from_data_frame}}
  #' @export fcor
  #' @examples
  #' \dontrun{
  #' set.seed(1234)
  #' data <- datasets::attitude
  #' cor <- fcor(data = data)
  #' }

  fcor <- function(data,
                   use = "everything",
                   method = "spearman",
                   mutualRank = TRUE,
                   mutualRank_mode = "unsigned",
                   pvalue = FALSE,
                   adjust = "BH",
                   flat = TRUE) {
    
    # Define and compile the rank_matrix function inline
    Rcpp::cppFunction('
  NumericMatrix rank_matrix(NumericMatrix mat, bool descending = true, bool use_abs = false) {
    int nrows = mat.nrow();
    int ncols = mat.ncol();
    NumericMatrix ranks(nrows, ncols);
    auto get_val = [use_abs](double x) { return use_abs ? fabs(x) : x; };
    auto comp = [descending](double va, double vb) { 
      if (descending) return va > vb;
      return va < vb;
    };
    for(int i = 0; i < nrows; i++) {
      std::vector<double> row_vec(ncols);
      for(int j=0; j<ncols; j++) row_vec[j] = mat(i,j);
      std::vector<size_t> idx(ncols);
      std::iota(idx.begin(), idx.end(), 0);
      std::sort(idx.begin(), idx.end(), 
        [&get_val, &comp, &row_vec](size_t a, size_t b){
          double va = get_val(row_vec[a]);
          double vb = get_val(row_vec[b]);
          if(va == vb) return a < b;
          return comp(va, vb);
        });
      size_t k = 0;
      while(k < ncols) {
        size_t start = k;
        double val = get_val(row_vec[idx[k]]);
        while(k < ncols && get_val(row_vec[idx[k]]) == val) ++k;
        size_t end = k - 1;
        double avg_rank = (static_cast<double>(start + 1) + static_cast<double>(end + 1)) / 2.0;
        for(size_t j = start; j <= end; ++j) {
          size_t id = idx[j];
          ranks(i, id) = avg_rank;
        }
      }
    }
    return ranks;
  }')
    
    # Define and compile the flatten_cor_matrix function inline
    Rcpp::cppFunction('
  List flatten_cor_matrix(NumericMatrix cormat, Nullable<NumericMatrix> mrmat = R_NilValue, 
                          Nullable<NumericMatrix> pmat = R_NilValue, Nullable<NumericMatrix> padjmat = R_NilValue,
                          CharacterVector row_names = CharacterVector::create()) {
    int m = cormat.nrow();
    bool has_mr = mrmat.isNotNull();
    bool has_p = pmat.isNotNull();
    bool has_pa = padjmat.isNotNull();
    bool has_names = row_names.size() == m;
    size_t num = (size_t)(m) * (m - 1LL) / 2;
    CharacterVector rows, cols;
    if(has_names) {
      rows = CharacterVector(num);
      cols = CharacterVector(num);
    }
    NumericVector cors(num);
    NumericVector mrs, ps, pas;
    if(has_mr) mrs = NumericVector(num);
    if(has_p) ps = NumericVector(num);
    if(has_pa) pas = NumericVector(num);
    NumericMatrix mr = has_mr ? NumericMatrix(mrmat) : NumericMatrix();
    NumericMatrix p = has_p ? NumericMatrix(pmat) : NumericMatrix();
    NumericMatrix pa = has_pa ? NumericMatrix(padjmat) : NumericMatrix();
    size_t cnt = 0;
    for(int i = 0; i < m-1; i++) {
      for(int j = i+1; j < m; j++) {
        if(has_names) {
          rows[cnt] = row_names[i];
          cols[cnt] = row_names[j];
        }
        cors[cnt] = cormat(i, j);
        if(has_mr) mrs[cnt] = mr(i, j);
        if(has_p) ps[cnt] = p(i, j);
        if(has_pa) pas[cnt] = pa(i, j);
        cnt++;
      }
    }
    List res;
    if(has_names) {
      res["row"] = rows;
      res["column"] = cols;
    }
    res["cor"] = cors;
    if(has_mr) res["mr"] = mrs;
    if(has_p) res["p"] = ps;
    if(has_pa) res["p.adj"] = pas;
    return res;
  }')
    
    #######################
    
    # Preserve dimnames
    var_names <- colnames(data)
    obs_names <- rownames(data)
    
    if(method == "spearman") {
      # Rank columns using C++ (ascending, no abs)
      data_t <- t(data)
      ranked_t <- rank_matrix(data_t, descending = FALSE, use_abs = FALSE)
      data <- t(ranked_t)
      colnames(data) <- var_names
      rownames(data) <- obs_names
    }
    
    #######################
    
    # Set initial NULL values
    r = NULL
    mutR = NULL
    p = NULL
    pa = NULL
    
    #######################
    
    # Perform correlation analysis using base::cor (faster and equivalent for pearson)
    r <- cor(x = data, use = use, method = "pearson")
    
    if(pvalue) {
      
      # Calculate n required for p-value measurement
      n <- nrow(data)
      
      # Calculate t required for p-value measurement
      t <- (r * sqrt(n - 2)) / sqrt(1 - r^2)
      
      # Calculate p-value
      p <- -2 * expm1(stats::pt(abs(t), (n - 2), log.p = TRUE))
      p[p > 1 | is.nan(p)] <- 1
      
      # Calculate adjusted p-value
      if (adjust != "none") {
        pa <- stats::p.adjust(p, adjust)
      }
    }
    
    # Calculate Mutual Rank
    if(mutualRank) {
      use_abs <- (mutualRank_mode == "unsigned")
      r_rank <- rank_matrix(r, descending = TRUE, use_abs = use_abs)
      dimnames(r_rank) <- dimnames(r)  # Add back dimnames
      mutR <- sqrt(r_rank * t(r_rank))
    }
    
    #######################
    
    # Flatten the results if requested
    if(flat) {
      flt_list <- flatten_cor_matrix(r, if(mutualRank) mutR else NULL,
                                     if(pvalue) p else NULL, if(pvalue && adjust != "none") pa else NULL,
                                     rownames(r))
      result <- data.frame(flt_list, stringsAsFactors = FALSE, row.names = NULL)
    } else {
      result <- list(r = r,
                     mr = mutR,
                     p = p,
                     p.adj = pa)
    }
    
    class(result) <- c(class(result), "fcor", "influential")
    return(result)
  }
  
#=============================================================================
#
#    Code chunk ∞: Required global variables
#
#=============================================================================

  utils::globalVariables(c("Feature",
                           "Node.name",
                           "P.adj",
                           "Rank",
                           "Type",
                           "i",
                           "X",
                           "X1",
                           "X2",
                           "Y",
                           "Y1",
                           "Y2",
                           "Z.score",
                           "s",
                           "h"
                           ))
