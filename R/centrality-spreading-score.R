#=============================================================================
#
#    Spreading score
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
