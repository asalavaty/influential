#=============================================================================
#
#    Hubness score
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
