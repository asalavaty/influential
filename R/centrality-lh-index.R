#=============================================================================
#
#    Local H-index (LH-index)
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
lh_index <- function(graph, vertices = V(graph), mode = "all",
                     ncores = "default", verbose = FALSE) {
  
  if(verbose) {
    message("Getting the first neighbors of each node")
  }
  
  # Getting the first neighbors of each node.
  # This preserves the original igraph::neighborhood() semantics.
  first.neighbors <- igraph::neighborhood(
    graph = graph,
    nodes = vertices,
    mode = mode
  )
  
  if(verbose) {
    message("Calculating H-index for all graph nodes once")
  }
  
  # Original function repeatedly called h_index() on neighbourhoods.
  # To preserve results while avoiding repeated recomputation, calculate h_index
  # once for all graph vertices, then sum within each neighbourhood.
  all_h_index <- h_index(
    graph = graph,
    vertices = igraph::V(graph),
    mode = mode,
    verbose = FALSE
  )
  
  graph_vertex_names <- igraph::as_ids(igraph::V(graph))
  names(all_h_index) <- graph_vertex_names
  
  if(verbose) {
    message("Calculating local H-index")
  }
  
  lhindex <- vapply(
    first.neighbors,
    function(vs) {
      ids <- igraph::as_ids(vs)
      sum(all_h_index[ids], na.rm = TRUE)
    },
    numeric(1)
  )
  
  # Getting the names of vertices
  if(inherits(vertices, "igraph.vs")) {
    node.names <- as.character(igraph::as_ids(vertices))
  } else {
    node.names <- as.character(vertices)
  }
  
  names(lhindex) <- node.names
  
  return(lhindex)
}
