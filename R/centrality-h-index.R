#=============================================================================
#
#   H-index
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
    message("Precomputing neighborhood sizes of all nodes")
  }
  
  # Names of all graph vertices
  graph_node_names <- as.character(igraph::as_ids(igraph::V(graph)))
  
  # Precompute the ego_size(order = 1) - 1
  all_neighbor_sizes <- igraph::ego_size(
    graph = graph,
    nodes = igraph::V(graph),
    mode = mode,
    order = 1
  ) - 1
  
  names(all_neighbor_sizes) <- graph_node_names
  
  if(verbose) {
    message("Getting the first neighbors of selected nodes")
  }
  
  first.neighbors <- igraph::neighborhood(
    graph = graph,
    nodes = vertices,
    mode = mode
  )
  
  if(verbose) {
    message("Calculating the H-index")
  }
  
  # internal helper
  .calc_h_index <- function(x) {
    
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    
    if(length(x) == 0) {
      return(0L)
    }
    
    if(max(x) == 0) {
      return(0L)
    }
    
    x <- sort(x, decreasing = TRUE)
    
    as.integer(utils::tail(which(x >= seq_along(x)), 1))
  }
  
  hindex <- vapply(
    first.neighbors,
    function(nbrs) {
      
      nbr_names <- as.character(igraph::as_ids(nbrs[-1]))
      
      if(length(nbr_names) == 0) {
        return(0L)
      }
      
      temp.neighbors.size <- all_neighbor_sizes[nbr_names]
      
      .calc_h_index(temp.neighbors.size)
    },
    integer(1)
  )
  
  # Getting the names of vertices
  if(inherits(vertices, "igraph.vs")) {
    node.names <- as.character(igraph::as_ids(vertices))
  } else {
    node.names <- as.character(vertices)
  }
  
  names(hindex) <- node.names
  
  return(hindex)
}
