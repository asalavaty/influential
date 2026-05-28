#=============================================================================
#
#    Neighborhood connectivity
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
neighborhood.connectivity <- function(graph, vertices = V(graph),
                                      mode = "all", verbose = FALSE) {
  
  # Getting the names of vertices
  if(inherits(vertices, "igraph.vs")) {
    node.names <- as.character(igraph::as_ids(vertices))
  } else {
    node.names <- as.character(vertices)
  }
  
  graph_vertex_names <- igraph::as_ids(igraph::V(graph))
  
  if(verbose) {
    message("Precomputing neighborhood sizes for all nodes")
  }
  
  all_neighborhood_size <- igraph::ego_size(
    graph = graph,
    nodes = igraph::V(graph),
    mode = mode,
    order = 1
  ) - 1
  
  names(all_neighborhood_size) <- graph_vertex_names
  
  if(verbose) {
    message("Getting the first neighbors of selected nodes")
  }
  
  # Faster igraph::neighbors() semantics.
  node.neighbors <- lapply(
    node.names,
    function(i) {
      as.character(
        igraph::as_ids(
          igraph::neighbors(
            graph = graph,
            v = i,
            mode = mode
          )
        )
      )
    }
  )
  
  names(node.neighbors) <- node.names
  
  if(verbose) {
    message("Calculating the neighborhood connectivity of nodes")
  }
  
  temp.nc <- vapply(
    node.neighbors,
    function(nbrs) {
      
      if(length(nbrs) == 0) {
        return(0)
      }
      
      mean(all_neighborhood_size[nbrs], na.rm = TRUE)
    },
    numeric(1)
  )
  
  temp.nc[!is.finite(temp.nc)] <- 0
  names(temp.nc) <- node.names
  
  return(temp.nc)
}