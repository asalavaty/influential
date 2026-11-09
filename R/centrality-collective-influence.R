#=============================================================================
#
#    Collective Influence (CI)
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
collective.influence <- function(graph, vertices = V(graph),
                                 mode = "all", d = 3, verbose = FALSE) {
  
  vertices.index <- .get_vertex_indices(graph, vertices)
  graph_names <- igraph::as_ids(igraph::V(graph))
  
  if(length(vertices.index) == 0) {
    return(numeric(0))
  }
  
  if(verbose) {
    message("Calculating collective influence")
  }
  
  # Reduced degrees for all nodes
  reduced_degrees_all <- igraph::degree(
    graph = graph,
    v = igraph::V(graph),
    mode = mode
  ) - 1
  
  names(reduced_degrees_all) <- graph_names
  
  # Nodes exactly at distance d
  nodes.at.distance <- igraph::neighborhood(
    graph = graph,
    nodes = graph_names[vertices.index],
    mode = mode,
    order = d,
    mindist = d
  )
  
  distance_reduced_degree_sum <- vapply(
    nodes.at.distance,
    function(vs) {
      ids <- igraph::as_ids(vs)
      sum(reduced_degrees_all[ids], na.rm = TRUE)
    },
    numeric(1)
  )
  
  ci <- reduced_degrees_all[graph_names[vertices.index]] *
    distance_reduced_degree_sum
  
  names(ci) <- graph_names[vertices.index]
  ci
}
