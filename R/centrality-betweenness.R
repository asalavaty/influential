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
#'                                     directed = FALSE,
#'                                     normalized = FALSE)
#' }
betweenness <- function(graph,
                        v = V(graph),
                        directed = TRUE,
                        weights = NULL,
                        normalized = FALSE,
                        ...) {
  igraph::betweenness(
    graph = graph,
    v = v,
    directed = directed,
    weights = weights,
    normalized = normalized,
    ...
  )
}