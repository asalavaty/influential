
#=============================================================================
#
#    Helper functions
#
#=============================================================================

.get_vertex_indices <- function(graph, vertices) {
  
  graph_names <- igraph::as_ids(igraph::V(graph))
  
  if(inherits(vertices, "igraph.vs")) {
    vertex_names <- igraph::as_ids(vertices)
  } else {
    vertex_names <- as.character(vertices)
  }
  
  stats::na.omit(match(vertex_names, graph_names))
}


.get_vertex_names <- function(graph, vertices) {
  
  if(inherits(vertices, "igraph.vs")) {
    as.character(igraph::as_ids(vertices))
  } else {
    as.character(vertices)
  }
}


.make_binary_neighbor_graph <- function(graph) {
  
  # Neighbourhood membership in igraph is based on reachable adjacent vertices,
  # not edge multiplicity. Therefore, this simplified graph is used only for
  # neighbour membership / adjacency operations.
  igraph::simplify(
    graph,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )
}


.range_normalize_1_100 <- function(x) {
  
  x[!is.finite(x)] <- 0
  
  if(length(x) > 1 && any(x > 0) && length(unique(x)) > 1) {
    x <- 1 + (((x - min(x)) * (100 - 1)) / (max(x) - min(x)))
  }
  
  x
}