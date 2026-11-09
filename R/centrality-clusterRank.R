#=============================================================================
#
#    ClusterRank
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
                        directed = FALSE, loops = TRUE,
                        ncores = "default", verbose = FALSE) {
  
  # Get selected vertex indices
  if(inherits(vids, "igraph.vs")) {
    vertices.index <- stats::na.omit(match(vids, igraph::V(graph)))
  } else {
    vertices.index <- stats::na.omit(match(vids, igraph::as_ids(igraph::V(graph))))
  }
  
  if(length(vertices.index) == 0) {
    return(numeric(0))
  }
  
  graph_names <- igraph::as_ids(igraph::V(graph))
  
  # -------------------------------------------------------------------------
  # Directed graph
  # -------------------------------------------------------------------------
  if(directed) {
    
    n_cores <- ifelse(ncores == "default", parallel::detectCores() - 1, ncores)
    n_cores <- max(1L, as.integer(n_cores))
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    vertex.transitivity <- vector(mode = "numeric")
    
    if(verbose) {
      message("Calculating the transitivity of nodes")
    }
    
    cl.Rank.mode <- "out"
    
    vertex.transitivity <- foreach::foreach(
      i = igraph::V(graph),
      .combine = c,
      .packages = "igraph"
    ) %dopar% {
      vertex.neighborhood <- igraph::neighborhood(
        graph = graph,
        order = 1,
        nodes = i,
        mode = cl.Rank.mode
      )[[1]][-1]
      
      if(length(vertex.neighborhood) < 2) {
        NaN
      } else {
        indc.subgraph <- igraph::induced.subgraph(graph, vertex.neighborhood)
        igraph::ecount(indc.subgraph) /
          (igraph::vcount(indc.subgraph) * (igraph::vcount(indc.subgraph) - 1))
      }
    }
    
    if(verbose) {
      message("Getting the neighborhood of selected nodes and calculating the ClusterRank")
    }
    
    cl.Rank <- foreach::foreach(
      i = igraph::V(graph)[vertices.index],
      .combine = c,
      .multicombine = TRUE,
      .packages = "igraph"
    ) %dopar% {
      if(is.nan(vertex.transitivity[i])) {
        NaN
      } else {
        selected.v.neighborhood <- igraph::neighborhood(
          graph = graph,
          order = 1,
          nodes = i,
          mode = cl.Rank.mode
        )[[1]][-1]
        
        temp.cl.Rank <- 0
        
        for(j in selected.v.neighborhood) {
          temp.cl.Rank <- temp.cl.Rank +
            igraph::degree(
              graph = graph,
              v = j,
              mode = cl.Rank.mode,
              loops = loops
            ) + 1
        }
        
        temp.cl.Rank * vertex.transitivity[i]
      }
    }
    
    if(igraph::is_named(graph)) {
      names(cl.Rank) <- graph_names[vertices.index]
    }
    
    return(cl.Rank)
  }
  
  # -------------------------------------------------------------------------
  # Undirected graph
  # -------------------------------------------------------------------------
  
  cl.Rank.mode <- "all"
  
  if(verbose) {
    message("Calculating the transitivity of nodes")
  }
  
  vertex.transitivity <- igraph::transitivity(
    graph = graph,
    type = "local"
  )
  
  if(verbose) {
    message("Precomputing node degrees")
  }
  
  degree_all <- igraph::degree(
    graph = graph,
    v = igraph::V(graph),
    mode = cl.Rank.mode,
    loops = loops
  )
  
  if(verbose) {
    message("Getting the neighborhood of selected nodes and calculating the ClusterRank")
  }
  
  selected_neighborhoods <- igraph::neighborhood(
    graph = graph,
    order = 1,
    nodes = igraph::V(graph)[vertices.index],
    mode = cl.Rank.mode
  )
  
  cl.Rank <- vapply(
    seq_along(selected_neighborhoods),
    function(k) {
      
      v_index <- vertices.index[k]
      
      if(is.nan(vertex.transitivity[v_index])) {
        return(NaN)
      }
      
      selected.v.neighborhood <- selected_neighborhoods[[k]][-1]
      
      if(length(selected.v.neighborhood) == 0) {
        temp.cl.Rank <- 0
      } else {
        neigh_idx <- as.integer(selected.v.neighborhood)
        temp.cl.Rank <- sum(degree_all[neigh_idx] + 1)
      }
      
      temp.cl.Rank * vertex.transitivity[v_index]
    },
    numeric(1)
  )
  
  if(igraph::is_named(graph)) {
    names(cl.Rank) <- graph_names[vertices.index]
  }
  
  return(cl.Rank)
}
