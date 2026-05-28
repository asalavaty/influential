#=============================================================================
#
#    IVI
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
                mode = "all", loops = TRUE, d = 3, scale = "range",
                ncores = "default", verbose = FALSE) {
  
  if(verbose) {
    message("Calculating IVI centrality components")
  }
  
  # If the graph has an edge weight and weights is NULL, igraph::betweenness()
  # will use the edge weight by default. Preserve that behaviour unless the user
  # explicitly provides weights.
  
  if(verbose) {
    message("Calculating degree centrality")
  }
  
  DC <- igraph::degree(
    graph = graph,
    v = vertices,
    mode = mode,
    loops = loops
  )
  
  if(verbose) {
    message("Calculating ClusterRank")
  }
  
  CR <- clusterRank(
    graph = graph,
    vids = vertices,
    directed = directed,
    loops = loops,
    ncores = ncores,
    verbose = verbose
  )
  
  if(verbose) {
    message("Calculating local H-index")
  }
  
  LH_index <- lh_index(
    graph = graph,
    vertices = vertices,
    mode = mode,
    ncores = ncores,
    verbose = verbose
  )
  
  if(verbose) {
    message("Calculating neighborhood connectivity")
  }
  
  NC <- neighborhood.connectivity(
    graph = graph,
    vertices = vertices,
    mode = mode,
    verbose = verbose
  )
  
  if(verbose) {
    message("Calculating betweenness centrality")
  }
  
  BC <- igraph::betweenness(
    graph = graph,
    v = vertices,
    directed = directed,
    weights = weights
  )
  
  if(verbose) {
    message("Calculating collective influence")
  }
  
  CI <- collective.influence(
    graph = graph,
    vertices = vertices,
    mode = mode,
    d = d,
    verbose = verbose
  )
  
  # -----------------------------------------------------------------------
  # Centrality vectors
  # -----------------------------------------------------------------------
  
  temp.DC <- DC
  temp.CR <- CR
  temp.LH_index <- LH_index
  temp.NC <- unlist(NC)
  temp.BC <- BC
  temp.CI <- CI
  
  temp.DC[!is.finite(temp.DC)] <- 0
  temp.CR[!is.finite(temp.CR)] <- 0
  temp.LH_index[!is.finite(temp.LH_index)] <- 0
  temp.NC[!is.finite(temp.NC)] <- 0
  temp.BC[!is.finite(temp.BC)] <- 0
  temp.CI[!is.finite(temp.CI)] <- 0
  
  # -----------------------------------------------------------------------
  # 1-100 normalization of centrality measures
  # -----------------------------------------------------------------------
  
  if(verbose) {
    message("Normalizing centrality measures")
  }
  
  temp.DC <- .range_normalize_1_100(temp.DC)
  temp.CR <- .range_normalize_1_100(temp.CR)
  temp.LH_index <- .range_normalize_1_100(temp.LH_index)
  temp.NC <- .range_normalize_1_100(temp.NC)
  temp.BC <- .range_normalize_1_100(temp.BC)
  temp.CI <- .range_normalize_1_100(temp.CI)
  
  # -----------------------------------------------------------------------
  # IVI calculation
  # -----------------------------------------------------------------------
  
  if(verbose) {
    message("Calculating IVI")
  }
  
  spreading.rank <- ((temp.NC + temp.CR) * (temp.BC + temp.CI))
  spreading.rank[!is.finite(spreading.rank) | spreading.rank == 0] <- 1
  
  hubness.rank <- (temp.DC + temp.LH_index)
  hubness.rank[!is.finite(hubness.rank) | hubness.rank == 0] <- 1
  
  temp.ivi <- hubness.rank * spreading.rank
  
  if(scale == "range") {
    
    if(verbose) {
      message("Normalizing IVI to 1-100 range")
    }
    
    if(length(unique(temp.ivi)) > 1) {
      temp.ivi <- 1 + (((temp.ivi - min(temp.ivi)) * (100 - 1)) /
                         (max(temp.ivi) - min(temp.ivi)))
    }
    
  } else if(scale == "z-scale") {
    
    if(verbose) {
      message("Z-score standardization of IVI")
    }
    
    temp.ivi <- base::scale(temp.ivi)
    temp.ivi.names <- rownames(temp.ivi)
    temp.ivi <- c(temp.ivi)
    names(temp.ivi) <- temp.ivi.names
  }
  
  temp.ivi
}
