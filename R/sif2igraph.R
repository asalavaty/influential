#=============================================================================
#
#    SIF to igraph
#
#=============================================================================

#' SIF to igraph
#'
#' This function imports and converts a SIF file from your local hard drive, cloud space,
#' or internet into a graph with an igraph class, which can then be used for the identification
#' of most influential nodes via the ivi function, for instance.
#' @param Path A string or character vector indicating the path to the desired SIF file. The SIF file
#' could be on your local hard drive, cloud space, or on the internet.
#' @param directed Logical scalar, whether or not to create a directed graph.
#' @return An igraph graph object.
#' @keywords SIF.to.igraph
#' @family network_reconstruction functions
#' @export sif2igraph
#' @examples
#' \dontrun{
#' MyGraph <- sif2igraph(Path = "/Users/User1/Desktop/mygraph.sif", directed=FALSE)
#' }
sif2igraph <- function(Path, directed=FALSE) {
  
  graph_from_data_frame(d = utils::read.delim(Path, header = FALSE)[,c(1,3)],directed = directed)
}
