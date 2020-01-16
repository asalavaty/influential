#=============================================================================
#
#    Code chunk 1: Documenting the coexpression.data data frame
#
#=============================================================================

#' Co-expression dataset
#'
#' A co-expression dataset of lncRNAs and mRNAs in lung adenocarcinoma
#'
#' @format A data frame with 2410 rows and 2 variables:
#' \describe{
#'   \item{lncRNA}{lncRNA symbol}
#'   \item{Coexpressed.Gene}{Co-expressed gene symbol}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/31211495/}
"coexpression.data"

#=============================================================================
#
#    Code chunk 2: Documenting the centrality.measures data frame
#
#=============================================================================

#' Centrality measures dataset
#'
#' The centrality measures of a co-expression network of lncRNAs and
#' mRNAs in lung adenocarcinoma
#'
#' @format A data frame with 794 rows and 4 variables:
#' \describe{
#'   \item{BetweennessCentrality}{Betweenness Centrality}
#'   \item{Degree}{Degree Centrality}
#'   \item{name}{Node (gene) name}
#'   \item{NeighborhoodConnectivity}{Neighborhood Connectivity}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/31211495/}
"centrality.measures"

#=============================================================================
#
#    Code chunk 3: Documenting the coexpression.adjacency data frame
#
#=============================================================================

#' Adjacency matrix
#'
#' The adjacency matrix of a co-expression network of lncRNAs and
#' mRNAs in lung adenocarcinoma that was generated using igraph functions
#'
#' @format A data frame with 794 rows and 794 variables:
#' \describe{
#'   \item{lncRNA}{lncRNA symbol}
#'   \item{lncRNA}{lncRNA symbol}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/31211495/}
"coexpression.adjacency"
