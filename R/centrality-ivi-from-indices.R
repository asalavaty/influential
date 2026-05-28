#=============================================================================
#
#    IVI from centrality measures
#
#=============================================================================

#' Integrated Value of Influence (IVI)
#'
#' This function calculates the IVI of the desired nodes from previously calculated centrality
#' measures. This function is not dependent to other packages and the required centrality
#' measures, namely degree centrality, ClusterRank, betweenness centrality, Collective Influence,
#' local H-index, and neighborhood connectivity could have been calculated by any means beforehand.
#' A shiny app has also been developed for the calculation of IVI as well as IVI-based network
#' visualization, which is accessible using the `influential::runShinyApp("IVI")` command.
#' You can also access the shiny app online at https://influential.erc.monash.edu/.
#' @param DC A vector containing the values of degree centrality of the desired vertices.
#' @param CR A vector containing the values of ClusterRank of the desired vertices.
#' @param LH_index A vector containing the values of local H-index of the desired vertices.
#' @param NC A vector containing the values of neighborhood connectivity of the desired vertices.
#' @param BC A vector containing the values of betweenness centrality of the desired vertices.
#' @param CI A vector containing the values of Collective Influence of the desired vertices.
#' @param scale Character string; the method used for scaling/normalizing the results. Options include 'range' (normalization within a 1-100 range), 
#' 'z-scale' (standardization using the z-score), and 'none' (no data scaling). The default selection is 'range'. Opting for the 'range' method is 
#' suitable when exploring a single network, allowing you to observe the complete spectrum and distribution of node influences. In this case, there is 
#' no intention to establish a specific threshold for the outcomes. However, it is possible to identify and present the top influential nodes 
#' based on their rankings. Conversely, the 'z-scale' option proves advantageous if the aim is to compare node influences across multiple networks or 
#' if there is a desire to establish a threshold (usually z-score > 1.645) for generating a list of the most influential nodes without manual intervention.
#' @param verbose Logical; whether the accomplishment of different stages of the algorithm should be printed (default is FALSE).
#' @return A numeric vector with the IVI values based on the provided centrality measures.
#' @aliases IVI.FI
#' @keywords ivi.from.indices
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{cent_network.vis}}
#' @export ivi.from.indices
#' @examples
#' \dontrun{
#' MyData <- centrality.measures
#' My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,
#'                                     CR = centrality.measures$CR,
#'                                     NC = centrality.measures$NC,
#'                                     LH_index = centrality.measures$LH_index,
#'                                     BC = centrality.measures$BC,
#'                                     CI = centrality.measures$CI)
#'                                     }
ivi.from.indices <- function(DC, CR, LH_index, NC, BC, CI, scale = "range", verbose = FALSE) {
  
  #Generating temporary measures
  
  temp.DC <- DC
  temp.CR <- CR
  temp.LH_index <- LH_index
  temp.NC <- NC
  temp.BC <- BC
  temp.CI <- CI
  
  #Removing the NAN and NA values
  
  temp.DC[c(which(is.nan(temp.DC)), which(is.na(temp.DC)))] <- 0
  temp.CR[c(which(is.nan(temp.CR)), which(is.na(temp.CR)))] <- 0
  temp.LH_index[c(which(is.nan(temp.LH_index)), which(is.na(temp.LH_index)))] <- 0
  temp.NC[c(which(is.nan(temp.NC)), which(is.na(temp.NC)))] <- 0
  temp.BC[c(which(is.nan(temp.BC)), which(is.na(temp.BC)))] <- 0
  temp.CI[c(which(is.nan(temp.CI)), which(is.na(temp.CI)))] <- 0
  
  if(verbose) {
    cat("1-100 normalization of centrality measures\n")
  }
  
  #1-100 normalization of centrality measures
  
  if(length(temp.DC) > 1 & any(temp.DC > 0)) {
    temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
  }
  
  if(length(temp.CR) > 1 & any(temp.CR > 0)) {
    temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
  }
  
  if(length(temp.LH_index) > 1 & any(temp.LH_index > 0)) {
    temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
  }
  
  if(length(temp.NC) > 1 & any(temp.NC > 0)) {
    temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
  }
  
  if(length(temp.BC) > 1 & any(temp.BC > 0)) {
    temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
  }
  
  if(length(temp.CI) > 1 & any(temp.CI > 0)) {
    temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
  }
  
  #Calculation of IVI
  
  if(verbose) {
    cat("Calculating the Spreading Rank\n")
  }
  
  spreading.rank <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))
  
  suppressWarnings(
    if(any(stats::na.omit(spreading.rank) == 0 | is.na(spreading.rank))) {
      spreading.rank[which(spreading.rank == 0 | is.na(spreading.rank))] <- 1
    }
  )
  
  if(verbose) {
    cat("Calculating the Hubness Rank\n")
  }
  
  hubness.rank <- (temp.DC+temp.LH_index)
  
  suppressWarnings(
    if(any(stats::na.omit(hubness.rank) == 0 | is.na(hubness.rank))) {
      hubness.rank[which(hubness.rank == 0 | is.na(hubness.rank))] <- 1
    }
  )
  
  
  if(verbose) {
    cat("Calculating the IVI\n")
  }
  
  temp.ivi <- (hubness.rank)*(spreading.rank)
  
  #1-100 normalization of IVI
  
  if(scale == "range") {
    
    if(verbose) {
      cat("1-100 normalization of IVI\n")
    }
    
    if(length(unique(temp.ivi)) > 1) {
      temp.ivi <- 1+(((temp.ivi-min(temp.ivi))*(100-1))/(max(temp.ivi)-min(temp.ivi)))
    }
  } else if(scale == 'z-scale') {
    if(verbose) {
      cat("Z-score standardization of IVI\n")
    }
    
    temp.ivi <- base::scale(temp.ivi)
    temp.ivi.names <- rownames(temp.ivi)
    temp.ivi <- c(temp.ivi)
    names(temp.ivi) <- temp.ivi.names
  }
  
  return(temp.ivi)
}
