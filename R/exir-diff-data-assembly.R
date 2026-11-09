#=============================================================================
#
#    Assembling the differential/regression data in a dataframe
#
#=============================================================================

#' Assembling the differential/regression data
#'
#' This function assembles a dataframe required for running the \emph{\strong{\code{ExIR}}} model. You may provide
#' as many differential/regression data as you wish. Also, the datasets should be filtered
#' beforehand according to your desired thresholds and, consequently, should only include the significant data.
#' Each dataset provided should be a dataframe with one or two columns.
#' The first column should always include differential/regression values
#' and the second one (if provided) the significance values. Please also note that the significance (adjusted P-value)
#' column is mandatory for differential datasets.
#' @param ... Desired datasets/dataframes.
#' @return A dataframe including the collective list of features in rows and all of the
#' differential/regression data and their statistical significance in columns with the same
#' order provided by the user.
#' @aliases DDA
#' @keywords diff_data.assembly
#' @seealso \code{\link[influential]{exir}}
#' @export diff_data.assembly
#' @examples
#' \dontrun{
#' my.Diff_data <- diff_data.assembly(Differential_data1,
#'                                    Differential_data2,
#'                                    Regression_data1)
#' }
diff_data.assembly <- function(...) {
  
  #Getting the list of all datasets provided
  datasets <- lapply(list(...), as.data.frame)
  
  #Getting the feature names
  feature.names <- unique(unlist(lapply(X = datasets, FUN = rownames)))
  
  #Creating the Diff_data dataset
  Diff_data <- data.frame(Diff_value1 = rep(0,length(feature.names)),
                          row.names = feature.names)
  
  for (i in 1:length(datasets)) {
    
    feature.names.index <- match(rownames(datasets[[i]]),
                                 rownames(Diff_data))
    
    if(ncol(datasets[[i]]) == 2) {
      
      Diff_data[,paste("Diff_value", i, sep = "")] <- 0
      Diff_data[,paste("Diff_value", i, sep = "")][feature.names.index] <- datasets[[i]][,1]
      
      Diff_data[,paste("Sig_value", i, sep = "")] <- 1
      Diff_data[,paste("Sig_value", i, sep = "")][feature.names.index] <- datasets[[i]][,2]
      
    } else if(ncol(datasets[[i]]) == 1) {
      
      Diff_data[,paste("Diff_value", i, sep = "")] <- 0
      Diff_data[,paste("Diff_value", i, sep = "")][feature.names.index] <- datasets[[i]][,1]
      
    }
  }
  
  return(Diff_data)
}
