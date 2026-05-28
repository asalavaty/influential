#=============================================================================
#
#    Calculation of the conditional probability of deviation from
#                  means in opposite directions
#
#=============================================================================

#' Conditional probability of deviation from means
#'
#' This function calculates the conditional probability of deviation of two
#' centrality measures (or any two other continuous variables) from their corresponding
#' means in opposite directions.
#' @param data A data frame containing the values of two continuous variables and the name of
#' observations (nodes).
#' @param nodes.colname The character format (quoted) name of the column containing
#' the name of observations (nodes).
#' @param Desired.colname The character format (quoted) name of the column containing
#' the values of the desired variable.
#' @param Condition.colname The character format (quoted) name of the column containing
#' the values of the condition variable.
#' @return A list of two objects including the conditional probability of deviation of two
#' centrality measures (or any two other continuous variables) from their corresponding
#' means in opposite directions based on both the entire network and the split-half random
#' sample of network nodes.
#' @aliases CPA
#' @keywords conditional_probability association_assessment
#' @family centrality association assessment functions
#' @export cond.prob.analysis
#' @examples
#' \dontrun{
#' MyData <- centrality.measures
#' My.conditional.prob <- cond.prob.analysis(data = MyData,
#'                                           nodes.colname = rownames(MyData),
#'                                           Desired.colname = "BC",
#'                                           Condition.colname = "NC")
#'                                           }
cond.prob.analysis <- function(data, nodes.colname, Desired.colname, Condition.colname) {
  
  #filtering the data to find those nodes meeting the conditions
  ncpositive <- data[data[, Condition.colname] >
                       mean(data[, Condition.colname]),]
  ncpositive.bcnegative <- ncpositive[ncpositive[, Desired.colname] <
                                        mean(data[, Desired.colname]),]
  ncpositive.bcnegative.prob <- (nrow(ncpositive.bcnegative)/nrow(ncpositive))*100
  
  ncnegative <- data[data[, Condition.colname] <
                       mean(data[, Condition.colname]),]
  ncnegative.bcpositive <- ncnegative[ncnegative[, Desired.colname] >
                                        mean(data[, Desired.colname]),]
  ncnegative.bcpositive.prob <- (nrow(ncnegative.bcpositive)/nrow(ncnegative))*100
  
  #calculation of conditional probability
  final.cond.prob <- sum(ncpositive.bcnegative.prob, ncnegative.bcpositive.prob)/2
  
  ##Reliability analysis based on split-half random sampling
  
  #split-half random sampling
  
  sample.data <- data[sample(1:nrow(data), size = round(nrow(data)/2), replace = FALSE),]
  
  #filtering the data to find those nodes meeting the conditions
  sample.ncpositive <- sample.data[sample.data[, Condition.colname] >
                                     mean(sample.data[, Condition.colname]),]
  sample.ncpositive.bcnegative <- sample.ncpositive[sample.ncpositive[, Desired.colname] <
                                                      mean(sample.data[, Desired.colname]),]
  sample.ncpositive.bcnegative.prob <- (nrow(sample.ncpositive.bcnegative)/nrow(sample.ncpositive))*100
  
  sample.ncnegative <- sample.data[sample.data[, Condition.colname] <
                                     mean(sample.data[, Condition.colname]),]
  sample.ncnegative.bcpositive <- sample.ncnegative[sample.ncnegative[, Desired.colname] >
                                                      mean(sample.data[, Desired.colname]),]
  sample.ncnegative.bcpositive.prob <- (nrow(sample.ncnegative.bcpositive)/nrow(sample.ncnegative))*100
  
  #calculation of conditional probability
  sample.final.cond.prob <- sum(sample.ncpositive.bcnegative.prob, sample.ncnegative.bcpositive.prob)/2
  
  results <- list(ConditionalProbability = final.cond.prob,
                  ConditionalProbability_split.half.sample = sample.final.cond.prob)
  
  return(results)
}
