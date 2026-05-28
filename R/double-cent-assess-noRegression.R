#=============================================================================
#
#    Assessment of innate features and associations of two network centrality
#         measures, without considering dependent and independent ones
#
#=============================================================================


#' Assessment of innate features and associations of two network centrality measures
#'
#' This function assesses innate features and the association of two centrality measures
#' (or any two other continuous variables) from the aspect of distribution mode, dependence,
#' linearity, partial-moments based correlation, and conditional probability of
#' deviating from corresponding means in opposite direction (centrality2 is used as the condition variable).
#' This function doesn't consider which variable is dependent and which one is
#' independent and no regression analysis is done.
#' Also, the correlation between two variables is assessed via non-linear non-parametric statistics (NNS).
#' For the conditional probability assessment, the centrality2 variable is considered
#' as the condition variable.
#' @param data A data frame containing the values of two continuous variables and the name of
#' observations (nodes).
#' @param nodes.colname The character format (quoted) name of the column containing
#' the name of observations (nodes).
#' @param centrality1.colname The character format (quoted) name of the column containing
#' the values of the Centrality_1 variable.
#' @param centrality2.colname The character format (quoted) name of the column containing
#' the values of the Centrality_2 variable.
#' @return A list of nine objects including:
#'
#' - Summary of the basic statistics of two centrality measures (or any two other continuous variables).
#'
#' - The results of normality assessment of two variable (p-value > 0.05 imply that the variable is normally distributed).
#'
#' - Description of the normality assessment of the centrality1 (first variable).
#'
#' - Description of the normality assessment of the centrality2 (second variable).
#'
#' - The Hoeffding's D Statistic of dependence (ranging from -0.5 to 1).
#'
#' - Description of the dependence significance.
#'
#' - Correlation between variables based on the NNS method.
#'
#' - The last two objects are the conditional probability of deviation of two
#' centrality measures from their corresponding means in opposite directions based
#' on both the entire network and the split-half random sample of network nodes.
#' @aliases DCANR
#' @keywords association_assessment dependence_assessment
#' @family centrality association assessment functions
#' @seealso \code{\link[nortest]{ad.test}} for Anderson-Darling test for normality,
#' \code{\link[Hmisc]{hoeffd}} for Matrix of Hoeffding's D Statistics, and
#' \code{\link[NNS]{NNS.dep}} for NNS Dependence
#' @export double.cent.assess.noRegression
#' @examples
#' \dontrun{
#' MyData <- centrality.measures
#' My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,
#'                                             nodes.colname = rownames(MyData),
#'                                             centrality1.colname = "BC",
#'                                             centrality2.colname = "NC")
#' }
double.cent.assess.noRegression <- function(data, nodes.colname,
                                            centrality1.colname,
                                            centrality2.colname) {
  
  if("parallel" %in% (.packages())) {
    detach("package:parallel", unload = TRUE)
    base::attachNamespace("parallel")
    parallel::detectCores(logical = TRUE)
  } else {
    base::attachNamespace("parallel")
    parallel::detectCores(logical = TRUE)
  }
  
  #Checking the availability of required packages
  
  if (nrow(data) >= 5000) { if(!requireNamespace(c("nortest", "Hmisc", "NNS"), quietly = TRUE)) {
    stop("The packages \"nortest\", \"Hmisc\" and \"NNS\" are required for this function to work.
    Please install the required packages before using this function.

  You can install the packages via one of the following options:

         install.packages(\"Package Name\")

         Or

         install.packages(\"BiocManager\")
         BiocManager::install(\"Package Name\")",
         call. = FALSE)
  }
  }
  
  if(!requireNamespace(c("Hmisc", "mgcv", "NNS"), quietly = TRUE)) {
    stop("The packages \"Hmisc\", \"mgcv\" and \"NNS\" are required for this function to work.
    Please install the required packages before using this function.

  You can install the packages via one of the following options:

         install.packages(\"Package Name\")

         Or

         install.packages(\"BiocManager\")
         BiocManager::install(\"Package Name\")",
         call. = FALSE)
  }
  
  #checking the normality of data
  summary.stat <- apply(data[, c(centrality1.colname, centrality2.colname)], 2, summary)
  
  if(length(unique(data[,centrality1.colname])) < 3 &
     length(unique(data[,centrality2.colname])) >= 3) {
    if(nrow(data) < 5000) {
      normality <- data.frame(p.value = c(NA, stats::shapiro.test(data[, centrality2.colname])$p.value))
    } else if (nrow(data) >= 5000) {
      normality <- data.frame(p.value = c(NA, nortest::ad.test(data[, centrality2.colname])$p.value))
    }
  } else if (length(unique(data[,centrality1.colname])) >= 3 &
             length(unique(data[,centrality2.colname])) < 3) {
    if(nrow(data) < 5000) {
      normality <- data.frame(p.value = c(stats::shapiro.test(data[, centrality1.colname])$p.value, NA))
    } else if (nrow(data) >= 5000) {
      normality <- data.frame(p.value = c(nortest::ad.test(data[, centrality1.colname])$p.value, NA))
    }
  } else if(length(unique(data[,centrality1.colname])) < 3 &
            length(unique(data[,centrality2.colname])) < 3) {
    normality <- data.frame(p.value = c(NA, NA))
  } else {
    if(nrow(data) < 5000) {
      normality <- apply(data[, c(centrality1.colname, centrality2.colname)], 2, stats::shapiro.test)
    } else if(nrow(data) >= 5000) {
      normality <- apply(data[, c(centrality1.colname, centrality2.colname)], 2, nortest::ad.test)
    }
    normality <- as.data.frame(sapply(normality, function(m) m[]$p.value))
    colnames(normality) <- "p.value"
  }
  
  if(is.na(normality[1,1])) {dependent.normality <- NA} else if(normality[1,1] < 0.05) {
    dependent.normality <- "Non-normally distributed"
  } else {dependent.normality <- "Normally distributed"}
  
  if(is.na(normality[2,1])) {independent.normality <- NA} else if(normality[2,1] < 0.05) {
    independent.normality <- "Non-normally distributed"
  } else {independent.normality <- "Normally distributed"}
  
  #calculation of Hoeffding’s D Statistics (Hoeffding Dependence Coefficient)
  hoeffd <- data.frame(D_statistic = as.data.frame(Hmisc::hoeffd(x = data[, centrality2.colname],
                                                                 y = data[, centrality1.colname])[1])[1,2],
                       P_value = as.data.frame(Hmisc::hoeffd(x = data[, centrality2.colname],
                                                             y = data[, centrality1.colname])[3])[1,2], row.names = "Results")
  
  if(hoeffd[1,2] < 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Significantly dependent",
                                          row.names = "Results")
  } else if(hoeffd[1,2] >= 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Not significantly dependent",
                                          row.names = "Results")
  } else if(hoeffd[1,2] < 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Significantly dependent",
                                          row.names = "Results")
  } else if (hoeffd[1,2] >= 0.05) {
    dependence.significance <- data.frame(Hoeffding = "Not significantly dependent",
                                          row.names = "Results")
  }
  
  
  ##assessment of descriptive non-linear non-parametric correlation/dependence between
  #dependent and independent variables
  
  nl.cor.dep <- NNS::NNS.dep(x = data[, centrality2.colname], y = data[, centrality1.colname],
                             print.map = FALSE, order = 2)
  
  nl.cor.dep <- data.frame(Correlation = unlist(nl.cor.dep)[1],
                           Dependence = unlist(nl.cor.dep)[2], row.names = "Results")
  
  ##assessment of conditional probability of deviation of BC and NC from their means in opposite directions
  #filtering the data to find those nodes meeting the conditions
  ncpositive <- data[data[, centrality2.colname] >
                       mean(data[, centrality2.colname]),]
  ncpositive.bcnegative <- ncpositive[ncpositive[, centrality1.colname] <
                                        mean(data[, centrality1.colname]),]
  ncpositive.bcnegative.prob <- (nrow(ncpositive.bcnegative)/nrow(ncpositive))*100
  
  ncnegative <- data[data[, centrality2.colname] <
                       mean(data[, centrality2.colname]),]
  ncnegative.bcpositive <- ncnegative[ncnegative[, centrality1.colname] >
                                        mean(data[, centrality1.colname]),]
  ncnegative.bcpositive.prob <- (nrow(ncnegative.bcpositive)/nrow(ncnegative))*100
  
  #calculation of conditional probability
  final.cond.prob <- sum(ncpositive.bcnegative.prob, ncnegative.bcpositive.prob)/2
  
  ##Reliability analysis based on split-half random sampling
  
  #split-half random sampling
  
  sample.data <- data[sample(1:nrow(data), size = round(nrow(data)/2), replace = FALSE),]
  
  #filtering the data to find those nodes meeting the conditions
  sample.ncpositive <- sample.data[sample.data[, centrality2.colname] >
                                     mean(sample.data[, centrality2.colname]),]
  sample.ncpositive.bcnegative <- sample.ncpositive[sample.ncpositive[, centrality1.colname] <
                                                      mean(sample.data[, centrality1.colname]),]
  sample.ncpositive.bcnegative.prob <- (nrow(sample.ncpositive.bcnegative)/nrow(sample.ncpositive))*100
  
  sample.ncnegative <- sample.data[sample.data[, centrality2.colname] <
                                     mean(sample.data[, centrality2.colname]),]
  sample.ncnegative.bcpositive <- sample.ncnegative[sample.ncnegative[, centrality1.colname] >
                                                      mean(sample.data[, centrality1.colname]),]
  sample.ncnegative.bcpositive.prob <- (nrow(sample.ncnegative.bcpositive)/nrow(sample.ncnegative))*100
  
  #calculation of conditional probability
  sample.final.cond.prob <- sum(sample.ncpositive.bcnegative.prob, sample.ncnegative.bcpositive.prob)/2
  
  results <- list(Summary_statistics = summary.stat,
                  Normality_results = normality,
                  Centrality1_Normality = dependent.normality,
                  Centrality2_Normality = independent.normality,
                  HoeffdingD_Statistic = hoeffd,
                  Dependence_Significance = dependence.significance,
                  NNS_dep_results = nl.cor.dep,
                  ConditionalProbability = final.cond.prob,
                  ConditionalProbability_split.half.sample = sample.final.cond.prob)
  
  return(results)
  
}
