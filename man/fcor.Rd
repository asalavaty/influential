% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{fcor}
\alias{fcor}
\title{Fast correlation and mutual rank analysis}
\usage{
fcor(
  data,
  use = "everything",
  method = "spearman",
  mutualRank = TRUE,
  mutualRank_mode = "unsigned",
  pvalue = FALSE,
  adjust = "BH",
  flat = TRUE
)
}
\arguments{
\item{data}{a numeric dataframe/matrix (features on columns and samples on rows).}

\item{use}{The NA handler, as in R's cov() and cor() functions. Options are "everything", "all.obs", and "complete.obs".}

\item{method}{a character string indicating which correlation coefficient is to be computed. One of "pearson" or "spearman" (default).}

\item{mutualRank}{logical, whether to calculate mutual ranks of correlations or not.}

\item{mutualRank_mode}{a character string indicating whether to rank based on "signed" or "unsigned" (default) correlation values. 
In the "unsigned" mode, only the level of a correlation value is important and not its sign (the function ranks the absolutes of correlations). 
Options are "unsigned", and "signed".}

\item{pvalue}{logical, whether to calculate p-values of correlations or not.}

\item{adjust}{p-value correction method (when pvalue = TRUE), a character string including any of "BH" (default),
"bonferroni", "holm", "hochberg", "hommel", or "none".}

\item{flat}{logical, whether to combine and flatten the result matrices or not.}
}
\value{
Depending on the input data, a dataframe or list including cor (correlation coefficients),
mr (mutual ranks of correlation coefficients), p (p-values of correlation coefficients), and p.adj (adjusted p-values).
}
\description{
This function calculates Pearson/Spearman correlations between all pairs of features in a matrix/dataframe much faster than the base R cor function.
It is also possible to simultaneously calculate mutual rank (MR) of correlations as well as their p-values and adjusted p-values.
Additionally, this function can automatically combine and flatten the result matrices. Selecting correlated features using an MR-based threshold
rather than based on their correlation coefficients or an arbitrary p-value is more efficient and accurate in inferring
functional associations in systems, for example in gene regulatory networks.
}
\examples{
\dontrun{
set.seed(1234)
data <- datasets::attitude
cor <- fcor(data = data)
}
}
\seealso{
\code{\link[coop]{pcor}}, \code{\link[stats]{p.adjust}},
and \code{\link[influential]{graph_from_data_frame}}
}
\keyword{fcor}
