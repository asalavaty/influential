% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\docType{package}
\name{influential-package}
\alias{influential}
\alias{influential-package}
\title{Influential package}
\description{
The goal of \emph{\strong{\code{influential}}} is to help identification of the most influential nodes in a network
as well as the classification and ranking of top candidate features.
This package contains functions for the classification and ranking of features,
reconstruction of networks from adjacency matrices and
data frames, analysis of the topology of the network and calculation of centrality measures
as well as a novel and powerful influential node ranking.
The \strong{Experimental data-based Integrative Ranking (ExIR)} is a sophisticated model
for classification and ranking of the top candidate features based on only the experimental data.
The first integrative method, namely the \strong{Integrated Value of Influence (IVI)},
that captures all topological dimensions of the network for
the identification of network most influential nodes is also provided as
a function. Also, neighborhood connectivity, H-index, local H-index, and collective
influence (CI), all of which required centrality measures for the calculation of IVI,
are for the first time provided in an R package. Additionally, a function is provided
for running \strong{SIRIR} model, which is the combination of leave-one-out cross validation
technique and the conventional SIR model, on a network to unsupervisedly rank the true
influence of vertices.Furthermore, some functions have been provided for the
assessment of dependence and correlation of two network centrality measures as well
as the conditional probability of deviation from their corresponding
means in opposite directions.

You may check the latest developmental version of the \emph{influential} package on its
\href{https://github.com/asalavaty/influential}{GitHub repository}

Also, a web-based \href{https://influential.erc.monash.edu/}{Influential Software Package} with a convenient
user-interface (UI) has been developed for the comfort of all users including those without a coding background.
}
\details{
\itemize{
  \item Package: influential
  \item Type: Package
  \item Version: 2.2.9
  \item Date: 01-12-2023
  \item License: GPL-3
}
}
\references{
\itemize{
  \item Fred Viole and David Nawrocki (2013, ISBN:1490523995).
  \item Csardi G, Nepusz T (2006). “The igraph software package for complex network research.”
InterJournal, Complex Systems, 1695. \url{https://igraph.org/}.
}

\strong{Note:} Adopted algorithms and sources are referenced in function document.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/asalavaty/influential}
  \item \url{https://asalavaty.github.io/influential/}
  \item Report bugs at \url{https://github.com/asalavaty/influential/issues}
}

}
\author{
Author: Adrian (Abbas) Salavaty

Advisors: Mirana Ramialison and Peter D. Currie


Maintainer: Adrian (Abbas) Salavaty \email{abbas.salavaty@gmail.com}


You may find more information on my personal website at \href{https://asalavaty.com/}{www.ASalavaty.com}
}
\keyword{internal}
