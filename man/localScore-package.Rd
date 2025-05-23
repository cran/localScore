% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localScore-package.R
\docType{package}
\name{localScore-package}
\alias{localScore}
\alias{localScore-package}
\title{localScore: Package for Sequence Analysis by Local Score}
\description{
Provides functionalities for:
\itemize{
\item calculating the local score
\item calculating statistical relevance (p-value) to find a local Score in a sequence of given distribution
}

Given a sequence of numerical score \eqn{X_1,\dots,X_n}, the local score is defined : \eqn{H_n = \max_{1 \leq i \leq j \leq n} \sum_{l = i}^j X_l}
This package find the value local score and the associated sub-sequence, and also sub-optimal local scores and segments. The complexity is linear with \eqn{n}. It can be viewed as a generalization of a sliding window method, considering all windows size. In order to be pertinent, the expectation of the scores \eqn{X_i} should be negative. Most of the methods concerning statistical relevance implemented in this package only applied on integer scores.
}
\details{
Please refer to the vignette of this package or the manual for details on how to use this package.
}
\references{
\itemize{
\item An Improved Approximation For Assessing The Statistical Significance of molecular Sequence Features, Mercier and al 2003
\item Exact distribution for the local score of one i.i.d. random sequence, Sabine Mercier and JJ Daudin, 2001
\item Limit Distributions of Maximal Segmental Score among Markov-Dependent Partial Sums, Karlin and Dembo 1992
\item Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes, Karlin and al 1990
\item Detection de courts segments inverses dans les genomes: methodes et applications, David Robelin 2005
\item A Linear Time Algorithm for Finding All Maximal Scoring Subsequences, Constable and Bates 1985
}
}
\author{
\strong{Maintainer}: David Robelin \email{david.robelin@inrae.fr}

Authors:
\itemize{
  \item Sebastian Simon
  \item Chris Verschelden
  \item Charly Marty
  \item Sabine Mercier \email{sabine.mercier@univ-tlse2.fr}
  \item Sebastien Dejean \email{sebastien.dejean@math.univ-toulouse.fr}
}

Other contributors:
\itemize{
  \item The authors of Eigen the library for the included version of Eigen [copyright holder]
}

}
\keyword{Package}
\keyword{analysis}
\keyword{by}
\keyword{for}
\keyword{internal}
\keyword{local}
\keyword{score}
\keyword{sequence}
