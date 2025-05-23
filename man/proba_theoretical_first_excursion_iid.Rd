% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/excursion.R
\name{proba_theoretical_first_excursion_iid}
\alias{proba_theoretical_first_excursion_iid}
\title{Probability \eqn{P(Q(1)\geq a)} that the height of the first excursion is greater or equal to \code{a} given a i.i.d. model on the letters sequence}
\usage{
proba_theoretical_first_excursion_iid(
  a,
  theta,
  theta_distribution,
  score_function
)
}
\arguments{
\item{a}{score strictly positive}

\item{theta}{vector containing the alphabet used}

\item{theta_distribution}{distribution vector of theta}

\item{score_function}{vector containing the scores of each letters of the alphabet (must be in the same order as theta)}
}
\value{
theoretical probability of reaching a score of a on the first excursion supposing an  i.i.d model on the letters sequence
}
\description{
Mathematical definition of an excursion of the Lindley process is based on the record times of the partial
sum sequence associated to the score sequence (see Karlin and Altschul 1990, Karlin and Dembo 1992) and
define the successive times where the partial sums are strictly decreasing. There must be distinguished
from the visual excursions of the Lindley sequence. The number \code{i} is the number of excursion in sequential order. Detailed definitions are given in the vignette.
}
\details{
Beware that a sequence beginning with a negative score gives a "flat" excursion, with score 0 are considered.
}
\examples{
proba_theoretical_first_excursion_iid(3, c("a","b","c","d"), 
                                      c(a=0.1,b=0.2,c=0.4,d=0.3), c(a=-3,b=-1,c=1,d=2))
}
