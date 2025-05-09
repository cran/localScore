% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rMethods.R
\name{lindley}
\alias{lindley}
\title{Lindley process}
\usage{
lindley(sequence)
}
\arguments{
\item{sequence}{numeric sequence of a Lindley process, eg service time per customer}
}
\value{
a vector with the Lindley process steps
}
\description{
Creates a sequence of a Lindley process, also called CUSUM process, on a given sequence.  For a sequence (X_k)k, the Lindley process is defined as follows: W_0:=0 and W_(k+1)=max(0,W_k+X_(k+1)). It defines positive excursions above 0.
}
\examples{
MySeq <- c(1,2,3,-4,1,-3,-1,2,3,-4,1)
lindley(MySeq)
plot(1:length(MySeq),lindley(MySeq),type='b')
}
