% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rMethods.R
\name{sequences2transmatrix}
\alias{sequences2transmatrix}
\title{Transition matrix from sequence(s)}
\usage{
sequences2transmatrix(sequences)
}
\arguments{
\item{sequences}{Sequences to be analyzed, can be a list of vectors, or a vector}
}
\value{
A list object containing
\item{transition_matrix}{Transition Matrix with row names and column names are the associated score/state}
\item{score_value}{a vector containing the score/state, ordered in the same way as the matrix columns and rows}
}
\description{
Calculates the transition matrix by counting occurrences of tuples in given vector list
}
\details{
In output, score_value is coerced as integer if possible. Else, it is a character vector containing the states
of the Markov chain.
}
\examples{
myseq <- sample(LETTERS[1:2], size=20,replace=TRUE)
sequences2transmatrix(myseq)
}
