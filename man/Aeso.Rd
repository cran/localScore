% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Aeso}
\alias{Aeso}
\title{Congenital oesophageal atresia data}
\format{
A matrix of 2191 lines and 2 columns. Each line is a day on the first column, and associated to a case (0/1) on the second column.
}
\source{
Dolk H. Secular pattern of congenital oesophageal atresia--George Knox, 1959. J Epidemiol Community Health. 1997;51(2):114-115. <doi:10.1136/jech.51.2.114>
}
\usage{
data(Aeso)
}
\description{
The data consists of individual dates of birth over n=35 cases of the birth defects oesophageal and tracheo-oesophagean fistula
observed in a hospital in Birmingham, U.K., over 2191 days from 1950 through 1955, with Day one set as 1 January 1950
}
\examples{
data(Aeso)
head(Aeso)
p <- sum(Aeso[,2]) / dim(Aeso)[1]
print(p)
}
\keyword{datasets}
