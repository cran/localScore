% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dico}
\alias{dico}
\title{Dictionnaire}
\format{
A score function for the 20 amino acid
}
\source{
Kyte & Doolittle (1982) J. Mol. Biol. 157, 105-132
}
\usage{
dico
}
\description{
Provides the score related to each base of the sequences.
}
\examples{
data(dico)
dico
data(MidSeq)
MidSeq
MidSeqScore=CharSequence2ScoreSequence(MidSeq,dico)
MidSeqScore[1:30]
localScoreC(MidSeqScore)$localScore
}
\keyword{datasets}
