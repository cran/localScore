#' Short protein sequence
#'
#' A short protein sequence of 31 amino acids corresponding to Q09FU3.fasta query in UniProt Data base.
#'
#' @format A character string with 31 characters "MLTITSYFGFLLAALTITSVLFIGLNKIRLI"
#' @source \url{https://www.uniprot.org/}
#' @examples
#' data(ShortSeq)
#' ShortSeq
#' nchar(ShortSeq)
#' data(dico) 
#' SeqScore=CharSequence2ScoreSequence(ShortSeq,dico)
#' SeqScore
#' localScoreC(SeqScore)$localScore
#' LS=localScoreC(SeqScore)$localScore[1]
#' prob1 = scoreSequences2probabilityVector(list(SeqScore))
#' daudin(localScore = LS, sequence_length = nchar(ShortSeq),
#'                score_probabilities = prob1,
#'                sequence_min = min(SeqScore),
#'                sequence_max = max(SeqScore))
#' score=-5:5
#' prob2=c(0.15,0.15,0.1,0.1,0.0,0.05,0.15,0.05,0.2,0.0,0.05)
#' sum(prob2*score)
#' karlin(localScore = LS, sequence_length = nchar(ShortSeq),
#' score_probabilities = prob2,
#' sequence_min = min(SeqScore),
#' sequence_max = max(SeqScore))

"ShortSeq"
#' 
#' Protein sequence
#'
#' A protein sequence.
#'
#' @format A character string with 219 characters corresponding to P49755.fasta query in UniProt Data base.
#' @source \url{https://www.uniprot.org/}
#' @examples
#' data(MidSeq)
#' MidSeq
#' nchar(MidSeq)
#' data(dico)
#' MidSeqScore=CharSequence2ScoreSequence(MidSeq,dico)
#' MidSeqScore[1:30]
#' localScoreC(MidSeqScore)$localScore
#' prob1 = scoreSequences2probabilityVector(list(MidSeqScore))
#' daudin(localScore = 52, sequence_length = nchar(MidSeq),
#'                score_probabilities = prob1,
#'                sequence_min = min(MidSeqScore),
#'                sequence_max = max(MidSeqScore))
#' score=-5:5
#' prob2=c(0.15,0.15,0.1,0.1,0.0,0.05,0.15,0.05,0.2,0.0,0.05)
#' daudin(localScore = 52, sequence_length = nchar(MidSeq),
#'        score_probabilities = prob2,
#'        sequence_min = min(MidSeqScore),
#'        sequence_max = max(MidSeqScore))
"MidSeq"
#' 
#' Long protein sequence
#'
#' A long protein sequence.
#'
#' @format A character string with 1093 characters corresponding to Q60519.fasta in UniProt Data base.
#' @source \url{https://www.uniprot.org/}
#' @examples
#' data(LongSeq)
#' LongSeq
#' nchar(LongSeq)
#' data(dico)
#' LongSeqScore=CharSequence2ScoreSequence(LongSeq,dico)
#' LongSeqScore[1:50]
#' localScoreC(LongSeqScore)$localScore
#' LS=localScoreC(LongSeqScore)$localScore[1]
#' prob1 = scoreSequences2probabilityVector(list(LongSeqScore))
#' daudin(localScore = LS, sequence_length = nchar(LongSeq),
#'                score_probabilities = prob1,
#'                sequence_min = min(LongSeqScore),
#'                sequence_max = max(LongSeqScore))
#' karlin(localScore = LS, sequence_length = nchar(LongSeq),
#' score_probabilities = prob1,
#' sequence_min = min(LongSeqScore),
#' sequence_max = max(LongSeqScore))
"LongSeq"
#' 
#' Dictionnaire
#'
#' Provides the score related to each base of the sequences.
#'
#' @format A score function for the 20 amino acid
#' @source Kyte & Doolittle (1982) J. Mol. Biol. 157, 105-132
#' @examples
#' data(dico)
#' dico
#' data(MidSeq)
#' MidSeq
#' MidSeqScore=CharSequence2ScoreSequence(MidSeq,dico)
#' MidSeqScore[1:30]
#' localScoreC(MidSeqScore)$localScore
"dico"
#'
#' Several sequences
#'
#' A vector of character strings
#'
#' @format A list of 285 character strings with their entry codes as names
#' @source Structural Classification Of Proteins database (SCOP). More precisely this data contain the 285 protein sequences of the data called "CF_scop2dom_20140205aa" with length from 31 to 404.
#' @examples
#' data(MySeqList)
#' head(MySeqList)
#' MySeqList[1]
#' nchar(MySeqList[1])
#' summary(sapply(MySeqList, nchar))
#' data(dico)
#' MySeqScoreList=lapply(MySeqList, FUN=CharSequence2ScoreSequence, dico)
#' head(MySeqScoreList)
#' AA=automatic_analysis(sequences=MySeqScoreList, model='iid')
#' AA[[1]]
#' # the p-value of the first 10 sequences 
#' sapply(AA, function(x){x$`p-value`})[1:10]
#' # the 20th smallest p-values
#' sort(sapply(AA, function(x){x$`p-value`}))[1:20]
#' which(sapply(AA, function(x){x$`p-value`})<0.05)
#' table(sapply(AA, function(x){x$`method`}))
#' # The maximum sequence length equals 404 so it here normal that the exact method is used for
#' # all the 606 sequences of the data base 
#' # Score distribution learnt on the data set
#' scoreSequences2probabilityVector(MySeqScoreList)
"MySeqList"