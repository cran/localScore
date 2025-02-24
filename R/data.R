#' Short protein sequence
#'
#' A short protein sequence of 31 amino acids corresponding to Q09FU3.fasta query in UniProt Data base.
#'
#' @format A character string with 31 characters "MLTITSYFGFLLAALTITSVLFIGLNKIRLI"
#' @usage data(Seq31)
#' @source \url{https://www.uniprot.org/}
#' @examples
#' data(Seq31)
#' Seq31
#' nchar(Seq31)
#' data(HydroScore) 
#' SeqScore <- CharSequence2ScoreSequence(Seq31,HydroScore)
#' SeqScore
#' localScoreC(SeqScore)$localScore
#' LS <- localScoreC(SeqScore)$localScore[1]
#' prob1 <- scoreSequences2probabilityVector(list(SeqScore))
#' daudin(local_score = LS, sequence_length = nchar(Seq31),
#'                score_probabilities = prob1,
#'                sequence_min = min(SeqScore),
#'                sequence_max = max(SeqScore))
#' score <- -5:5
#' prob2 <- c(0.15,0.15,0.1,0.1,0.0,0.05,0.15,0.05,0.2,0.0,0.05)
#' sum(prob2*score)
#' karlin(local_score = LS, sequence_length = nchar(Seq31),
#'        score_probabilities = prob2,
#'        sequence_min = min(SeqScore),
#'        sequence_max = max(SeqScore))
"Seq31"
#' 
#' Protein sequence
#'
#' A protein sequence.
#'
#' @format A character string with 219 characters corresponding to P49755.fasta query in UniProt Data base.
#' @usage data(Seq219)
#' @source \url{https://www.uniprot.org/}
#' @examples
#' data(Seq219)
#' Seq219
#' nchar(Seq219)
#' data(HydroScore)
#' seqScore=CharSequence2ScoreSequence(Seq219,HydroScore)
#' seqScore[1:30]
#' localScoreC(seqScore)$localScore
#' prob1 <- scoreSequences2probabilityVector(list(seqScore))
#' daudin(local_score = 52, sequence_length = nchar(Seq219),
#'                score_probabilities = prob1,
#'                sequence_min = min(seqScore),
#'                sequence_max = max(seqScore))
#' score <- -5:5
#' prob2 <- c(0.15,0.15,0.1,0.1,0.0,0.05,0.15,0.05,0.2,0.0,0.05)
#' daudin(local_score = 52, sequence_length = nchar(Seq219),
#'        score_probabilities = prob2,
#'        sequence_min = min(seqScore),
#'        sequence_max = max(seqScore))
"Seq219"
#' 
#' Long protein sequence
#'
#' A long protein sequence.
#'
#' @format A character string with 1093 characters corresponding to Q60519.fasta in UniProt Data base.
#' @usage data(Seq1093)
#' @source \url{https://www.uniprot.org/}
#' @examples
#' data(Seq1093)
#' Seq1093
#' nchar(Seq1093)
#' data(HydroScore)
#' seqScore <- CharSequence2ScoreSequence(Seq1093,HydroScore)
#' seqScore[1:50]
#' localScoreC(seqScore)$localScore
#' LS <- localScoreC(seqScore)$localScore[1]
#' prob1 <- scoreSequences2probabilityVector(list(seqScore))
#' daudin(local_score = LS, sequence_length = nchar(Seq1093),
#'        score_probabilities = prob1,
#'        sequence_min = min(seqScore),
#'        sequence_max = max(seqScore))
#' karlin(local_score = LS, sequence_length = nchar(Seq1093),
#'        score_probabilities = prob1,
#'        sequence_min = min(seqScore),
#'        sequence_max = max(seqScore))
"Seq1093"
#' 
#' Dictionary
#'
#' Provides integer scores related to an hydrophobicity level of each amino acid. This score function is inspired by the Kyte and Doolittle (1982) scale.
#'
#' @format A score function for the 20 amino acid
#' @usage data(HydroScore)
#' @source Kyte & Doolittle (1982) J. Mol. Biol. 157, 105-132
#' @examples
#' data(HydroScore)
#' HydroScore
#' data(Seq219)
#' Seq219
#' seqScore <- CharSequence2ScoreSequence(Seq219,HydroScore)
#' seqScore[1:30]
#' localScoreC(seqScore)$localScore
"HydroScore"
#'
#' Several sequences
#'
#' A vector of character strings
#'
#' @format A list of 285 character strings with their entry codes as names
#' @usage data(SeqListSCOPe)
#' @source Structural Classification Of Proteins database (SCOP). More precisely this data contain the 285 protein sequences of the data called "CF_scop2dom_20140205aa" with length from 31 to 404.
#' @examples
#' data(SeqListSCOPe)
#' head(SeqListSCOPe)
#' SeqListSCOPe[1]
#' nchar(SeqListSCOPe[1])
#' summary(sapply(SeqListSCOPe, nchar))
#' data(HydroScore)
#' MySeqScoreList <- lapply(SeqListSCOPe, FUN=CharSequence2ScoreSequence, HydroScore)
#' head(MySeqScoreList)
#' AA <- automatic_analysis(sequences=MySeqScoreList, model='iid')
#' AA[[1]]
#' # the p-value of the first 10 sequences 
#' sapply(AA, function(x){x$`p-value`})[1:10]
#' # the 20th smallest p-values
#' sort(sapply(AA, function(x){x$`p-value`}))[1:20]
#' which(sapply(AA, function(x){x$`p-value`})<0.05)
#' table(sapply(AA, function(x){x$`method`}))
#' # The maximum sequence length equals 404 so it here normal that the exact method is used for
#' # all the 606 sequences of the data base 
#' # Score distribution learned on the data set
#' scoreSequences2probabilityVector(MySeqScoreList)
"SeqListSCOPe"
#'
#' Congenital oesophageal atresia data
#'
#' The data consists of individual dates of birth over n=35 cases of the birth defects oesophageal and tracheo-oesophagean fistula
#' observed in a hospital in Birmingham, U.K., over 2191 days from 1950 through 1955, with Day one set as 1 January 1950
#'
#' @format A matrix of 2191 lines and 2 columns. Each line is a day on the first column, and associated to a case (0/1) on the second column.
#' @usage data(Aeso)
#' @source Dolk H. Secular pattern of congenital oesophageal atresia--George Knox, 1959. J Epidemiol Community Health. 1997;51(2):114-115. <doi:10.1136/jech.51.2.114>
#' @examples
#' data(Aeso)
#' head(Aeso)
#' p <- sum(Aeso[,2]) / dim(Aeso)[1]
#' print(p)
"Aeso"
#'
#' Stevens-Johnson syndrome data
#'
#' The Stevens-Johnson syndrome is an acute and serious dermatological disease due to a drug allergy. 
#' The syndrome appearance is life-threatening emergency. They are very rare, around 2 cases per million people per year. 
#'
#' @format A data.frame of 824 lines, each describing a syndrome appearance described by 15 covariates:
#' Case ID,
#' Initial FDA Received Date,
#' days since last fda,
#' Event Date,
#' Latest FDA Received Date,
#' Suspect Product Names,
#' Suspect Product Active Ingredients,
#' Reason for Use,
#' Reactions,
#' Serious,
#' Outcomes,
#' Sex,
#' Patient Age,
#' Sender,
#' Concomitant Product Names
#' The third column correspond to the number of days between two adverse events.
#' @usage data(SJSyndrome.data)
#' @source FDA open data.
#' @examples
#' data(SJSyndrome)
#' summary(SJSyndrome)
#' 
"SJSyndrome"
#'
#' Deprecated
#' 
#' Some datasets are deprecated and will be removed in next version of the package.
#' Please use instead:
#' \itemize{
#'   \item data(aeso.data) -> data(Aeso)
#'   \item data(SJSyndrome.data) -> data(SJSyndrome)
#'   \item data(MySeqList) -> data(SeqListSCOPe)
#'   \item data(dico) -> data(HydroScore)
#'   \item data(LongSeq) -> data(Seq1093)
#'   \item data(MidSeq) -> data(Seq219)
#'   \item data(ShortSeq) -> data(Seq31)
#' }
#' @usage data(Aeso) ; data(SJSyndrome) ;
#'        data(SeqListSCOPe) ; data(HydroScore) ;
#'        data(Seq1093) ; data(Seq219) ; data(Seq31)
#' @format Refers to the "See also" links below to obtain a description of each dataset.
#' @aliases SJSyndrome.data MySeqList dico LongSeq MidSeq ShortSeq
#' @seealso \code{\link{Aeso}} \code{\link{SJSyndrome}} \code{\link{SeqListSCOPe}} 
#'          \code{\link{HydroScore}} \code{\link{Seq1093}} \code{\link{Seq219}} \code{\link{Seq31}}
#' @keywords internal
"aeso.data"
