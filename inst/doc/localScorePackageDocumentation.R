## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(localScore)

## -----------------------------------------------------------------------------
library(localScore)
help(localScore)
ls(pos = 2)
sequence <- sample(-2:2, 100, replace = TRUE, prob = c(0.5, 0.3, 0.05, 0.1, 0.05))
sequence
localScoreC(sequence)

## -----------------------------------------------------------------------------
library(localScore)
sequence <- sample(c(-3, 2, 0, 1, 5), 100, replace = TRUE, prob = c(0.5, 0.3, 0.05, 0.1, 0.05))
head(sequence)
localScoreC(sequence)

## -----------------------------------------------------------------------------
score_reels <- c(-1, -0.5, 0, 0.5, 1)
proba_score_reels <- c(0.2, 0.3, 0.1, 0.2, 0.2)
sample_from_model <- function(score.sple, proba.sple, length.sple) {
  sample(score.sple,
    size = length.sple, prob = proba.sple, replace = TRUE
  )
}
seq.essai <- sample_from_model(score.sple = score_reels, proba.sple = proba_score_reels, length.sple = 10)
localScoreC_double(seq.essai)

## -----------------------------------------------------------------------------
# Loading a fasta protein
data(Seq219)
Seq219
# or using your own fasta sequence
# MySeqAA_P49755 <- as.character(read.table(file="P49755.fasta",skip=1)[,1])
# MySeqAA_P49755
# Loading a scoring function
data(HydroScore)
?HydroScore
# or using your own scoring function
# HydroScoreKyte<- loadScoreFromFile("Kyte1982.txt")
# Transforming the amino acid sequence into a score sequence with the score function in HydroScore file
SeqScore_P49755 <- CharSequence2ScoreSequence(Seq219, HydroScore)
head(SeqScore_P49755)
length(SeqScore_P49755)
# Computing the local score
localScoreC(SeqScore_P49755)

## ----  fig.width=7.2, fig.height=4--------------------------------------------
monteCarlo(local_score = 10, FUN = function(x) {
  return(sample(x = x, size = length(x), replace = TRUE))
}, x = sequence)

## -----------------------------------------------------------------------------
# Example
score_reels <- c(-1, -0.5, 0, 0.5, 1)
proba_score_reels <- c(0.2, 0.3, 0.1, 0.2, 0.2)
sample_from_model <- function(score.sple, proba.sple, length.sple) {
  sample(score.sple,
    size = length.sple, prob = proba.sple, replace = TRUE
  )
}
monteCarlo_double(5.5,
  FUN = sample_from_model, plot = TRUE, score.sple = score_reels, proba.sple = proba_score_reels,
  length.sple = 100, numSim = 1000
)

## ----  fig.width=7.2, fig.height=4--------------------------------------------
fu <- function(n, size, prob, shift) {
  rbinom(size = size, n = n, prob = prob) + shift
}
karlinMonteCarlo(12,
  FUN = fu, n = 10000, size = 8, prob = 0.2, shift = -2,
  sequence_length = 1000000, simulated_sequence_length = 10000
)

## ----  fig.width=7.2, fig.height=4--------------------------------------------
score_reels <- c(-1, -0.5, 0, 0.5, 1)
proba_score_reels <- c(0.2, 0.3, 0.1, 0.2, 0.2)
fu <- function(score.sple, proba.sple, length.sple) {
  sample(score.sple,
    size = length.sple, prob = proba.sple, replace = TRUE
  )
}
karlinMonteCarlo(85.5,
  FUN = fu, score.sple = score_reels, proba.sple = proba_score_reels,
  length.sple = 10000, numSim = 10000, sequence_length = 100000,
  simulated_sequence_length = 10000
)

## -----------------------------------------------------------------------------
daudin(
  localScore = 15, sequence_length = 500, score_probabilities =
    c(0.2, 0.3, 0.3, 0.0, 0.1, 0.1), sequence_min = -3, sequence_max = 2
)

## -----------------------------------------------------------------------------
score_reels <- c(-1, -0.5, 0, 0.5, 1)
proba_score_reels <- c(0.2, 0.3, 0.1, 0.2, 0.2)
sample_from_model <- function(score.sple, proba.sple, length.sple) {
  sample(score.sple,
    size = length.sple, prob = proba.sple, replace = TRUE
  )
}
seq.essai <- sample_from_model(score.sple = score_reels, proba.sple = proba_score_reels, length.sple = 100)
localScoreC_double(seq.essai)
C <- 10 # homothetic coefficient
localScoreC(C * seq.essai)

## -----------------------------------------------------------------------------
RealScores2IntegerScores(score_reels, proba_score_reels, coef = C)
M.s.r <- RealScores2IntegerScores(score_reels, proba_score_reels, coef = C)$ExtendedIntegerScore
M.s.prob <- RealScores2IntegerScores(score_reels, proba_score_reels, coef = C)$ProbExtendedIntegerScore
M.SL <- localScoreC(C * seq.essai)$localScore[1]
M.SL
pval.E <- daudin(
  localScore = M.SL, sequence_length = 100, score_probabilities = M.s.prob,
  sequence_min = -10, sequence_max = 10
)
pval.E

## -----------------------------------------------------------------------------
SL.real <- localScoreC_double(seq.essai)$localScore[1]
SL.real
pval.MC <- monteCarlo_double(
  local_score = SL.real, FUN = sample_from_model,
  score.sple = score_reels, proba.sple = proba_score_reels, length.sple = 100, plot = TRUE, numSim = 10000
)
pval.MC

## -----------------------------------------------------------------------------
score.v <- -2:1
score.p <- c(0.3, 0.2, 0.2, 0.3)
mean(score.v * score.p)
karlin(
  localScore = 14, sequence_length = 100000, sequence_min = -2, sequence_max = 1,
  score_probabilities = c(0.3, 0.2, 0.2, 0.3)
)
karlin(
  localScore = 14, sequence_length = 1000, sequence_min = -2, sequence_max = 1,
  score_probabilities = c(0.3, 0.2, 0.2, 0.3)
)

## -----------------------------------------------------------------------------
# With missing score values
karlin(
  localScore = 14, sequence_length = 1000, sequence_min = -3, sequence_max = 1,
  score_probabilities = c(0.3, 0.2, 0.0, 0.2, 0.3)
)

## -----------------------------------------------------------------------------
mcc(
  localScore = 14, sequence_length = 1000, sequence_min = -3, sequence_max = 2,
  score_probabilities = c(0.2, 0.3, 0.3, 0.0, 0.1, 0.1)
)

daudin(
  localScore = 14, sequence_length = 1000, score_probabilities =
    c(0.2, 0.3, 0.3, 0.0, 0.1, 0.1), sequence_min = -3, sequence_max = 2
)
karlin(
  localScore = 14, sequence_length = 1000, sequence_min = -3, sequence_max = 2,
  score_probabilities = c(0.2, 0.3, 0.3, 0.0, 0.1, 0.1)
)

## -----------------------------------------------------------------------------
automatic_analysis(sequences = list("x1" = c(1, -2, 2, 3, -2, 3, -3, -3, -2)), model = "iid")

## -----------------------------------------------------------------------------
score <- c(-2, -1, 0, 1, 2)
proba_score <- c(0.2, 0.3, 0.1, 0.2, 0.2)
sum(score * proba_score)
sample_from_model <- function(score.sple, proba.sple, length.sple) {
  sample(score.sple,
    size = length.sple, prob = proba.sple, replace = TRUE
  )
}
seq.essai <- sample_from_model(score.sple = score, proba.sple = proba_score, length.sple = 5000)
MyAnalysis <- automatic_analysis(
  sequences = list("x1" = seq.essai),
  distribution = proba_score, score_extremes = c(-2, 2), model = "iid"
)$x1
MyAnalysis$"p-value"
MyAnalysis$"method applied"
MyAnalysis$localScore$localScore

## -----------------------------------------------------------------------------
score_reels <- c(-1, -0.5, 0, 0.5, 1)
proba_score_reels <- c(0.2, 0.3, 0.1, 0.2, 0.2)
sample_from_model <- function(score.sple, proba.sple, length.sple) {
  sample(score.sple,
    size = length.sple, prob = proba.sple, replace = TRUE
  )
}
seq.essai <- sample_from_model(score.sple = score_reels, proba.sple = proba_score_reels, length.sple = 1000)

# Homothetie
RealScores2IntegerScores(score_reels, proba_score_reels, coef = C)
M.s.r <- RealScores2IntegerScores(score_reels, proba_score_reels, coef = C)$ExtendedIntegerScore
M.s.prob <- RealScores2IntegerScores(score_reels, proba_score_reels, coef = C)$ProbExtendedIntegerScore
# The analysis
MyAnalysis <- automatic_analysis(
  sequences = list("x1" = 10 * seq.essai), model = "iid",
  distribution = M.s.prob, score_extremes = range(M.s.r)
)
MyAnalysis$x1$"p-value"
MyAnalysis$x1$"method applied"

# Without the homothethie, the function gives a wrong result
MyAnalysis2 <- automatic_analysis(sequences = list("x1" = seq.essai), model = "iid")
MyAnalysis2$x1$"p-value"
MyAnalysis2$x1$"method applied"

## ----fig.width=7.2------------------------------------------------------------
set.seed(1)
sequence <- sample(c(-3, 2, 0, 1, 5), 100, replace = TRUE, prob = c(0.5, 0.3, 0.05, 0.1, 0.05))
lindley(sequence)
plot(lindley(sequence), type = "s")

## -----------------------------------------------------------------------------
seq1 <- sample(7:8, size = 10, replace = TRUE)
seq2 <- sample(2:3, size = 15, replace = TRUE)
l <- list(seq1, seq2)
r <- scoreSequences2probabilityVector(l)
r
length(r)

## -----------------------------------------------------------------------------
data(Seq219)
MySeq <- Seq219
data(HydroScore)
SeqScore <- CharSequence2ScoreSequence(MySeq, HydroScore)
n <- length(SeqScore)
n

## -----------------------------------------------------------------------------
LS <- localScoreC(SeqScore)$localScore[1]
LS

## -----------------------------------------------------------------------------
prob <- scoreSequences2probabilityVector(list(SeqScore))
prob

## -----------------------------------------------------------------------------
time.daudin <- system.time(
  res.daudin <- daudin(
    localScore = LS, sequence_length = n,
    score_probabilities = prob,
    sequence_min = min(SeqScore),
    sequence_max = max(SeqScore)
  )
)
res.daudin

## -----------------------------------------------------------------------------
time.karlin <- system.time(
  res.karlin <- karlin(
    localScore = LS, sequence_length = n,
    score_probabilities = prob,
    sequence_min = min(SeqScore),
    sequence_max = max(SeqScore)
  )
)
res.karlin

## -----------------------------------------------------------------------------
time.mcc <- system.time(
  res.mcc <- mcc(
    localScore = LS, sequence_length = n,
    score_probabilities = prob,
    sequence_min = min(SeqScore),
    sequence_max = max(SeqScore)
  )
)
res.mcc

## -----------------------------------------------------------------------------
time.MonteCarlo1 <- system.time(
  res.MonteCarlo1 <- monteCarlo(
    local_score = LS,
    FUN = function(x) {
      return(sample(
        x = x, size = length(x),
        replace = TRUE
      ))
    },
    x = SeqScore, numSim = 200
  )
)
res.MonteCarlo1

## -----------------------------------------------------------------------------
time.MonteCarlo2 <- system.time(
  res.MonteCarlo2 <- monteCarlo(
    local_score = LS,
    FUN = function(x) {
      return(sample(
        x = x, size = length(x),
        replace = TRUE
      ))
    },
    x = SeqScore, numSim = 10000
  )
)
res.MonteCarlo2

## -----------------------------------------------------------------------------
res.pval <- c(
  Daudin = res.daudin, Karlin = res.karlin, MCC = res.mcc,
  MonteCarlo1 = res.MonteCarlo1, MonteCarlo1 = res.MonteCarlo2
)
names(res.pval) <- c("Exact", "Approximation", "Improved appx", "MonteCarlo1", "MonteCarlo2")
res.pval

## -----------------------------------------------------------------------------
rbind(time.daudin, time.karlin, time.mcc, time.MonteCarlo1, time.MonteCarlo2)

## -----------------------------------------------------------------------------
data(Seq31)
MySeq.Short <- Seq31
SeqScore.Short <- CharSequence2ScoreSequence(MySeq.Short, HydroScore)
n.short <- length(SeqScore.Short)
n.short

## -----------------------------------------------------------------------------
SeqScore.S <- SeqScore.Short
LS.S <- localScoreC(SeqScore.S)$localScore[1]
prob.S <- scoreSequences2probabilityVector(list(SeqScore.S))

LS.S
prob.S

## -----------------------------------------------------------------------------
time.daudin <- system.time(
  res.daudin. <- daudin(
    localScore = LS.S, sequence_length = n.short,
    score_probabilities = prob.S,
    sequence_min = min(SeqScore.S),
    sequence_max = max(SeqScore.S)
  )
)

time.karlin <- system.time(
  res.karlin <- try(karlin(
    localScore = LS.S, sequence_length = n.short,
    score_probabilities = prob.S,
    sequence_min = min(SeqScore.S),
    sequence_max = max(SeqScore.S)
  ))
)

time.mcc <- system.time(
  res.mcc <- try(mcc(
    localScore = LS.S, sequence_length = n.short,
    score_probabilities = prob.S,
    sequence_min = min(SeqScore.S),
    sequence_max = max(SeqScore.S)
  ))
)

time.karlinMonteCarlo <- system.time(
  res.karlinMonteCarlo <-
    karlinMonteCarlo(
      local_score = LS.S, plot = FALSE,
      sequence_length = n.short,
      simulated_sequence_length = 1000,
      FUN = sample, x = min(SeqScore.S):max(SeqScore.S),
      size = 1000, prob = prob.S, replace = TRUE,
      numSim = 100000
    )
)

time.MonteCarlo <- system.time(
  res.MonteCarlo <- monteCarlo(
    local_score = LS.S, plot = FALSE,
    FUN = function(x) {
      return(sample(
        x = x, size = length(x),
        replace = TRUE
      ))
    },
    x = SeqScore.S, numSim = 10000
  )
)

## -----------------------------------------------------------------------------
res.pval <- c(Daudin = res.daudin, MonteCarlo = res.MonteCarlo)
names(res.pval) <- c("Daudin", "MonteCarlo")
res.pval
rbind(time.daudin, time.MonteCarlo)

## -----------------------------------------------------------------------------
set.seed(1)
prob.bis <- dnorm(-5:5, mean = -0.5, sd = 1)
prob.bis <- prob.bis / sum(prob.bis)
names(prob.bis) <- -5:5
# Score Expectation
sum((-5:5) * prob.bis)

time.mcc <- system.time(
  res.mcc <- mcc(
    localScore = LS.S, sequence_length = n.short,
    score_probabilities = prob.bis,
    sequence_min = min(SeqScore.S),
    sequence_max = max(SeqScore.S)
  )
)

time.daudin <- system.time(
  res.daudin <- daudin(
    localScore = LS.S, sequence_length = n.short,
    score_probabilities = prob.bis,
    sequence_min = min(SeqScore.S),
    sequence_max = max(SeqScore.S)
  )
)

simu <- function(n, p) {
  return(sample(x = -5:5, size = n, replace = TRUE, prob = p))
}
time.MonteCarlo <- system.time(
  res.MonteCarlo <-
    monteCarlo(
      local_score = LS.S, plot = FALSE,
      FUN = simu, n.short, prob.bis, numSim = 100000
    )
)

## -----------------------------------------------------------------------------
res.pval <- c(MCC = res.mcc, Daudin = res.daudin, MonteCarlo = res.MonteCarlo)
names(res.pval) <- c("MCC", "Daudin", "MonteCarlo")
res.pval
rbind(time.mcc, time.daudin, time.MonteCarlo)

## -----------------------------------------------------------------------------
data(Seq1093)
MySeq.Long <- Seq1093
SeqScore.Long <- CharSequence2ScoreSequence(MySeq.Long, HydroScore)
n.Long <- length(SeqScore.Long)
n.Long
SeqScore.Long <- CharSequence2ScoreSequence(MySeq.Long, HydroScore)
LS.L <- localScoreC(SeqScore.Long)$localScore[1]
LS.L
prob.L <- scoreSequences2probabilityVector(list(SeqScore.Long))
prob.L
sum(prob.L * as.numeric(names(prob.L)))

## ----echo=FALSE---------------------------------------------------------------
time.daudin.L <- system.time(
  res.daudin.L <- daudin(
    localScore = LS.L, sequence_length = n.Long,
    score_probabilities = prob.L,
    sequence_min = min(SeqScore.Long),
    sequence_max = max(SeqScore.Long)
  )
)
res.daudin.L

time.karlin.L <- system.time(
  res.karlin.L <- karlin(
    localScore = LS.L, sequence_length = n.Long,
    score_probabilities = prob.L,
    sequence_min = min(SeqScore.Long),
    sequence_max = max(SeqScore.Long)
  )
)

time.mcc.L <- system.time(
  res.mcc.L <- mcc(
    localScore = LS.L, sequence_length = n.Long,
    score_probabilities = prob.L,
    sequence_min = min(SeqScore.Long),
    sequence_max = max(SeqScore.Long)
  )
)

time.MonteCarlo.L <- system.time(
  res.MonteCarlo.L <-
    monteCarlo(
      local_score = LS.L,
      FUN = sample, x = min(SeqScore.Long):max(SeqScore.Long),
      size = n.Long, prob = prob.L, replace = TRUE, plot = FALSE,
      numSim = 100000
    )
)
res.MonteCarlo.L

time.karlinMonteCarlo.L <- system.time(
  res.karlinMonteCarlo.L <-
    karlinMonteCarlo(
      local_score = LS.L,
      sequence_length = n.Long,
      simulated_sequence_length = 1000,
      FUN = sample, x = min(SeqScore.Long):max(SeqScore.Long),
      size = 1000, prob = prob.L, replace = TRUE,
      numSim = 100000
    )
)
res.karlinMonteCarlo.L

## -----------------------------------------------------------------------------
res.pval.L <- c(res.daudin.L, res.mcc.L, res.karlin.L, res.karlinMonteCarlo.L$`p-value`, res.MonteCarlo.L)
names(res.pval.L) <- c("Daudin", "MCC", "Karlin", "KarlinMonteCarlo", "MonteCarlo")
res.pval.L

## -----------------------------------------------------------------------------
rbind(
  time.daudin.L, time.karlin.L, time.mcc.L, time.karlinMonteCarlo.L,
  time.MonteCarlo.L
)

## -----------------------------------------------------------------------------
MySeqsList <- list(MySeq, MySeq.Short, MySeq.Long)
names(MySeqsList) <- c("Q09FU3.fasta", "P49755.fasta", "Q60519.fasta")
MySeqsScore <- lapply(MySeqsList, FUN = CharSequence2ScoreSequence, HydroScore)
AA <- automatic_analysis(MySeqsScore, model = "iid")
AA$Q09FU3.fasta
AA$Q09FU3.fasta$`method applied`
AA$Q60519.fasta$`method applied`

## -----------------------------------------------------------------------------
cbind(prob, prob.S, prob.L, "3 sequences" = scoreSequences2probabilityVector(MySeqsScore))

## -----------------------------------------------------------------------------
daudin.bis <- daudin(localScore = LS.S, sequence_length = n.short, score_probabilities = scoreSequences2probabilityVector(MySeqsScore), sequence_max = 5, sequence_min = -5)
daudin.bis
AA$P49755.fasta$`p-value`

# automatic_analysis(sequences=list('MySeq.Short'=MySeq.Short), model='iid', distribution=proba.S)

## -----------------------------------------------------------------------------
library(localScore)
data(HydroScore)
data(SeqListSCOPe)
MySeqScoreList <- lapply(SeqListSCOPe, FUN = CharSequence2ScoreSequence, HydroScore)
head(MySeqScoreList)
AA <- automatic_analysis(sequences = MySeqScoreList, model = "iid")
AA[[1]]
# the p-value of the first 10 sequences
sapply(AA, function(x) {
  x$`p-value`
})[1:10]
# the 20th smallest p-values
sort(sapply(AA, function(x) {
  x$`p-value`
}))[1:20]
which(sapply(AA, function(x) {
  x$`p-value`
}) < 0.05)
table(sapply(AA, function(x) {
  x$`method`
}))
# The maximum sequence length equals 404 so it here normal that the exact method is used for all the 606 sequences of the data base
scoreSequences2probabilityVector(MySeqScoreList)

## -----------------------------------------------------------------------------
MyTransMat <-
  matrix(c(
    0.3, 0.1, 0.1, 0.1, 0.4, 0.2, 0.2, 0.1, 0.2, 0.3, 0.3, 0.4, 0.1, 0.1, 0.1, 0.3, 0.3, 0.1, 0.0, 0.3,
    0.1, 0.1, 0.2, 0.3, 0.3
  ), ncol = 5, byrow = TRUE)

MySeq.CM <- transmatrix2sequence(matrix = MyTransMat, length = 150, score = -2:2)
MySeq.CM
AA.CM <- automatic_analysis(sequences = list("x1" = MySeq.CM), model = "markov")
AA.CM

## -----------------------------------------------------------------------------
Ls.CM <- AA.CM$x1$localScore[[1]][1]
monteCarlo(
  local_score = Ls.CM,
  FUN = transmatrix2sequence, matrix = MyTransMat,
  length = 150, score = -2:2,
  plot = FALSE, numSim = 10000
)

