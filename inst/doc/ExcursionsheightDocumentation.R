## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE--------------------------------------------------------
library(localScore)

## -----------------------------------------------------------------------------
seq <- c(-1,2,-3,-4,-1,3,-1,-2,3,-4,2)
lindley(seq)
plot(1:length(seq), lindley(seq), type = 'l')
localScoreC(seq)

## -----------------------------------------------------------------------------
recordTimes(seq)

## -----------------------------------------------------------------------------
proba_theoretical_ith_excursion_markov(a = 5,
                                       theta = c("A", "B"),
                                       lambda = matrix(c(0.5, 0.5, 0.8, 0.2),
                                                       ncol = 2, byrow = TRUE),
                                       score_function = c(-1,1),
                                       i = 2,
                                       epsilon = 1e-16,
                                       prob0 = c(0.5, 0.5))

## -----------------------------------------------------------------------------
proba_theoretical_ith_excursion_markov(a = 5,
                                       theta = c("K", 12, 1.54),
                                       lambda = matrix(c(0.9, 0, 0.1, 0.9, 0.1, 0, 0, 0.8, 0.2),
                                                       ncol = 3, byrow = TRUE),
                                       score_function = c(-3,-1,2),
                                       i = 4)

## -----------------------------------------------------------------------------
transition_matrix <- matrix(runif(400, min = 0, max = 1), nrow = 20)
transition_matrix <- t(apply(transition_matrix, 1, function(x) x/sum(x)))
theta <- letters[1:20]
score_f <- c(-3,-1,2,-3,-1,2,-3,-1,2,-1,
              -3,-1,2,-3,-1,2,-3,-1,2,-1)
sum(stationary_distribution(transition_matrix)*score_f) #score expectation (stationary)
system.time(pv1 <- proba_theoretical_ith_excursion_markov(a=5, theta, transition_matrix, score_f, i = 4))
pv1$proba_q_i_geq_a

## -----------------------------------------------------------------------------
library(localScore)
data("Seq1093")
Seq1093
nchar(Seq1093)
data(HydroScore)
LongSeqScore <- CharSequence2ScoreSequence(Seq1093, HydroScore)
table(LongSeqScore)
localScoreC(LongSeqScore)

head(sort(localScoreC(LongSeqScore)$suboptimalSegmentScores[,1], decreasing = TRUE))
plot(1:1093,lindley(LongSeqScore), type = 'l')
plot(1:400,lindley(LongSeqScore)[1:400], type = 'l')

## -----------------------------------------------------------------------------
LS <- localScoreC(LongSeqScore)$localScore["value"]
prob1 <- scoreSequences2probabilityVector(list(LongSeqScore))
prob1
daudin(local_score = LS, sequence_length = length(LongSeqScore),
       score_probabilities = prob1,
       sequence_min = min(LongSeqScore),
       sequence_max = max(LongSeqScore))

## -----------------------------------------------------------------------------
tmpMarkovParameters <- sequences2transmatrix(LongSeqScore)
lambda <- tmpMarkovParameters$transition_matrix
apply(lambda, 1, sum) #to check stochasticity
prob0 <- stationary_distribution(lambda)
print(prob0)
score_values <- tmpMarkovParameters$score_value
print(score_values)

## -----------------------------------------------------------------------------
exact_mc(LS, lambda, sequence_length = length(LongSeqScore), score_values = score_values)

## -----------------------------------------------------------------------------
subOptSegment <- localScoreC(LongSeqScore)$suboptimalSegmentScores
o <- order(subOptSegment[,1], decreasing = TRUE)
subOptSegment <- subOptSegment[o,] # reordering segments by decreasing score values
print(subOptSegment)

## -----------------------------------------------------------------------------
lindley(LongSeqScore)
recordTimes(LongSeqScore)
LongSeqScore[1] > 0

## -----------------------------------------------------------------------------
subOptSegment["1",] #First excursion
subOptSegment[as.character(dim(subOptSegment)[1]),] #Last excursion

## -----------------------------------------------------------------------------
theta <- letters[1:length(score_values)] # arbitrary
score_function <- score_values           # defined earlier 
a <- 40 
i <- 1
system.time(pv2<-proba_theoretical_ith_excursion_markov(a, theta, lambda, 
                                       score_function,i)$proba_q_i_geq_a)
pv2

a <- 6 
i <- 77
system.time(pv3<-proba_theoretical_ith_excursion_markov(a, theta, lambda, 
                                       score_function, i)$proba_q_i_geq_a)
pv3

## -----------------------------------------------------------------------------
a <- 6 
i <- 20
system.time(pv4 <- proba_theoretical_ith_excursion_markov(a, theta, 
                                       lambda, 
                                       score_function,i)$proba_q_i_geq_a)
pv4

i <-  10
system.time(pv5 <- proba_theoretical_ith_excursion_markov(a, theta, 
                                       lambda, 
                                       score_function,i)$proba_q_i_geq_a)
pv5

## -----------------------------------------------------------------------------
LongSeqScore.inv <- rev(LongSeqScore)
localScoreC(LongSeqScore.inv)
head(sort(localScoreC(LongSeqScore.inv)$suboptimalSegmentScores[,1], decreasing = TRUE))
plot(1:1093,lindley(LongSeqScore.inv), type = 'l')

## -----------------------------------------------------------------------------
markov_parameters <- sequences2transmatrix(LongSeqScore.inv)
lambda.inv <- markov_parameters$transition_matrix
system.time(pv5<-proba_theoretical_ith_excursion_markov(a = 38, theta, 
                                       lambda.inv, 
                                       score_function, i = 96
                                       )$proba_q_i_geq_a)
pv5
proba_theoretical_ith_excursion_markov(a = 6, theta, lambda.inv,
                                       score_function,i = 1
                                       )$proba_q_i_geq_a

## -----------------------------------------------------------------------------
a <- 65 
i <- 30
proba_theoretical_ith_excursion_markov(a, theta, 
                                       lambda, score_function, 
                                       i)$proba_q_i_geq_a


## -----------------------------------------------------------------------------
subOptSegment.inv <- localScoreC(LongSeqScore.inv)$suboptimalSegmentScores
which.max(subOptSegment.inv[,1])
print(subOptSegment.inv[16,])

a <- 65 
i <- 16
proba_theoretical_ith_excursion_markov(a, theta, lambda.inv, 
                                       score_function,i)$proba_q_i_geq_a


