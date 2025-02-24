#' @keywords internal
#' Package for sequence analysis by local score
#'
#' @md
#' @description
#' Provides functionalities for:
#'   * calculating the local score
#'   * calculating statistical relevance (p-value) to find a local Score in a sequence of given distribution
#' 
#' Given a sequence of numerical score \eqn{X_1,\dots,X_n}, the local score is defined : \eqn{H_n = \max_{1 \leq i \leq j \leq n} \sum_{l = i}^j X_l}
#' This package find the value local score and the associated sub-sequence, and also sub-optimal local scores and segments. The complexity is linear with \eqn{n}. It can be viewed as a generalization of a sliding window method, considering all windows size. In order to be pertinent, the expectation of the scores \eqn{X_i} should be negative. Most of the methods concerning statistical relevance implemented in this package only applied on integer scores.
#' 
#' @details
#' Please refer to the vignette of this package or the manual for details on how to use this package.
#' 
#' @references 
#' * An Improved Approximation For Assessing The Statistical Significance of molecular Sequence Features, Mercier and al 2003
#' * Exact distribution for the local score of one i.i.d. random sequence, Sabine Mercier and JJ Daudin, 2001
#' * Limit Distributions of Maximal Segmental Score among Markov-Dependent Partial Sums, Karlin and Dembo 1992
#' * Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes, Karlin and al 1990
#' * Detection de courts segments inverses dans les genomes: methodes et applications, David Robelin 2005
#' * A Linear Time Algorithm for Finding All Maximal Scoring Subsequences, Constable and Bates 1985
#' 
"_PACKAGE"
  
## usethis namespace: start
## usethis namespace: end
NULL
