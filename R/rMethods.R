#' @useDynLib localScore
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline barplot par plot
#' @importFrom stats lm ecdf knots density 
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar
#' 
#' @rdname localScoreC 
#' @export
localScoreC_double <- function(v, suppressWarnings = FALSE) {
  mode(v) <- "numeric"
  return(localScoreC(v, suppressWarnings))
}
#' 
#' @rdname localScoreC
#' @export
localScoreC_int <- function(v, suppressWarnings = FALSE) {
  v <- as.integer(v)
  return(localScoreC(v, suppressWarnings))
}
#' 
#' @title Monte Carlo method [p-value]
#' @description Calculates an empirical p-value based on Monte Carlo simulations.
#' Perfect for small sequences (both markov chains and identically and independently distributed) with length ~ 10^3. 
#' @param local_score local score observed in a segment.
#' @param FUN function to simulate similar sequences with.
#' @param ... parameters for FUN
#' @param plot boolean value if to display plots for cumulated function and density
#' @param numSim number of sequences to generate during simulation
#' @param keepSimu Boolean, default to FALSE. If TRUE, the simulated local scores are returned as the localScores element of the output list.
#' @return If \code{keepSimu} is FALSE, returns a numeric value corresponding to the probability to obtain a local score with value greater or equal to the parameter \code{local_score}. \cr
#' If \code{keepSimu} is TRUE, returns a list containing:\tabular{ll}{
#'     \code{p_value} \tab Floating value corresponding to the probability to obtain a local score with a value greater or equal to the parameter \code{local_score} \cr
#'     \code{localScores} \tab Vector of size \code{numSim} containing the simulated local scores
#' }
#' @details
#' Be careful that the parameters names of the function FUN should differ from those of monteCarlo function.\cr
#' The density plot produced by \code{plot == TRUE} depends on the type of the simulated local scores: 
#' if they are integer, a barplot of relative frequency is used, else \code{plot(density(...))} is used. \cr
#' This function calls \code{\link{localScoreC}} which type of the output depends on the type of the input. 
#' To be efficient, be aware to use a simulating function \code{FUN} that return a vector of adequate type ("integer" or "numeric"). Warning: in R, \code{typeof(c(1,3,4,10)) == "double"}. You can set a type of a vector with \code{mode()} or \code{as.integer()} functions for example. \cr 
#' \code{monteCarlo_double()} is deprecated. At this point, it is just a call to \code{monteCarlo()} function.
#' @seealso \code{\link{karlinMonteCarlo}} \code{\link{localScoreC}}
#' @examples
#' \donttest{
#' monteCarlo(120, FUN = rbinom, n = 100, size = 5, prob=0.2)
#' }
#' \donttest{
#' mySeq <- sample(-7:3, replace = TRUE, size = 1000)
#' monteCarlo(local_score = 18, FUN = function(x) {return(sample(x = x, 
#'            size = length(x), replace = TRUE))}, x = mySeq)
#' }
#' \donttest{
#' #Examples of non integer score function
#' mySeq2 <- sample(-7:6 - 0.5, replace = TRUE, size = 1000)
#' monteCarlo(local_score = 50.5, FUN = function(x) {return(sample(x = x, 
#'            size = length(x), replace = TRUE))}, x = mySeq2)
#' }
#' \donttest{
#' #Examinating simulated local scores
#' mySeq2 <- sample(-7:6, replace = TRUE, size = 1000)
#' simu <- monteCarlo(local_score = 50.5, FUN = function(x) {return(sample(x = x, 
#'            size = length(x), replace = TRUE))}, x = mySeq2, keepSimu = TRUE)
#' hist(simu$localScores)
#' }
#' \donttest{
#' # Markovian example
#' MyTransMat <-
#'     matrix(c(0.3,0.1,0.1,0.1,0.4, 0.2,0.2,0.1,0.2,0.3, 0.3,0.4,0.1,0.1,0.1, 0.3,0.3,0.1,0.0,0.3,
#'              0.1,0.1,0.2,0.3,0.3), ncol = 5, byrow=TRUE)
#' monteCarlo(local_score = 50,
#'           FUN = transmatrix2sequence, matrix = MyTransMat,
#'           length=150, score = c(-2,-1,0,2,3), plot=FALSE, numSim = 5000)
#' }
#' @export
monteCarlo <- function(local_score, FUN, ..., plot = TRUE, numSim = 1000, keepSimu = FALSE){
  epsilon <- 1e-8 #Awfull trick used to correct CdF in case of discrete score values
  if (local_score < 0)
    stop("[Invalid Input] local scores are never strictly negative.")
  if (local_score == 0) #Trivial case
    return(1.0)
  if (numSim <= 0)
    stop("[Invalid Input] Please supply a positive number of simulations.")
  mcParametersNames <- c("local_score", "FUN", "numSim", "plot", "keepSimu")
  FUNparametersName <- names(formals(FUN))
  if (any(FUNparametersName %in% mcParametersNames))
    stop("[Invalid Input] Parameters of FUN should differ from parameters of monteCarlo function")

  #simulate sequences to get some local scores to work with
  maxScores <- sim_local_scores(numSim = numSim, FUN = FUN, ...)
  if (!all(!is.na(maxScores)))
    warning(paste0("Some simulated local scores are NA. Number of NA: ",sum(is.na(maxScores))))
  
  #Empirical cumulated distribution of found local Scores
  scoreCdf <- ecdf(maxScores)
  if (all(maxScores %% 1 == 0)) { # is.integer(maxScores) can lead to misleading 
    scoreD <- table(maxScores)/length(maxScores)
  } else if (is.numeric(maxScores)) {
    scoreD <- density(maxScores)
  } else {
    stop("The local score associated to simulating function should be either an integer or real number.")
  }
  #if the user wants to display the distribution of local scores, print it to the plot screen
  if (plot) {
    #set plot window
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    par(mfrow = c(1,2))
    plot(scoreD, main = "Distribution of local scores\n for given sequence", ylab = "Probability", xlab = "Local score")
    abline(v = local_score, col = "red", lty = 3) 
    plot(scoreCdf, main = "Cumulative Distribution Function", ylab = "cumulative distribution(local_scores)))", xlab = "Local score")
    abline(v = local_score, col = "red", lty = 3)  #greater or equal, thus set marker at localScore -1
    abline(h = 1 - scoreCdf(local_score * (1 - epsilon)), col = "red", lty = 3)
  }
  result <- 1 - scoreCdf(local_score * (1 - epsilon))
  names(result) <- c("p_value")
  if (keepSimu) {
    return(list("p_value" = result, localScores = maxScores))
  } else {
    return(result) 
  }
}

#' @title Simulate a sequence of local score [p-value]
#' @description Simulate a sequence of local score based on a simulation function of individual sequence of score.
#' @param FUN function to simulate similar sequences with.
#' @param ... parameters for FUN
#' @param numSim number of sequences to generate during simulation
#' @return vector containing a local score sample of size numSim
#' @details
#' If all simulated local score are integer (in mathematical sense), the returned value type is set to \code{integer}
#' @noRd
sim_local_scores <- function(FUN, ..., numSim = 1000) {
  #maxScores <- vector(mode = "numeric", length = numSim)
  maxScores <- rep(NA, length = numSim)
  if (interactive())
    progressBar <- txtProgressBar(min = 0, max = numSim, initial = 0, char = "=",
                                  width = NA, 
                                  paste0("progress in local scores simulation (", numSim, ")"),
                                  "progressbar", style = 3, file = "")
  for (i in 1:numSim) {
      simSequence <- FUN(...)
      #print(paste0("length simSequence: ",length(simSequence)))
      maxScores[i] <- localScoreC(simSequence, suppressWarnings = TRUE)$localScore["value"]
      if (interactive())
        setTxtProgressBar(progressBar, i)
  }
  if (all(maxScores %% 1 == 0))  # is.integer(maxScores) is misleading
    mode(maxScores) <- "integer"
  return(maxScores)
}

#' @describeIn monteCarlo Monte-Carlo function for double [deprecated]
#' @export
monteCarlo_double <- function(local_score, FUN, ..., plot = TRUE, numSim = 1000, keepSimu = FALSE) {
  return(monteCarlo(local_score, FUN, ..., plot = plot, numSim = numSim, keepSimu = keepSimu))
}

#' @title Lindley process
#' @description Creates a sequence of a Lindley process, also called CUSUM process, on a given sequence.  For a sequence (X_k)k, the Lindley process is defined as follows: W_0:=0 and W_(k+1)=max(0,W_k+X_(k+1)). It defines positive excursions above 0.
#' @param sequence numeric sequence of a Lindley process, eg service time per customer
#' @return a vector with the Lindley process steps
#' @examples
#' MySeq <- c(1,2,3,-4,1,-3,-1,2,3,-4,1)
#' lindley(MySeq)
#' plot(1:length(MySeq),lindley(MySeq),type='b')
#' @export
lindley <- function(sequence){
  seq <- rep(NA, length(sequence))
  sum <- 0
  i <- 1
  for (score in sequence) {
    sum <- max(0, score + sum)
    seq[i] <- sum
    i <- i + 1
  }
  return(seq)
}

#' @title Calculate the record times of a sequence
#' @description For a given sequence of real numbers, this function returns the
#'   index of the values where Lindley/CUSUM process are <0
#' @param sequence numeric sequence of a Lindley process, eg service time per
#'   customer
#' @return a vector with the record times. If no record time are found, return
#'   empty vector \code{integer(0)}
#' @details Note that the first record times which is always 0 is not included
#' in the returned vector
#' @examples
#' ####This example should return this vector: c(1,3,4,5,10)
#' seq1 <- c(-1,2,-2.00001,-4,-1,3,-1,-2,3,-4,2)
#' recordTimes(seq1)
#' ####This example should return integer(0) because there is no record times in this sequence
#' seq2 <- c(4,1,0,-1,-2,5,-5,0,4)
#' recordTimes(seq2)
#' @export
recordTimes <- function(sequence) {
  record_times <- rep(NA, length(sequence))
  record_times <- as.vector(record_times, mode = "integer")
  sum <- 0
  i <- 1
  j <- 1
  for (score in sequence) {
    if ( sum + score < 0) {
      record_times[j] <- i
      j <- j + 1
    }
    sum = max(0, sum + score)
    i <- i + 1
  }
  record_times <- record_times[!is.na(record_times)]
  return(record_times)
}

#' @title Monte Carlo - Karlin [p-value]
#' @description Estimates p-value of the local score based on a Monte Carlo
#'   estimation of Gumble parameters from simulations of smaller sequences with
#'   same distribution. Appropriate for great sequences with length > 10^3, for
#'   i.i.d and markovian sequence models.
#' @param local_score local score observed in a sequence.
#' @param sequence_length length of the sequence
#' @param simulated_sequence_length length of simulated sequences produced by
#' @param FUN function to simulate similar sequences with.
#' @param ... parameters for FUN
#' @param numSim number of sequences to create for estimation
#'   FUN
#' @param plot boolean value if to display plots for cumulated function, density
#'   and linearization of cumulative density function
#' @param keepSimu Boolean, default to FALSE. If TRUE, the simulated local
#'   scores are returned as the localScores element of the output list.
#' @return If \code{keepSimu} is FALSE, returns a list containing:\tabular{ll}{
#'   \code{p_value} \tab Probability to obtain a local score with a value
#'   greater or equal to the parameter \code{local_score} \cr \tab \cr \code{K*}
#'   \tab Parameter \eqn{K^*} defined in Karlin and Dembo (1990) \cr \tab
#'   \cr \code{lambda} \tab Parameter \eqn{\lambda} defined in Karlin and
#'   Dembo (1990) \cr }
#' If \code{keepSimu} is TRUE, returns a list containing:\tabular{ll}{
#'     \code{p_value} \tab Probability to obtain a local score with a value greater or equal to the parameter \code{local_score} \cr
#'     \tab \cr
#'     \code{K*} \tab Parameter \eqn{K^*} defined in Karlin and Dembo (1990) \cr
#'     \tab \cr
#'     \code{lambda} \tab Parameter \eqn{\lambda} defined in Karlin and Dembo (1990) \cr
#'     \code{localScores} \tab Vector of size \code{numSim} containing the simulated local scores for sequence size of \code{simulated_sequence_length}
#' }
#' @details The length of the simulated sequences is an argument specific to the
#'   function provided for simulation. Thus, it has to be provided also in the
#'   parameter \code{simulated_sequence_length} in the arguments of the "Monte
#'   Carlo - Karlin" function. It is a crucial detail as it influences precision
#'   and computation time of the result. Note that to get an appropriate
#'   estimation, the given average score must be non-positive. Be careful that
#'   the parameters names of the function \code{FUN} should differ from those of
#'   \code{karlinMonteCarlo} function. \cr Methods - Parameters \eqn{K^\star}
#'   and \eqn{\lambda} of Karlin and Dembo (1990) are estimated by a linear
#'   regression on the log(-log(cumulative distribution function of the local
#'   scores)) on shorter sequences (size \code{simulated_sequence_length}). The
#'   formula used are : \eqn{\hat{\lambda} = -\hat{b}} and \eqn{\hat{K^\star} =
#'   exp(\hat{a})/simulated\_sequence\_length} where \eqn{\hat{a}} is the
#'   intercept of the regression and \eqn{\hat{b}} is the slope of the
#'   regression. Then p-value is given by \eqn{p = exp(-K^\star * exp(-\lambda*x
#'   ))} where \eqn{x = local\_score - \log(sequence\_length)/\lambda}. \cr The
#'   density plot produced by \code{plot == TRUE} depends on the type of the
#'   simulated local scores: if they are integer, a barplot of relative
#'   frequency is used, else \code{plot(density(...))} is used. \cr This
#'   function calls \code{\link{localScoreC}} which type of the output depends
#'   on the type of the input. To be efficient, be aware to use a simulating
#'   function \code{FUN} that return a vector of adequate type ("integer" or
#'   "numeric"). Warning: in R, \code{typeof(c(1,3,4,10)) == "double"}. You can
#'   set a type of a vector with \code{mode()} or \code{as.integer()} functions
#'   for example. \cr \code{karlinMonteCarlo_double()} is deprecated. At this
#'   point, it is just a call to \code{karlinMonteCarlo()} function.
#' @seealso \code{\link{monteCarlo}} \code{\link{localScoreC}}
#' @examples
#' \donttest{
#' mySeq <- sample(-7:6, replace = TRUE, size = 100000)
#' #MonteCarlo taking random sample from the input sequence itself
#' karlinMonteCarlo(local_score = 160, sequence_length = 100000,
#'                simulated_sequence_length = 1000,
#'                FUN = function(x, sim_length) {
#'                         return(sample(x = x,
#'                                size = sim_length,
#'                                replace = TRUE))
#'                      },
#'                x = mySeq,
#'                sim_length = 1000,
#'                numSim = 1000)
#' }
#' \donttest{
#' #Markovian example (longer computation)
#' MyTransMat_reels <-  matrix(c(0.3, 0.1, 0.1, 0.1, 0.4,
#'                               0.2, 0.2, 0.1, 0.2, 0.3,
#'                               0.3, 0.4, 0.1, 0.1, 0.1,
#'                               0.3, 0.3, 0.1, 0.2, 0.1,
#'                               0.2, 0.1, 0.2, 0.4, 0.1),
#'                               ncol = 5, byrow=TRUE)
#' karlinMonteCarlo(local_score = 18.5, sequence_length = 200000,
#'                  simulated_sequence_length = 1500,
#'                  FUN = transmatrix2sequence,
#'                  matrix = MyTransMat_reels,
#'                  score =c(-1.5,-0.5,0,0.5,1), length = 1500,
#'                  plot=TRUE, numSim = 1500)
#' }
#' @export
karlinMonteCarlo <- function(local_score, sequence_length, simulated_sequence_length, FUN,  ..., numSim = 1000, plot = TRUE, keepSimu = FALSE) {
  if (local_score < 0)
    stop("[Invalid Input] local Scores are not negative.")
  if (numSim <= 0)
    stop("[Invalid Input] Please supply a positive number of simulations.")
  if (missing(sequence_length))
    stop("[Missing Argument] sequence_length is a required parameter")
  if (missing(simulated_sequence_length))
    stop("[Missing Argument] simulated_sequence_length is a required parameter")
  kmcParametersNames <- c("local_score", "sequence_length", "simulated_sequence_length", "FUN", "numSim", "plot")
  FUNparametersName <- names(formals(FUN))
  if (any(FUNparametersName %in% kmcParametersNames))
    stop("[Invalid Input] Parameters of FUN should differ from parameters of karlinMonteCarlo function")
  #simulate sequences to get some local scores to work with
  #maxScores <- sim_local_scores(numSim = numSim, FUN = FUN, ...) # bug voir #11238 sur forge DGA
  maxScores <- sim_local_scores(numSim = numSim, FUN = FUN, ...)

  #Empirical distribution of found local Scores on (short) simulated sequences
#  scoreCdf <- ecdf(maxScores)
  scoreShortD <- table(maxScores)/length(maxScores)

  # if (all(maxScores %% 1 == 0)) { # is.integer(maxScores) can lead to misleading 
  # } else if (is.numeric(maxScores)) {
  #   scoreShortD <- density(maxScores)
  # } else {
  #   stop("The local score associated to simulating function should be either an integer or real number.")
  # }
  
  #Cumulative distribution function of local score probabilities
  log_cumsum <- log(-log(cumsum(scoreShortD))) 
  #clean up eventual inf values
  log_cumsum <- log_cumsum[is.finite(log_cumsum)]
  
  #linear regression (Warning : at least two different simulated values of MaxScores needed)
  if (length(log_cumsum) < 2) {
    p_value <- NA
    k_star <- NA
    lambda <- NA
  } else {
    reg <- lm(log_cumsum~as.numeric(names(log_cumsum)))
    lambda <- -reg$coefficients[2]      #Inclination
    names(lambda) <- NULL
    k_star <- exp(reg$coefficients[1])/simulated_sequence_length #SIMULATED SEQUENCE LENGTH
    names(k_star) <- NULL
    
    x <- local_score - log(sequence_length)/lambda - 1
    #Karlin: for n great, P( ln(n)/lambda+x>= M) = exp(-K_star*exp(-lambda*x))
    #thus we set ln(n)/lambda+x = local_score and obtain for x = local_score - ln(n)/lambda
    #x = local_score - log(sequence_length)/lambda#-1 ? XXX a priori oui, il faut le '-1' (cf le code C++)
    
    #now we calculate p-value with our approximate K star and lambda
    p_value <- exp(-k_star*exp(-lambda*x))
  }
  if (keepSimu) {
    result <- list("p_value" = 1 - p_value, "K*" = as.numeric(k_star), 
                   "lambda" = as.numeric(lambda), localScores = maxScores)
  } else {
    result <- list("p_value" = 1 - p_value, "K*" = as.numeric(k_star),
                   "lambda" = as.numeric(lambda))
  }
  if (plot) {
    if (is.numeric(maxScores) && !all(maxScores %% 1 == 0))  # real scores so density plot
      scoreShortD <- density(maxScores)
    
    #absice transformations
    # xnew <- as.numeric(names(scoreShortD)) + log(sequence_length/simulated_sequence_length)/lambda
    # names(scoreShortD) <- xnew
    #set plot window to contain 3 slots
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    par(mfrow = (c(1,3)))
    # lab <- rep("", length(xnew))
    # lab[seq.int(from = 1, to = length(lab), by = trunc(length(lab)/2))] <- round(xnew[seq.int(from = 1, to = length(lab), by = trunc(length(lab)/2))],1)
    #lab <- round(xnew, 1) # Une décimale pour les labels
    plot(scoreShortD, main = paste0("Distribution of simulated local scores\n (n=",simulated_sequence_length,")"), ylab = "Empirical density", xlab = "Local score values", frame = TRUE) #, xaxt = "n")
    #axis(1, at = lab, las = 2)
    
    plot(cumsum(table(maxScores)/length(maxScores))~as.numeric(names(table(maxScores))), main = paste0("Cumulative Distribution Function\n (n=",simulated_sequence_length,")"), ylab = "cumulative distribution(local_scores)))", xlab = "Local Score Values")
    plot(log_cumsum~as.numeric(names(log_cumsum)), main = "Linear Regression", ylab = "ln(-ln(cumulative destribution(local_scores)))", xlab = "Local score values")
    abline(reg, col = c("green"))
  }
  return(result)
}

#' @rdname karlinMonteCarlo
#' @export
karlinMonteCarlo_double <- function(local_score, sequence_length, simulated_sequence_length, FUN,  ..., numSim = 1000, plot = TRUE, keepSimu = FALSE){
  return(karlinMonteCarlo(local_score, sequence_length, simulated_sequence_length, FUN,  ..., numSim = numSim, plot = plot, keepSimu = keepSimu))
}

#' @title Transition matrix from sequence(s)
#' @description Calculates the transition matrix by counting occurrences of tuples in given vector list
#' @param sequences Sequences to be analyzed, can be a list of vectors, or a vector
#' @return A list object containing
#' \item{transition_matrix}{Transition Matrix with row names and column names are the associated score/state}
#' \item{score_value}{a vector containing the score/state, ordered in the same way as the matrix columns and rows}
#' @details
#' In output, score_value is coerced as integer if possible. Else, it is a character vector containing the states
#' of the Markov chain.
#' @examples
# seq1 <- sample(-1:1, size = 20, replace = TRUE)
# seq2 <- sample(-6:1, size = 20, replace = TRUE)
# seq3 <- sample(3:6, size = 50, replace = TRUE)
# sequences2transmatrix(list(seq1, seq2, seq3))
#' myseq <- sample(LETTERS[1:2], size=20,replace=TRUE)
#' sequences2transmatrix(myseq)
#' @export
sequences2transmatrix <- function(sequences){
  if (is.list(sequences))
    #get score in 1 vector separate by NA values
    allseq <- unlist(lapply(sequences, FUN = function(x){c(x,NA)}))
  else if (is.vector(sequences))
    allseq <- c(sequences,NA)
  else
    stop("The parameter should be either a list of vector or a vector")
  
  # Count Nuv/Nu
  p <- prop.table(table(allseq[-length(allseq)],allseq[-1]),margin = 1)
  # Treat p name as numeric or not
  state <- suppressWarnings(as.numeric(rownames(p)))
  if (all(!is.na(state))) #case all names are numeric
    score_value <- state
  else
    score_value <- rownames(p)
  if ((length(rownames(p)) != length(colnames(p))) || !all(rownames(p) == colnames(p)))
    stop("The score beginning the sequence appears only one time.")
  if (!all(!is.na(p))) # there is at least 1 NA is the matrix
    stop(("The score ending the sequence appears only one time."))
  return(list("transition_matrix" = p, "score_value" = score_value))
}

#' @title Sampling function for Markov chains
#' @description Creates Markov chains based on a transition matrix. Can be used as parameter for the Monte Carlo function.
#' @details The transition matrix is considered representing the transition from one score to another such that the score in the first row is the lowest 
#' and the last row are the transitions from the highest score to another. The matrix must be stochastic (no rows filled up with only '0' values).
#' @param matrix transition matrix of Markov process
#' @param length length of sequence to be sampled
#' @param initialIndex (optional) index of matrix which should be initial value of sequence. If none supplied, a value from the stationary distribution is sampled as initial value.
#' @param score (optional) a vector representing the scores (in ascending order) of the matrix index. If supplied, the result will be a vector of these values.
#' @return a Markov chain sampled from the transition matrix 
#' @details
#' It is possible to have the same score for different states of the markov chain. 
#' If no score supplied, the function returns a markov chain with state in 1:ncol(matrix)
#' @examples
#' B <-  matrix (c(0.2, 0.8, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' transmatrix2sequence(B, length = 10)
#' transmatrix2sequence(B, length = 10, score = c(-2,1))
#' transmatrix2sequence(B, length = 10, score = c("A","B"))
#' MyTransMat <-
#'    matrix(c(0.3,0.1,0.1,0.1,0.4,
#'             0.2,0.2,0.1,0.2,0.3,
#'             0.3,0.4,0.1,0.1,0.1,
#'             0.3,0.3,0.1,0.0,0.3,
#'             0.1,0.1,0.2,0.3,0.3),
#'             ncol = 5, byrow=TRUE)
#' MySeq.CM <- transmatrix2sequence(matrix = MyTransMat,length=90, score =c(-2,-1,0,2,3))
#' MySeq.CM
#' @export
transmatrix2sequence <- function(matrix, length, initialIndex, score) {
  if (!is.matrix(matrix))
    stop("[Invalid Input] transition matrix must be of type matrix.")
  if (nrow(matrix) != ncol(matrix))
    stop("[Invalid Input] transition matrix must be a square matrix.")
  if (!missing(initialIndex) && (initialIndex > nrow(matrix)))
    stop("[Invalid Input] Initial index must be within transition matrix size.")
  if (length <= 0)
    stop("[Invalid Input] Please supply a positive length for the sequence to be built.")
  if (!missing(score) && (length(score) != ncol(matrix)) )
    stop("[invalid input] the size of the transition matrix and the dimension of the score vector should match.")
  #check if initial index exists: if not, use stationary distribution to dermine the first element
  #if stationary fails, abort
  #1 Create Initial Element
  markovChain <- vector(mode = "integer", length = length)
  if (missing(initialIndex)) {
    init_distribution <-  stationary_distribution(matrix)
    first_elem <-  sample(1:length(init_distribution), size = 1, prob = init_distribution)
    markovChain[1] <- first_elem
  } else {
    markovChain[1] <- initialIndex
  }
  #2 Create Sequence with Transition Matrix
  for (i in 2:length) {
    markovChain[i] <- sample(1:ncol(matrix), size = 1, prob = matrix[markovChain[i - 1],])
  }
  #3 Map Score to Index if parameter 'score' is specified
  if (!missing(score)) {
    markov_sequence <- score[markovChain] 
    if (is.numeric(score) && all(score %% 1 == 0)) #integer values
      storage.mode(markov_sequence) <- "integer"
    return(markov_sequence)
  } else
    return(markovChain)
}

loadCharSequencesFromFile <- function(filepath){
  if (missing(filepath))
    filepath = file.choose()
  return(File2CharSequences(filepath))
}

#' @title Load score from file
#' @description Reads a csv file containing an alphabet, the associated scores, and eventually the associated probabilities.
#' @param filepath optional: location of file on disk. If not provided, a file picker dialogue box will be opened.
#' @param header does the file contain a header line (defaut to \code{TRUE})
#' @param ... optional: use arguments from read.csv
#' @details
#' The file should contains a header and 2 or 3 columns : first column the letters, second column the associated scores, optional third column associated probabilities. Letters should be unique and probabilities should sum to 1.
#' 
#' @return A data.frame. Rownames correspond to the first column, usually Letters. Associated numerical scores are in the second column. If probabilities are provided, 
#' they will be loaded too and presumed to be in the third column
#' @export
loadScoreFromFile <- function(filepath, header = TRUE, ...){
  #Pour info : do.call(rbind,dic) renvoie un dataframe classique à 1 ou 2 colonnes avec la lettre comme 'rowname'
  if (!missing(filepath) && !is.character(filepath))
    stop("[Invalid Input] Filepath must be a string.")
  if (!missing(filepath) && !file.exists(filepath))
    stop("[Invalid Input] File does not exist.")
  if (missing(filepath))
    filepath <- file.choose()
  dic <- read.csv(filepath, header = header, ...)
  #File content checks
  if (!is.numeric(dic[,2]))
    stop("This file is incompatible with the required format. Second column should contains numerical values")
  if (length(dic) > 2) {
    if (!sum(dic[,3]) == 1)
      stop("This file is incompatible with the required format. Third column should sum to 1 (probabilities)")
    if (!all(dic[,3] >= 0) && !all(dic[,3] <= 1))
      stop("This file is incompatible with the required format. Third column values should be between 0 and 1 (probabilities)")
  }
  if (anyDuplicated(dic[,1]) != 0)
    stop("This file is incompatible with the required format. First column should not contain duplicated values.")
  rownames(dic) <- dic[,1]
  dic[,1] <- NULL
  return(dic)
}

#' @title Loads matrix from csv-File
#' @description Reads a csv file without header and returns the matrix. For file formats please see section "File Formats" in vignette.
#' @param filepath optional: Location of file on disk. If not provided, a file picker dialog will be opened.
#' @return A Matrix Object
#' @export
loadMatrixFromFile <- function(filepath){
  if (!missing(filepath) && !is.character(filepath))
    stop("[Invalid Input] Filepath must be a string.")
  if (!missing(filepath) && !file.exists(filepath))
    stop("[Invalid Input] File does not exist.")
  if (missing(filepath))
    filepath <- file.choose()
  obj <- read.csv(filepath, header = FALSE)
  return(data.matrix(obj))
}

#' @title Convert a character sequence into a score sequence
#' @description Convert a character sequence into a score sequence. See CharSequences2ScoreSequences() function for several sequences
#' @param charseq a character sequence, given as a string
#' @param dictionary a data.frame with rownames containing letters, first column containing associated scores, optional second column containing associated probabilities
#' @return a vector of a score sequence
#' @seealso \code{\link{CharSequences2ScoreSequences}}
#' @examples
#' data(Seq31)
#' Seq31
#' data(HydroScore)
#' CharSequence2ScoreSequence(Seq31,HydroScore)
#' @export
CharSequence2ScoreSequence <- function(charseq, dictionary){
  string <- strsplit(charseq, "")[[1]]
  if (!is.data.frame(dictionary)) {
    if (!is.list(dictionary)) {
      stop("[Error Input] Bad dictionary format")
    } else {
      dictionary <- old2newDico(dictionary)
    }
  }
  unname(unlist(dictionary[string,1]))
}

# internal Utilitary function used to convert old format of dictionary to new format
# old format : named list of list of 1 or 2 element
# new format : dataframe with row names qual to letters, first column equal to score,
# and optionnal second column equal to probability of scores
# example : localScore:::old2newDico(Hydroscore)
old2newDico <- function(oldDico)
{
  return(as.data.frame(do.call(rbind, oldDico)))
}

#' @title Convert several character sequences into score sequences
#' @description Convert several character sequence into score sequences. For only one sequence see CharSequence2ScoreSequence() function.
#' @param sequences a list of character sequences given as string
#' @param dictionary a data.frame with rownames containing letters, first column containing associated scores, optional second column containing associated probabilities
#' @return a list of score sequences
#' @seealso \code{\link{CharSequence2ScoreSequence}}
#' @examples
#' data(Seq31)
#' Seq31
#' data(Seq219)
#' Seq219
#' data(HydroScore)
#' MySequences=list("A1"=Seq31,"A2"=Seq219)
#' CharSequences2ScoreSequences(MySequences,HydroScore)
#' @export
CharSequences2ScoreSequences <- function(sequences, dictionary){
  result <- list()
  if (is.null(names(sequences)))
    names(sequences) <- 1:length(sequences)
  for (name  in names(sequences))
    result[[name]] <- CharSequence2ScoreSequence(sequences[[name]], dictionary)
  return(result)
}
#reads file from filepath as Fasta
File2CharSequences <- function(filepath) {
  sequenceList <- list()
  add <- FALSE
  title <- ""
  con <- file(filepath, "r")
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    if (add) {
      sequenceList[[title]] <- line
      add <- FALSE
    }
    if (substr(line, 1, 1) == ">") {
      add <- TRUE
      title <- substr(line, 2, nchar(line))
      title <- strsplit(title, ' ')[[1]][[1]]
      #cut comment if any:
      
    }
  }
  close(con)
  return(sequenceList)
}

#' @title Empirical distribution from sequences
#' @description Builds empirical distribution from a list of numerical sequences
#' @details By determining the extreme scores in the sequences, this function creates a vector of probabilities
#' including values that do not occur at all. In this it differs from table(). For example, two sequences containing
#' values from 1:2 and 5:6 will produce a vector of size 6. 
#' @param sequences list of numerical sequences
#' @return empirical distribution from minimum score to maximum score as a vector of floating numbers.
#' @examples
#' seq1 = sample(7:8, size = 10, replace = TRUE)
#' seq2 = sample(2:3, size = 15, replace = TRUE)
#' l = list(seq1, seq2)
#' scoreSequences2probabilityVector(l)
#' @export
scoreSequences2probabilityVector <- function(sequences) {
  #1 concat sequences
  master_seq <- unlist(sequences)
  #2 find max and min
  max <- max(master_seq)
  min <- min(master_seq)
  result <- rep(0, length = max - min + 1)
  names(result) <- min:max
  #3 count occurence of each
  prob <- prop.table(table(master_seq))
  result[names(prob)] <- prob
  return(result)
}

#' @title Automatic analysis
#' @description Calculates local score and p-value for sequence(s) with integer scores.
#' @details This method picks the adequate p-value method for your input.\cr 
#' If no sequences are passed to this function, it will let you pick a FASTA file.\cr 
#' If this is the case, and if you haven't provided any score system 
#' (as you can do by passing a named list with the appropriate scores for each character),
#' the second file dialog which will pop up is for choosing a file containing the score 
#' (and if you provide an extra column for the probabilities, they will be used, too - see
#' section File Formats in the vignette for details).\cr 
#' The function then either uses empirical distribution based on your input - or if you provided
#' a distribution, then yours - to calculate the p-value based on the length of each of the sequences
#' given as input. \cr 
#' You can influence the choice of the method by providing the modelFunc argument. In this case, the
#' function uses exclusively simulation methods (monteCarlo, karlinMonteCarlo).  \cr 
#' By setting the method_limit you can further decide to which extent computation-intensive methods (daudin, exact_mc)
#' should be used to calculate the p-value.
#' Remark that the warnings of the localScoreC() function have be deleted when called by automatic_analysis() function
#' @param sequences sequences to be analyzed (named list)
#' @param scores vector of minimum and maximum score range
#' @param transition_matrix if the sequences are markov chains, this is their transition matrix
#' @param distribution vector of probabilities in ascending score order (iid sequences). Note that names of the vector must be the associated scores.
#' @param modelFunc function to create similar sequences. In this case, Monte Carlo is used to calculate p-value
#' @param method_limit limit length from which on computation-intensive exact calculation methods for p-value are replaced by approximative methods
#' @param score_extremes a vector with two elements: minimal score value, maximal score value
#' @param simulated_sequence_length if a modelFunc is provided and the sequence happens to be longer than method_limit, the method karlinMonteCarlo
#' is used. This method requires the length of the sequences that will be created by the modelFunc for estimation of Gumble parameters.
#' @param ... parameters for modelFunc
#' @param model the underlying model of the sequence (either "iid" for identically independently distributed variable or "markov" for Markov chains)
#' @return A list object containing
#' \item{Local score}{local score...}
#' \item{p-value}{p-value ...}
#' \item{Method}{the method used for the calculus of the p-value}
#' @examples
#' # Minimal example
#' l = list()
#' seq1 = sample(-2:1, size = 3000, replace = TRUE)
#' seq2 = sample(-3:1, size = 150, replace = TRUE)
#' l[["hello"]] = seq1
#' l[["world"]] = seq2
#' automatic_analysis(l, "iid")
#' # Example with a given distribution 
#' automatic_analysis(l,"iid",scores=-3:1,distribution=c(0.3,0.3,0.1,0.1,0.2))
#' #Equivalent call
#' automatic_analysis(l,"iid",score_extremes=c(-3,1),distribution=c(0.3,0.3,0.1,0.1,0.2))
#' # forcing the exact method for the longest sequence
#' aa1=automatic_analysis(l,"iid")
#' aa1$hello$`method applied`
#' aa1$hello$`p-value`
#' aa2=automatic_analysis(l,"iid",method_limit=3000)
#' aa2$hello$`method applied`
#' aa2$hello$`p-value`
#' # Markovian example 
#' MyTransMat <-
#' matrix(c(0.3,0.1,0.1,0.1,0.4, 0.3,0.2,0.2,0.2,0.1, 0.3,0.4,0.1,0.1,0.1, 0.3,0.3,0.3,0.0,0.1, 
#'         0.1,0.3,0.2,0.3,0.1), ncol = 5, byrow=TRUE)
#' MySeq.CM=transmatrix2sequence(matrix = MyTransMat,length=150, score =-2:2)
#' MySeq.CM2=transmatrix2sequence(matrix = MyTransMat,length=110, score =-2:2)
#' automatic_analysis(sequences = list("x1" = MySeq.CM, "x2" = MySeq.CM2), model = "markov")
#' @export
automatic_analysis <- function(sequences, model, scores, transition_matrix, distribution, method_limit = 2000, score_extremes, modelFunc, simulated_sequence_length = 1000, ...){
  #1 Check if the model of the sequences can be deduced
  if ((missing(model) && missing(transition_matrix) && missing(distribution)) || (missing(model) && !missing(transition_matrix) && !missing(distribution)))
    stop("unclear which model is used for sequence creation. Please provide either the model to be used, a transition matrix or a distribution vector.")
  #2 if model is missing, deduce it from perhaps present transition matrix
  if (missing(model)) {
    if (missing(transition_matrix))
      model <- "iid"
    else
      model <- "markov"
  }
  #3 sequences provided? If not, ask user to pick file
  if (missing(sequences)) {
    sequence_filepath <- file.choose()
    sequences <- loadCharSequencesFromFile(sequence_filepath)
  }
  #4 Sequences numeric? If yes, no scores needed, else pick file and use it
  if (!is.numeric(sequences[[1]]) && missing(scores)) {
    score_filepath <- file.choose()
    scores <- loadScoreFromFile(score_filepath)
    sequences <- CharSequences2ScoreSequences(sequences, scores)
  } 
  else if (!is.numeric(sequences[[1]]) && !missing(scores)) {
    #use scores provided in argument
    sequences <- CharSequences2ScoreSequences(sequences, scores)
  }
  
  #4.1 declare result list
  results <- list()
  #5 if extremes of scores are missing, extract them from sequences
  if (missing(score_extremes)) {
    max <- sequences[[1]][[1]]
    min <- sequences[[1]][[1]]
    for (sequence in sequences) {
      if (max(sequence) > max)
        max <- max(sequence)
      if (min(sequence) < min)
        min <- min(sequence)
    }
    score_extremes <- c(min, max)
  }
  #6 Setup Progress Bar
  if (interactive())
  progressBar <- txtProgressBar(min = 0, max = length(sequences), initial = 0, char = "=",
                               width = NA, "progress in sequence analysis", "progressbar", style = 3, file = "")
  progressCounter = 0;
  #7 model specific behaviour
  if (model == "markov") {
    #7.1 check if transition matrix available: If not, build from sequences
    if (missing(transition_matrix)) {
      tmp <- sequences2transmatrix(sequences)
      transition_matrix <- tmp$transition_matrix
      score_values <- tmp$score_values
      }
    #7.2 check if modelFunc provided
    if (missing(modelFunc)) {
      #use exact Method if sequences are not too long
      for (name in names(sequences)) {
        #calculate local score
        localscore <- localScoreC(sequences[[name]],suppressWarnings = TRUE)
        sequence_length <- length(sequences[[name]])
        if (sequence_length <= method_limit) {
          method <- "Exact Method"
          p_value <- exact_mc(localscore$localScore[1], transition_matrix, sequence_length, score_values)
#          p_value <- exact_mc(transition_matrix, localscore$localScore[1], sequence_length, sequence_min = score_extremes[1], sequence_max = score_extremes[2])
          results[[name]] <- list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        } else {
          warning(paste(c("[Ignoring Sequence] sequence length ", name, " is of length ", sequence_length), collapse = ""))
          results[[name]] <- list("localScore" = localscore)
        }
        progressCounter <- progressCounter + 1
        if (interactive())
          setTxtProgressBar(progressBar, progressCounter)
      }
    } else {
      #Function supplied, using simulation methods for Markov Chain p value determination
      for (name in names(sequences)) {
        #calculate local score
        localscore <- localScoreC(sequences[[name]],suppressWarnings = TRUE)
        sequence_length <- length(sequences[[name]])
        if (sequence_length <= method_limit) {
          method <- "Monte Carlo"
          p_value <- monteCarlo(localscore$localScore[1], modelFunc, ...)
          results[[name]] <- list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        } else {
          if (mean(sequence[[name]]) >= 0) {
            warning(paste0(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length)))
            results[[name]] <- list("localScore" = localscore)
          }
          else{
            method <- "Asymptotic Method by Monte Carlo"
            res <- karlinMonteCarlo(local_score = localscore$localScore[1], FUN = modelFunc, sequence_length = sequence_length, simulated_sequence_length = simulated_sequence_length , ...)
            results[[name]] = list("p-value" = res$p_value, "method applied" = method, "localScore" = localscore, "K*" = res$`K*`, "lambda" = res$lambda)
          }
        }
        progressCounter <- progressCounter + 1
        if (interactive())
        setTxtProgressBar(progressBar, progressCounter)
      }
    }
  } 
  else if (model == "iid") {
    #7.1 check if probability vector provided or in score (loading from file)
    if (missing(distribution)) {
      if (!missing(scores) && dim(scores)[2] > 1 ) {
        #distribution provided with score: adjust extremes if file has greater extreme values
        distribution <- scoreDictionnary2probabilityVector(scores)
        score_extremes <- c(min(as.numeric(names(distribution))), max(as.numeric(names(distribution))))
      } else {
        #learn from sequences
        distribution <- scoreSequences2probabilityVector(sequences)
      }
    } else if (!missing(scores)) {
      names(distribution) <- scores
    } else if (!missing(score_extremes)) {
      names(distribution) <- score_extremes[1]:score_extremes[2]
    } else {
      stop("[Bad Input] no score values associated to distribution. Please provide one either by parameter 'scores' or 'score_extreme'.")
    }
  
    # Moyenne des scores sur l'ensemble de la liste
    MeanDistribution <- sum(as.numeric(names(distribution))*distribution)
    
    if (missing(modelFunc)) {
      for (name in names(sequences)) {
        #calculate local score
        localscore <- localScoreC(sequences[[name]], suppressWarnings = TRUE)
        sequence_length <- length(sequences[[name]]) 
        if (sequence_length <= method_limit) {
          method <- "Exact Method Daudin et al"
          p_value <- daudin(localscore$localScore[1], sequence_length, distribution, score_extremes[1], score_extremes[2])
          results[[name]] <- list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        }
        else if (sequence_length <= method_limit*10) {
          if (MeanDistribution < 0) {
            method <- "Asymptotic Method Karlin et al"
            p_value <- karlin(localscore$localScore[1], sequence_length, distribution, score_extremes[1], score_extremes[2])
            results[[name]] <- list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
          } else {
            warning(paste0(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length)))
            results[[name]] <- list("localScore" = localscore)
          }
        } else {
          if (MeanDistribution < 0) {
            method <- "Improved Asymptotic Method Mercier et al"
            p_value <- mcc(localscore$localScore[1], sequence_length, distribution, score_extremes[1], score_extremes[2])
            results[[name]] <- list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
          } else {
            warning(paste0(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length)))
            results[[name]] <- list("localScore" = localscore)
          }
        }
        progressCounter <- progressCounter + 1
        if (interactive())
        setTxtProgressBar(progressBar, progressCounter)      
      }
    } else {
      #Function supplied, using simulation methods for iid p value determination
      for (name in names(sequences)) {
        #calculate local score
        localscore <- localScoreC(sequences[[name]],suppressWarnings = TRUE)
        sequence_length <- length(sequences[[name]])
        if (sequence_length <= method_limit) {
          method <- "Monte Carlo"
          p_value <- monteCarlo(localscore$localScore[1], FUN = modelFunc, ...)
          results[[name]] <- list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        } else {
          if (MeanDistribution >= 0) {
            warning(paste(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length), collapse = ""))
            results[[name]] <- list("localScore" = localscore)
          } else {
            method <- "Asymptotic Method by Monte Carlo"
            res <- karlinMonteCarlo(local_score = localscore$localScore[1], FUN = modelFunc, sequence_length = sequence_length, ...)
            results[[name]] <- list("p-value" = res$p_value, "method applied" = method, "localScore" = localscore, "K*" = res$`K*`, "lambda" = res$lambda)
          }
        }
        progressCounter = progressCounter + 1
        if (interactive())
        setTxtProgressBar(progressBar, progressCounter)
      }
    }
  } 
  else {
    stop("unknown model. Please mention either iid or markov.")
  }
  if (interactive())
  close(progressBar)
  return(results)
}



#' @title Check for missing scores values in the score distribution
#' @description Get extremes scores, then create a complete list and set a probability equal to zeros for not present scores
#' @param listScore list of vector containing each (score value, probability)
#' @param score_extremes (optional) vector of min and max value of scores (used for truncation).
#' @return vector containing all score values between extremes and the probability equal to 0 for missing score
#' @examples
#' Mylist <- list("x1"=c(-2,0.1),"x2"=c(0,0.7),"x3"=c(1,0.2))
#' scoreDictionnary2probabilityVector(list = Mylist, score_extremes = c(-2,1))
#' @noRd
scoreDictionnary2probabilityVector <- function(listScore, score_extremes = NULL){
  score <- do.call(rbind, listScore)
  probs <- score[,2]
  names(probs) <- score[,1]

 #check for missing values and set zeros if not present
 #get extremes, then create list and check difference
 missingScores <- setdiff(min(score[,1]):max(score[,1]), score[,1])
 missingProbs <- rep(0,length(missingScores))
 names(missingProbs) <- missingScores
 probs <- c(probs, missingProbs)
 #order from lowest to highest score
 probs <- probs[order(as.numeric(names(probs)))]
 #si la fonction est appele avec des valeurs "limitantes" on coupe ce qui depasse ces valeurs
 if (!is.null(score_extremes)) {
   probs <- subset(probs, as.numeric(names(probs)) >= score_extremes[1] & as.numeric(names(probs)) <= score_extremes[2])
 }
 return(probs)
}

#' @title Convert a real scores vector into an integer scores vector
#' @description Convert real scores into integer scores
#' @details Convert real scores into integer scores by multiplying real
#' scores by a coefficient (default 10) and then assigning probability to corresponding
#' extended (from the minimum to the maximum) integer scores
#' @param RealScore vector of real scores
#' @param ProbRealScore vector of probability
#' @param coef coefficient
#' @return list containing ExtendedIntegerScore and ProbExtendedIntegerScore
#' @examples
#' score <- c(-1,-0.5,0,0.5,1)
#' prob.score <- c(0.2,0,0.4,0.1,0.3)
#' (res1 <- RealScores2IntegerScores(score, prob.score, coef=10))
#' prob.score.err <- c(0.1,0,0.4,0.1,0.3)
#' (res2 <- RealScores2IntegerScores(score, prob.score.err, coef=10))
#' # When coef=1, the function can handle integer scores
#' ex.integer.score <- c(-3,-1,0,1, 5)
#' (res3 <- RealScores2IntegerScores(ex.integer.score, prob.score, coef=1))
#' @export
RealScores2IntegerScores <- function(RealScore, ProbRealScore, coef = 10)
{
  if (sum(ProbRealScore) != 1) warning("Sum of probabilities is different from 1")
  
  IntegerScore <- RealScore * coef
  
  MinIntegerScore <- min(IntegerScore)
  MaxIntegerScore <- max(IntegerScore)
  
  ExtendedIntegerScore <- as.integer(MinIntegerScore:MaxIntegerScore)
  
  ProbExtendedIntegerScore <- rep(0, times = length(ExtendedIntegerScore))
  
  ProbExtendedIntegerScore[ExtendedIntegerScore %in% IntegerScore] <- ProbRealScore
  names(ProbExtendedIntegerScore) <- ExtendedIntegerScore
  
  list(ExtendedIntegerScore = ExtendedIntegerScore,
       ProbExtendedIntegerScore = ProbExtendedIntegerScore)
}

#' @title Tool function used in calculus of probability of excursion (Markov case). Return the vector z
#'
#' @param a score strictly positive
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix of size, theta x theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#'
#' @return Return the vector z
#' @examples
#' theta = c("a","b","c","d")
#' lambda = matrix(c(0.1,0.2,0.4,0.3,
#'                   0.3,0.1,0.5,0.1,
#'                   0.2,0.6,0.1,0.1,
#'                   0.4,0.1,0.1,0.4),
#'                 ncol=4, byrow = TRUE, dimnames = list(theta,theta))
#' all(apply(lambda,1,sum)==1) #TRUE for markov transition matrix
#' stationary_dist = stationary_distribution(lambda))
#' score_function = c(a=-3,b=-1,c=1,d=2)
#' sum(score_function*stationary_dist) #Score expectation (should be <0)
#' proba_excursion_markov_calcul_z(3, theta, lambda, score_function)
#' @noRd
proba_excursion_markov_calcul_z <- function(a, theta, lambda, score_function) {
  # RM théoriquement plus utilisée (le calcul de z est effectué directement dans proba_theoretical_first_excursion_alpha)
  list.theta <- as.list(theta)
  list.a <- as.list(c(-1:a))
  combination <- as.vector(outer(list.theta,list.a, paste, sep = ":"))
  
  P <- create_matrix_P(a, theta, lambda, score_function, combination)
  
  Q_C <- P[,-(length(combination) - length(theta) + 1):-(length(combination))]
  Q_C <- Q_C[,-1:-length(theta)]
  Q_B <- P[,1:length(theta)]                                            
  Q_A <- P[,(length(combination) - length(theta) + 1):length(combination)]
  
  return(z = (solve(diag(ncol(Q_C)) - Q_C)) %*% (Q_A %*% rep(1, ncol(Q_A))))
}

#' @title Used to generate markov's sequences (can output sequence of cumulative scores or simple sequences of scores)
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix of size, theta x theta
#' @param prob0 stationary distribution of lambda 
#' @param n length of generated Markov's sequence
#' @return vector containing the generated Markov's score sequence
#' @examples
#' hist_S = generate_markov_hist_score(theta, lambda, score_function, prob0, n)
#' @noRd
generate_markov_hist_score <- function(theta, lambda, score_function, prob0, n){
  card <- length(theta)
  previous_index <- sample(card, size = 1, replace = TRUE, prob = prob0)
  hist_score <- rep(NA, length = (n + 1))
  hist_score[1] <- 0
  hist_score[2] <- score_function[previous_index]
  
  for (x in 3:n ) {
    previous_index <- sample(1:card, size = 1, replace = TRUE, prob = lambda[previous_index,])
    hist_score[x] <- hist_score[x - 1] + score_function[previous_index]
  }
  
  return(hist_score)
}

#' @title Used generate a markov score sequence and extract its first excursion
#' @param a score
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix of size, theta x theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param prob0 stationary distribution of lambda 
#' @param n length of generated Markov's chain
#' @return list containing three elements : hist_S      : cumulative score sequence
#'                                          hist_S_star : cumulative score sequence filled after the first excurtion with -1 if the score is not reached andor else with the max reached score 
#' @examples
#' hist_s_star = generation_empirical_processes(a, theta, lambda, score_function, prob0, n)$hist_S_star
#' @noRd
generation_empirical_processes <- function(a, theta, lambda, score_function, prob0 = NULL, n) {
  if (is.null(prob0)){
    prob0 = stationary_distribution(lambda)
  }
  
  hist_S = generate_markov_hist_score(theta, lambda, score_function, prob0, n)
  
  hist_S_star = rep(NA, length=(n))
  hist_S_star[1] = 0
  
  for (i in 2:(length(hist_S))) {
    #print(hist_S[i])
    if (0 <= hist_S[(i-1)] & hist_S[(i-1)] < a) {
      hist_S_star[(i-1)] = hist_S[(i-1)]
    } else if ( hist_S[(i-1)] >= a ) {
      hist_S_star[(i-1):n] = rep(a, length = ( (n+1) - (i-1)) )
      return(list(hist_S = hist_S, hist_S_star = hist_S_star))
    } else if ( hist_S[(i-1)] < 0  ) {
      hist_S_star[(i-1):n] = rep(-1, length = ( (n+1) - (i-1)) )
      return(list(hist_S = hist_S, hist_S_star = hist_S_star))
    }
  }
  return(list(hist_S = hist_S, hist_S_star = hist_S_star))
}

#' @title return the empirical probability of reaching a score a on the first excursion of a markov's score sequence
#' @param a score
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix of size, theta x theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param prob0 stationary distribution of lambda
#' @param n length of Markov's chain
#' @param N number of repetition of the process
#' @return reel number between 0 and 1 representing the empirical probability of reaching a score of a
#'  on the first excursion 
#' @examples
#' probability = localScore::proba_empirical_first_excursion(5, c("K","L","M"), matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), c(-2,-1,2), c(0.4444444, 0.2962963, 0.2592593), n = 1000, N = 10000)
#' @noRd
proba_empirical_first_excursion_markov <- function(a, theta, lambda, score_function, prob0, n, N) {
  s = 0
  for (x in 1:N) {
    hist_s_star = generation_empirical_processes(a, theta, lambda, score_function, prob0, n)$hist_S_star
    if (hist_s_star[n] == a) {
      s = s + 1
    }
  }
  return(s/N)
}

#' @title return the empirical probability of reaching a score a on the ith (with i >= 2) excursion of a markov's score sequence
#' @param a score
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix of size, theta x theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param prob0 stationary distribution of lambda
#' @param n length of Markov's chain
#' @param N number of repetition of the process
#' @param i wanted excursion
#' @return reel number between 0 and 1 representing the empirical probability of reaching a score of a on the ith excursion 
#' @examples
#' probability = localScore::proba_empirical_ith_excursion(c("K","L","M"), 5, matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), c(-2,-1,2), c(0.4444444, 0.2962963, 0.2592593), i = 2, n = 1000, N = 10000)
#' @noRd
proba_empirical_ith_excursion_markov <- function(a, theta, lambda, score_function, prob0, n, N, i) {
  c = 0
  for ( it in 1:(N+1) ) {
    partial_sum = generation_empirical_processes(a, theta, lambda, score_function, prob0, n)
    partial_sum = partial_sum$hist_S
    index = 0
    exc   = 1
    val_prev_exc = 0
    
    while ( ( exc < i ) & index < ( n - 1 ) ) {
      if ( partial_sum[index+1] >= val_prev_exc ) {
        index = index + 1
      } else {
        exc = exc + 1
        val_prev_exc = partial_sum[index+1]
      }
    }
    
    while ( partial_sum[index+1] >= val_prev_exc & index < ( n - 1 ) & partial_sum[index+1] < ( a + val_prev_exc ) ) {
      index = index + 1
    }
    
    if ( partial_sum[index+1] >= ( a + val_prev_exc ) & exc == i ) {
      c = c + 1
    }
    
  }
  return (c/N)
}

##########################################################
