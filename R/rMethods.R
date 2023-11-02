#' @useDynLib localScore
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline barplot par plot
#' @importFrom stats lm ecdf knots
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar
#' 
#' 
#' @title Monte Carlo method [p-value]
#' @description Calculates an empirical p-value based on simulations of similar integer sequences of the same length.
#' Perfect for small sequences (both markov chains and identically and independantly distributed) with length ~ 10^3. 
#' See function monteCarlo_double() for possible real scores.
#' @param local_score local score observed in a segment.
#' @param FUN function to simulate similar sequences with.
#' @param ... parameters for FUN
#' @param plot boolean value if to display plots for cumulated function and density
#' @param numSim number of sequences to generate during simulation
#' @return Floating value corresponding to the probability to obtain a local score with value greater or equal to the parameter
#' @examples
#' \donttest{
#' monteCarlo(120, FUN = rbinom, n = 100, size = 5, prob=0.2)
#' }
#' new = sample(-7:6, replace = TRUE, size = 1000)
#' #MonteCarlo taking random sample from the input sequence itself
#' \donttest{
#' monteCarlo(local_score = 20, FUN = function(x) {return(sample(x = x, 
#' size = length(x), replace = TRUE))}, x=new)
#' }
#' # Markovian example
#' MyTransMat <-
#' +     matrix(c(0.3,0.1,0.1,0.1,0.4, 0.2,0.2,0.1,0.2,0.3, 0.3,0.4,0.1,0.1,0.1, 0.3,0.3,0.1,0.0,0.3,
#' +              0.1,0.1,0.2,0.3,0.3), ncol = 5, byrow=TRUE)
#' \donttest{
#' monteCarlo(local_score = 50,
#'           FUN = transmatrix2sequence, matrix = MyTransMat,
#'           length=150, score = c(-2,-1,0,2,3), plot=FALSE, numSim = 5000)
#' }
#' @export
monteCarlo <- function(local_score, FUN, ..., plot = TRUE, numSim = 1000){
  epsilon <- 1e-8 #Awfull trick used to correct CdF in cas of discrete score values
  if (local_score < 0)
    stop("[Invalid Input] local scores are never negative.")
  if (numSim <= 0)
    stop("[Invalid Input] Please supply a positive number of simulations.")
  #simulate sequences to get some local scores to work with
  maxScores <- sim_local_scores(numSim = numSim, func = FUN, ...)
  #empirical distribution of found local Scores
  scoreCdf <- ecdf(maxScores)
  #Compute Density
  scoreD = NULL
  smin <- min(knots(scoreCdf))
  smax <- max(knots(scoreCdf))
  for (k in smin:smax) { 
    scoreD[k - smin + 1] <- scoreCdf(k) - scoreCdf(k - 1)
  }
  names(scoreD) = smin:smax
  #if the user wants to display the distribution of local scores, print it to the plot screen
  if (plot) {
    #set plot window
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    par(mfrow = c(1,2))
    colors <- rep("gray", smax - smin + 1)
    if (local_score >= smin & local_score <= smax)
      colors[as.character(local_score)] = "red"
    barplot(scoreD, main = "Distribution of local scores [iid]\n for given sequence", ylab = "Probability", xlab = "local Scores", col = colors, border = "white", space = 0.1)
    plot(scoreCdf, main = "Cumulative Distribution Function", ylab = "cumulative destribution(local_scores)))", xlab = "Local Score Values")
    abline(v = local_score, col = "red",lty = 3)  #greater or equal, thus set marker at localScore -1
    abline(h = scoreCdf(local_score), col = "red",lty = 3)
  }
  result <- 1 - scoreCdf(local_score * (1 - epsilon))
  names(result) <- c("p-value")
  return(result)
}

sim_local_scores <- function(func, ..., numSim = 1000) {
  # creates a vector of numSim localScores of numSim simulated sequences by the function passed
  maxScores <- vector(mode = "numeric", length = numSim)

  for (i in 1:numSim) {
    simSequence <- func(...)
    maxScores[i] <- localScoreC(simSequence, supressWarnings = TRUE)$localScore["value"]
  }
  return(maxScores)
}


#' @title Monte Carlo method for real score case [p-value]
#' @description Calculates an empirical p-value based on simulations of similar sequences of the same length.
#' Perfect for small sequences (both markov chains and identically and independantly distributed) with length ~ 10^3. Function dedicated for real score case.
#' @param local_score local score observed in a segment.
#' @param FUN function to simulate similar sequences with.
#' @param ... parameters for FUN
#' @param plot Boolean value if to display plots for cumulated function and density
#' @param numSim number of sequences to generate during simulation
#' @return Floating value corresponding to the probability to obtain a local score with value greater or equal to the parameter
#' @examples
#'
#' score_reels=c(-1,-0.5,0,0.5,1)
#' proba_score_reels=c(0.2,0.3,0.1,0.2,0.2)
#' sample_from_model <- function(score.sple,proba.sple, length.sple){sample(score.sple,
#'                                   size=length.sple, prob=proba.sple, replace=TRUE)}
#' \donttest{
#' monteCarlo_double(5.5,FUN=sample_from_model, plot = TRUE, 
#' score.sple=score_reels,proba.sple=proba_score_reels, length.sple=100, numSim = 1000)
#' }
#' @export
monteCarlo_double <- function (local_score, FUN, ..., plot = TRUE, numSim = 1000){
  epsilon = 1e-8 #Awfull trick used to correct CdF in case of discrete score values
  if(local_score<0)
    stop("[Invalid Input] local Scores are newer negative.")
  if(numSim<=0)
    stop("[Invalid Input] Please supply a positive number of simulations.")
  #simulate sequences to get some local scores to work with
  maxScores = sim_local_scores_double(numSim = numSim, func = FUN, ...)
  #empirical distribution of found local Scores
  scoreCdf = ecdf(maxScores)
  #Compute Density
  scoreD = NULL
  smin = min(knots(scoreCdf))
  smax = max(knots(scoreCdf))
  for (k in smin:smax) { 
    scoreD[k-smin+1]=scoreCdf(k)-scoreCdf(k-1)
  }
  names(scoreD) = smin:smax
  #if the user wants to display the distribution of local scores, print it to the plot screen
  if(plot){
    #set plot window
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    par(mfrow =(c(1,2)))
    colors <- rep("gray", smax-smin+1)
    if(local_score >= smin & local_score<=smax)
      colors[as.character(local_score)] = "red"
    barplot(scoreD, main = "Distribution of local scores [iid]\n for given sequence", ylab = "Probability", xlab = "local Scores", col=colors, border = "white", space = 0.1)
    plot(scoreCdf, main="Cumulative Distribution Function", ylab="cumulative destribution(local_scores)))", xlab="Local Score Values")
    abline(v=local_score, col="red",lty=3)  #greater or equal, thus set marker at localScore -1
    abline(h=scoreCdf(local_score), col="red",lty=3)
  }
  result = 1-scoreCdf(local_score*(1-epsilon))
  names(result) = c("p-value")
  return (result)
}

sim_local_scores_double <- function (func, ..., numSim = 1000){
  #creates a vector of numSim localScores of numSim simulated sequences by the function passed
  maxScores = vector(mode="numeric", length=numSim)
  
  for(i in 1:numSim){
    simSequence = func(...)
    maxScores[i] <- localScoreC_double(simSequence,supressWarnings = TRUE)$localScore["value"]
  }
  return (maxScores)
}


#' @title Lindley process
#' @description Creates a sequence of a Lindley process, also called CUSUM process, on a given sequence.  For a sequence (X_k)k, the Lindley process is defined as follows: W_0:=0 and W_(k+1)=max(0,W_k+X_(k+1)). It defines positive excursions above 0.
#' @param sequence numeric sequence of a Lindley process, eg service time per customer
#' @return a vector with the Lindley process steps
#' @examples
#' seq=c(1,2,3,-4,1,-3,-1,2,3,-4,1)
#' lindley(seq)
#' plot(1:length(seq),lindley(seq),type='b')
#' @export
lindley <- function(sequence){
  seq = rep(NA,length(sequence))
  i=1
  sum = 0
  for( score in sequence){
    if(sum+score<=0)
      sum=0
    else
      sum=score+sum
    seq[i]= sum
    i=i+1
  }
  return (seq)
}

#' @title Monte Carlo - Karlin [p-value]
#' @description Estimates p-value, for integer scores, based on a Monte Carlo estimation of Gumble
#' parameters from simulations of smaller sequences with same distribution. 
#' Appropriate for great sequences with length > 10^3, for i.i.d and markovian sequence models.
#' @details The length of the simulated sequences is an argument specific to the function provided for simulation. Thus, it has to be
#' provided also in the parameter simulated_sequence_length in the arguments of the "Monte Carlo - Karlin" function. 
#' It is a crucial detail as it influences precision and computation time of the result. Note that to get an appropriate 
#' estimation, the given average score must be non-positive.
#' @param local_score local score observed in a segment.
#' @param FUN function to simulate similar sequences with.
#' @param ... parameters for FUN
#' @param plot boolean value if to display plots for cumulated function and density
#' @param numSim number of sequences to create for estimation
#' @param sequence_length length of the sequence
#' @param simulated_sequence_length length of simulated sequences produced by FUN 
#' @return Floating value corresponding to the probability to obtain a local score with a value greater or equal to the parameter local_score
#' @examples
#' new = sample(-7:6, replace = TRUE, size = 1000) 
#' #MonteCarlo taking random sample from the input sequence itself
#' \donttest{
#' karlinMonteCarlo(local_score = 66, sequence_length = 1000,  
#'                FUN = function(x, simulated_sequence_length) {return(sample(x = x, 
#'                size = simulated_sequence_length, replace = TRUE))}, 
#'                x=new, simulated_sequence_length = 1000,  numSim = 1000)
#' }
#'                
#' @export
karlinMonteCarlo <- function (local_score, sequence_length, simulated_sequence_length, FUN,  ..., numSim = 1000, plot = TRUE){
  if(local_score<0)
    stop("[Invalid Input] local Scores are not negative.")
  if(numSim<=0)
    stop("[Invalid Input] Please supply a positive number of simulations.")
  if(missing(sequence_length))
    stop("[Missing Argument] sequence_length is a required parameter")
  if(missing(simulated_sequence_length))
    stop("[Missing Argument] simulated_sequence_length is a required parameter")
  #simulate sequences to get some local scores to work with
  maxScores = sim_local_scores(numSim = numSim, func = FUN, ...)
  
  #Cumulative distribution function of local score probabilities
  log_cumsum = log(-log(cumsum(table(maxScores)/length(maxScores))))
  
  #clean up eventual inf values
  log_cumsum <- log_cumsum[is.finite(log_cumsum)]
  
  #linear regression (Warning : at least two different simulated values of MaxScores needed)
  if (length(log_cumsum) < 2) {
    p_value = NA
    k_star = NA
    lambda = NA
  } else {
    reg <- lm(log_cumsum~as.numeric(names(log_cumsum)))
    lambda <- -reg$coefficients[2]      #Inclination
    names(lambda) <- NULL
    k_star <- exp(reg$coefficients[1])/simulated_sequence_length #SIMULATED SEQUENCE LENGTH
    names(k_star) <- NULL
    x <- local_score - log(sequence_length)/lambda
    #Karlin: for n great, P( ln(n)/lambda+x>= M) = exp(-K_star*exp(-lambda*x))
    #thus we set ln(n)/lambda+x = local_score and obtain for x = local_score - ln(n)/lambda
    #x = local_score - log(sequence_length)/lambda#-1 ? XXX a priori oui, il faut le '-1' (cf le code C++)
    
    #now we calculate p-value with our approximate K star and lambda
    p_value <- exp(-k_star*exp(-lambda*x))
  }
  result <- list("p-value" = 1 - p_value, "K*" = as.numeric(k_star), "lambda" = as.numeric(lambda))
  if (plot){
    #set plot window to contain 3 slots
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    par(mfrow =(c(1,3)))
    colors <- rep("gray", length(unique(maxScores)))
    if(local_score %in% as.numeric(names(table(maxScores))))
      colors[match(local_score, as.numeric(names(table(maxScores)))):max(maxScores)] = "red"
    barplot(table(maxScores), main = "Distribution of local scores [iid]\n for given sequence", ylab = "Number of occurences", xlab = "Local score values", col=colors, border = "white", space = 0.1)
    plot(cumsum(table(maxScores)/length(maxScores))~as.numeric(names(table(maxScores))), main="Cumulative Distribution Function", ylab="cumulative distribution(local_scores)))", xlab="Local Score Values")
    abline(v=local_score-1, col="red")
    plot(log_cumsum~as.numeric(names(log_cumsum)), main="Linear Regression", ylab="ln(-ln(cumulative destribution(local_scores)))", xlab="Local score values")
    abline(reg, col=c("green"))
  }
  return (result)
}


#' @title Monte Carlo - Karlin for real scores[p-value]
#' @description Estimates p-value, for integer scores, based on a Monte Carlo estimation of Gumble
#' parameters from simulations of smaller sequences with same distribution. 
#' Appropriate for great sequences with length > 10^3, for i.i.d and markovian sequence models.
#' @details The length of the simulated sequences is an argument specific to the function provided for simulation. Thus, it has to be
#' provided also in the parameter simulated_sequence_length in the arguments of the "Monte Carlo - Karlin" function. 
#' It is a crucial detail as it influences precision and computation time of the result. Note that to get an appropriate 
#' estimation, the given average score must be non-positive.
#' @param local_score local score observed in a segment.
#' @param FUN function to simulate similar sequences with.
#' @param ... parameters for FUN
#' @param plot boolean value if to display plots for cumulated function and density
#' @param numSim number of sequences to create for estimation
#' @param sequence_length length of the sequence
#' @param simulated_sequence_length length of simulated sequences produced by FUN 
#' @return Floating value corresponding to the probability to obtain a local score with value greater or equal the parameter
#' @examples
#' new = sample(-7:6, replace = TRUE, size = 1000) 
#' #MonteCarlo taking random sample from the input sequence itself
#' \donttest{
#' karlinMonteCarlo_double(local_score = 66, sequence_length = 1000,  
#'                FUN = function(x, simulated_sequence_length) {return(sample(x = x, 
#'                size = simulated_sequence_length, replace = TRUE))}, 
#'                x=new, simulated_sequence_length = 1000,  numSim = 1000)
#'  }     
#'  #Markovian example (longer computation)
#'  MyTransMat_reels <-  matrix(c(0.3,0.1,0.1,0.1,0.4, 0.2,0.2,0.1,0.2,0.3, 0.3,0.4,0.1,0.1,0.1, 
#' 0.3,0.3,0.1,0.2,0.1, 0.2,0.1,0.2,0.4,0.1),ncol = 5, byrow=TRUE)
#' \donttest{
#' karlinMonteCarlo(local_score = 10.5,sequence_length=200,FUN = transmatrix2sequence, 
#' matrix = MyTransMat_reels, score =c(-1,-0.5,0,0.5,1),length=1500, 
#' plot=FALSE, numSim = 1500, simulated_sequence_length =1500)
#' }
#' @export
karlinMonteCarlo_double <- function (local_score, sequence_length, simulated_sequence_length, FUN,  ..., numSim = 1000, plot = TRUE){
  if(local_score<0)
    stop("[Invalid Input] local scores are not negative.")
  if(numSim<=0)
    stop("[Invalid Input] Please supply a positive number of simulations.")
  if(missing(sequence_length))
    stop("[Missing Argument] sequence_length is a required parameter")
  if(missing(simulated_sequence_length))
    stop("[Missing Argument] simulated_sequence_length is a required parameter")
  #simulate sequences with real scores to get some local scores to work with
  maxScores = sim_local_scores_double(numSim = numSim, func = FUN, ...)
  
  #Cumulative distribution function of local score probabilities
  log_cumsum = log(-log(cumsum(table(maxScores)/length(maxScores))))
  
  #clean up eventual inf values
  log_cumsum <- log_cumsum[is.finite(log_cumsum)]
  
  #linear regression
  reg <- lm(log_cumsum~as.numeric(names(log_cumsum)))
  lambda <- -reg$coefficients[2]      #Inclination
  names(lambda) <- NULL
  k_star <- exp(reg$coefficients[1])/simulated_sequence_length #SIMULATED SEQUENCE LENGTH
  names(k_star) <- NULL
  x <- local_score - log(sequence_length)/lambda
  #Karlin: for n great, P( ln(n)/lambda+x>= M) = exp(-K_star*exp(-lambda*x))
  #thus we set ln(n)/lambda+x = local_score and obtain for x = local_score - ln(n)/lambda
  #x = local_score - log(sequence_length)/lambda#-1 ?
  
  #now we calculate p-value with our approximate K star and lambda
  p_value <- exp(-k_star*exp(-lambda*x))
  result <- list("p-value" = 1 - p_value, "K*" = as.numeric(k_star), "lambda" = as.numeric(lambda))
  if(plot){
    #set plot window to contain 3 slots
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    par(mfrow =(c(1,3)))
    colors <- rep("gray", length(unique(maxScores)))
    if(local_score %in% as.numeric(names(table(maxScores))))
      colors[match(local_score, as.numeric(names(table(maxScores)))):max(maxScores)] = "red"
    barplot(table(maxScores), main = "Distribution of local scores [iid]\n for given sequence", ylab = "Number of Occurences", xlab = "Local score values", col=colors, border = "white", space = 0.1)
    plot(cumsum(table(maxScores)/length(maxScores))~as.numeric(names(table(maxScores))), main="Cumulative Distribution Function", ylab="cumulative destribution(local_scores)))", xlab="Local Score Values")
    abline(v=local_score-1, col="red")
    plot(log_cumsum~as.numeric(names(log_cumsum)), main="Linear Regression", ylab="ln(-ln(cumulative destribution(local_scores)))", xlab="Local score values")
    abline(reg, col=c("green"))
  }
  return (result)
}

#' @title Transition matrix from sequence(s)
#' @description Calculates the transition matrix by counting occurences of tuples in given vector list
#' @param sequences Sequences to be analyzed, as list of numeric vectors
#' @return A list object containing
#' \item{Transition Matrix}{Transition Matrix}
#' \item{Min Value}{minimal score found in supplied sequences}
#' \item{Max Value}{maximal score found in supplied sequences}  
#' @details  The transition matrix will be structured so that the lowest score corresponds to the first column and row
#' and the highest score corresponds to the last column and row. Note that the resulting matrix is not stochastic because it can occur rows filled up with only 0 for not observed score in Min Value Max Value interval.
#' @examples
#' seq = sample(-1:1, size = 20, replace = TRUE)
#' seq2 = sample(-6:1, size = 20, replace = TRUE)
#' seq3 = sample(3:6, size = 50, replace = TRUE)
#' sequences2transmatrix(list(seq, seq2, seq3))
#' @export
sequences2transmatrix <- function(sequences){
  if(!is.list(sequences))
    stop("[Invalid Input] Sequences must be numeric vectors in a list.")
  #1 get score range
  unique_values = unique(sequences[[1]])
  score = c(min(unique_values), max(unique_values))
  if (length(sequences)>1) {
  for(x in 2:length(sequences)){
    unique_values = unique(sequences[[x]])
    if(min(unique_values)<score[1])
      score[1] = min(unique_values)
    if(max(unique_values)>score[2])
      score[2] = max(unique_values)
  }
  }
  #2 create matrix
  p <- matrix(nrow = score[2]-score[1]+1, ncol = score[2]-score[1]+1, 0)
  #3 count tuples
  for (sequence in sequences){
    for (t in 1:(length(sequence) - 1)) 
      p[sequence[t]-score[1]+1, sequence[t + 1]-score[1]+1] <- p[sequence[t]-score[1]+1, sequence[t + 1]-score[1]+1] + 1  
  }
  #4 normalize rows = 1
  for (i in 1:nrow(p)){ 
    if(sum(p[i, ]>0))
      p[i, ] <- p[i, ] / sum(p[i, ])
  }
  #5 return result
  return(list("Transition Matrix" = p, "Min Value" = score[1], "Max Value" = score[2]))
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
#' @examples
#' B = t(matrix (c(0.2, 0.8, 0.4, 0.6), nrow = 2))
#' transmatrix2sequence(B, length = 10)
#' MyTransMat <-
#'    matrix(c(0.3,0.1,0.1,0.1,0.4, 0.2,0.2,0.1,0.2,0.3, 0.3,0.4,0.1,0.1,0.1, 0.3,0.3,0.1,0.0,0.3,
#'     0.1,0.1,0.2,0.3,0.3), ncol = 5, byrow=TRUE)
#' MySeq.CM=transmatrix2sequence(matrix = MyTransMat,length=90, score =c(-2,-1,0,2,3))
#' MySeq.CM
#' @export
transmatrix2sequence <- function(matrix, length, initialIndex, score){
  if(!is.matrix(matrix))
    stop("[Invalid Input] transition matrix must be of type matrix.")
  if(nrow(matrix)!=ncol(matrix))
    stop("[Invalid Input] transition matrix must be a square matrix.")
  if(!missing(initialIndex) && initialIndex>nrow(matrix))
    stop("[Invalid Input] Initial index must be within transition matrix size.")
  if(length<=0)
    stop("[Invalid Input] Please supply a positive length for the sequence to be built.")
  #check if initial index exists: if not, use stationary distribution to dermine the first element
  #if stationary fails, abort
  #1 Create Initial Element
  markovChain = c()
  if(missing(initialIndex)){
    init_distribution = stationary_distribution(matrix)
    first_elem = sample(1:length(init_distribution), size = 1, prob = init_distribution)
    markovChain <- append(markovChain, first_elem)
  } else {
    markovChain <- append(markovChain, initialIndex)
  }
  #2 Create Sequence with Transition Matrix
  for(i in 2:length){
    markovChain <- append(markovChain, sample(1:ncol(matrix), size = 1, prob = matrix[markovChain[i-1],]))
  }
  #3 Map Score to Index if necessary
  if(!missing(score)){
    markov_sequence = c()
    for(x in markovChain)
      markov_sequence <- append(markov_sequence, score[x])
    return (markov_sequence)
  } else
    return(markovChain)
}

loadCharSequencesFromFile <- function(filepath){
  if(missing(filepath))
    filepath = file.choose()
  return (File2CharSequences(filepath))
}

#' @title Load score from file
#' @description Reads a csv file with 2-3 columns and returns it as a list object of vectors, with names corresponding to the first column of the file. For details view section "File Formats" in vignette.
#' @param filepath optional: location of file on disk. If not provided, a file picker dialog will be opened.
#' @param ... optional: use arguments from read.csv
#' @return A List Object - Names correspond to the first column, usually Letters. Associated numerical scores are in the second column. If probabilities are provided, 
#' they will be loaded too and presumed to be in the third column
#' @export
loadScoreFromFile <- function(filepath, ...){
  if(!missing(filepath) && !is.character(filepath))
    stop("[Invalid Input] Filepath must be a string.")
  if(!missing(filepath) && !file.exists(filepath))
    stop("[Invalid Input] File does not exist.")
  if(missing(filepath))
    filepath = file.choose()
  obj = read.csv(filepath, ...)
  #File content checks
  if (!is.numeric(obj[,2]))
    stop("This file is incompatible with the required format. Second column should contains numerical values")
  if (length(obj) > 2) {
    if (!sum(obj[,3]) == 1)
      stop("This file is incompatible with the required format. Third column should sum to 1 (probabilities)")
    if (!all(obj[,3] >= 0) && !all(obj[,3] <= 1))
      stop("This file is incompatible with the required format. Third column values should be between 0 and 1 (probabilities)")
  }
    
  dic = list()
  if(length(obj) > 2){
    for(i in 1:nrow(obj))
      dic[[as.character(obj[i,1])]] <- list(obj[i,2], obj[i,3])
  }
  else if(length(obj) == 2){
    for(i in 1:nrow(obj))
      dic[[as.character(obj[i,1])]]=list(obj[i,2])
  }  
  else
    warning("This file is incompatible with the required format. Please check the documentation, chapter FILE FORMATS")
  return(dic)
}

#' @title Loads matrix from csv-File
#' @description Reads a csv file without header and returns the matrix. For file formats please see section "File Formats" in vignette.
#' @param filepath optional: Location of file on disk. If not provided, a file picker dialog will be opened.
#' @return A Matrix Object
#' @export
loadMatrixFromFile <- function(filepath){
  if(!missing(filepath) && !is.character(filepath))
    stop("[Invalid Input] Filepath must be a string.")
  if(!missing(filepath) && !file.exists(filepath))
    stop("[Invalid Input] File does not exist.")
  if(missing(filepath))
    filepath = file.choose()
  obj = read.csv(filepath, header = FALSE)
  return (data.matrix(obj))
}

#' @title Convert a character sequence into a score sequence
#' @description Convert a character sequence into a score sequence. See CharSequences2ScoreSequences() fonction for several sequences
#' @param sequence a character sequence
#' @param dictionary a dictionary
#' @return a vector of a score sequence
#' @examples
#' data(ShortSeq)
#' ShortSeq
#' data(dico)
#' CharSequence2ScoreSequence(ShortSeq,dico)
#' @export
CharSequence2ScoreSequence <- function(sequence, dictionary){
  scoresequence = c()
  string <- strsplit(sequence, "")[[1]]
  for(char in string)
    scoresequence = append(scoresequence, dictionary[[char]][[1]])
  return (scoresequence)
}

#' @title Convert several character sequences into score sequences
#' @description Convert several character sequence into score sequences. For only one sequence see CharSequence2ScoreSequence() function.
#' @param sequences a list of character sequences
#' @param dictionary a dictionary
#' @return a list of score sequences
#' @examples
#' data(ShortSeq)
#' ShortSeq
#' data(MidSeq)
#' MidSeq
#' data(dico)
#' MySequences=list("A1"=ShortSeq,"A2"=MidSeq)
#' CharSequences2ScoreSequences(MySequences,dico)
#' @export
CharSequences2ScoreSequences <- function(sequences, dictionary){
  result = list()
  for(name  in names(sequences))
    result[[name]] = CharSequence2ScoreSequence(sequences[[name]], dictionary)
  return (result)
}
#reads file from filepath as Fasta
File2CharSequences = function(filepath) {
  sequenceList = list()
  add = FALSE
  title = ""
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0) {
      break
    }
    if(add){
      sequenceList[[title]] = line
      add = FALSE
    }
    if(substr(line, 1, 1)==">"){
      add = TRUE
      title = substr(line, 2, nchar(line))
      title = strsplit(title, ' ')[[1]][[1]]
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
scoreSequences2probabilityVector <- function(sequences){
  #1 concat sequences
  result = c()
  master_seq = c()
  for(l in sequences)
    master_seq = append(master_seq, l)
  #2 find max and min
  max = max(master_seq)
  min = min(master_seq)
  #3 count occurence of each
  for(x in min:max)
    result = append(result, length(which(master_seq == x)))
  names(result)<-min:max
  
  return (result/length(master_seq))
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
#' @param sequences sequences to be analysed (named list)
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
#' automatic_analysis(l,"iid",scores=c(-3,1),distribution=c(0.3,0.3,0.1,0.1,0.2))
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
  if((missing(model) && missing(transition_matrix) && missing(distribution)) || (missing(model) && !missing(transition_matrix) && !missing(distribution)))
    stop("unclear which model is used for sequence creation. Please provide either the model to be used, a transition matrix or a distribution vector.")
  #2 if model is missing, deduce it from perhaps present transition matrix
  if(missing(model)){
    if(missing(transition_matrix))
      model = "iid"
    else
      model = "markov"
  }
  #3 sequences provided? If not, ask user to pick file
  if(missing(sequences)){
    sequence_filepath = file.choose()
    sequences = loadCharSequencesFromFile(sequence_filepath)
  }
  #4 Sequences numeric? If yes, no scores needed, else pick file and use it
  if(!is.numeric(sequences[[1]]) && missing(scores)){
    score_filepath = file.choose()
    scores = loadScoreFromFile(score_filepath)
    sequences = CharSequences2ScoreSequences(sequences, scores)
  } 
  else if (!is.numeric(sequences[[1]]) && missing(scores)){
    #use scores provided in argument
    sequences = CharSequences2ScoreSequences(sequences, scores)
  }
  
  #4.1 declare result list
  results = list()
  #5 if extremes of scores are missing, extract them from sequences
  if(missing(score_extremes)){
    max = sequences[[1]][[1]]
    min = sequences[[1]][[1]]
    for(sequence in sequences){
      if(max(sequence)>max)
        max = max(sequence)
      if(min(sequence)<min)
        min = min(sequence)
    }
    score_extremes = c(min, max)
  }
  #6 Setup Progress Bar
  if (interactive())
    progressBar = txtProgressBar(min = 0, max = length(sequences), initial = 0, char = "=",
                                 width = NA, "progress in sequence analysis", "progressbar", style = 3, file = "")
  progressCounter = 0;
  #7 model specific behaviour
  if(model == "markov"){
    #7.1 check if transition matrix available: If not, build from sequences
    if(missing(transition_matrix))
      transition_matrix = sequences2transmatrix(sequences)[["Transition Matrix"]]
    #7.2 check if modelFunc provided
    if(missing(modelFunc)){
      #use exact Method if sequences are not too long
      for (name in names(sequences)){
        #calculate local score
        localscore = localScoreC(sequences[[name]],supressWarnings = TRUE)
        sequence_length = length(sequences[[name]])
        if(sequence_length<=method_limit){
          method = "Exact Method"
          p_value = exact_mc(localscore$localScore[1], transition_matrix, sequence_length, score_extremes[1]:score_extremes[2])
#          p_value = exact_mc(transition_matrix, localscore$localScore[1], sequence_length, sequence_min = score_extremes[1], sequence_max = score_extremes[2])
          results[[name]] = list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        } else {
          warning(paste(c("[Ignoring Sequence] sequence length ", name, " is of length ", sequence_length), collapse=""))
          results[[name]] = list("localScore" = localscore)
        }
        progressCounter = progressCounter + 1
        if (interactive())
          setTxtProgressBar(progressBar, progressCounter)
      }
    } else {
      #Function supplied, using simulation methods for Markov Chain p value determination
      for (name in names(sequences)){
        #calculate local score
        localscore = localScoreC(sequences[[name]],supressWarnings = TRUE)
        sequence_length = length(sequences[[name]])
        if(sequence_length<=method_limit){
          method = "Monte Carlo"
          p_value = monteCarlo(localscore$localScore[1], modelFunc, ...)
          results[[name]] = list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        } else {
          if(mean(sequence[[name]])>=0){
            warning(paste(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length), collapse=""))
            results[[name]] = list("localScore" = localscore)
          }
          else{
            method = "Asymptotic Method by Monte Carlo"
            res = karlinMonteCarlo(local_score = localscore$localScore[1], FUN = modelFunc, sequence_length = sequence_length, simulated_sequence_length = simulated_sequence_length , ...)
            results[[name]] = list("p-value" = res[["p-value"]], "method applied" = method, "localScore" = localscore, "K*" = res[["K*"]], "lambda" = res[["lambda"]])
          }
        }
        progressCounter = progressCounter + 1
        if (interactive())
          setTxtProgressBar(progressBar, progressCounter)
      }
    }
  } 
  else if (model == "iid"){
    #7.1 check if probability vector provided or in score (loading from file)
    if(missing(distribution)){
      if(!missing(scores) && length(scores)>2){
        #distribution provided with score: adjust extremes if file has greater extreme values
        distribution = scoreDictionnary2probabilityVector(scores)
        score_extremes = c(min(as.numeric(names(distribution))), max(as.numeric(names(distribution))))
      }
      else{
        #learn from sequences
        distribution = scoreSequences2probabilityVector(sequences)
      }
    }
    
    # Moyenne des scores sur l'ensemble de la liste
    MeanDistribution = sum(as.integer(names(distribution))*distribution)
    
    if(missing(modelFunc)){
      for (name in names(sequences)){
        #calculate local score
        localscore = localScoreC(sequences[[name]], supressWarnings = TRUE)
        sequence_length = length(sequences[[name]]) 
        if(sequence_length<=method_limit){
          method = "Exact Method Daudin et al"
          p_value = daudin(localscore$localScore[1], sequence_length, distribution, score_extremes[1], score_extremes[2])
          results[[name]] = list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        }
        else if(sequence_length<=method_limit*10){
          if(MeanDistribution<0){
            method = "Asymptotic Method Karlin et al"
            p_value = karlin(localscore$localScore[1], sequence_length, distribution, score_extremes[1], score_extremes[2])
            results[[name]] = list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
          } else {
            warning(paste(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length), collapse=""))
            results[[name]] = list("localScore" = localscore)
          }
        } else {
          if(MeanDistribution<0){
            method = "Improved Asymptotic Method Mercier et al"
            p_value = mcc(localscore$localScore[1], sequence_length, distribution, score_extremes[1], score_extremes[2])
            results[[name]] = list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
          } else {
            warning(paste(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length), collapse=""))
            results[[name]] = list("localScore" = localscore)
          }
        }
        progressCounter = progressCounter + 1
        if (interactive())
          setTxtProgressBar(progressBar, progressCounter)      
      }
    } else {
      #Function supplied, using simulation methods for iid p value determination
      for (name in names(sequences)){
        #calculate local score
        localscore = localScoreC(sequences[[name]],supressWarnings = TRUE)
        sequence_length = length(sequences[[name]])
        if(sequence_length<=method_limit){
          method = "Monte Carlo"
          p_value = monteCarlo(localscore$localScore[1], FUN = modelFunc, ...)
          results[[name]] = list("p-value" = p_value, "method applied" = method, "localScore" = localscore)
        } else {
          if(MeanDistribution>=0){
            warning(paste(c("[Ignoring Sequence] sequence mean >= 0 and sequence length = ", sequence_length), collapse=""))
            results[[name]] = list("localScore" = localscore)
          }else{
            method = "Asymptotic Method by Monte Carlo"
            res = karlinMonteCarlo(local_score = localscore$localScore[1], FUN = modelFunc, sequence_length = sequence_length, ...)
            results[[name]] = list("p-value" = res[["p-value"]], "method applied" = method, "localScore" = localscore, "K*" = res[["K*"]], "lambda" = res[["lambda"]])
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
#' @param list vector of scores
#' @param score_extremes vector of probability
#' @return vector containing all score values between extremes and the probability equal to 0 for missing score
#' @examples
#' Mylist=list("x1"=c(-2,0.1),"x2"=c(0,0.7),"x3"=c(1,0.2))
#' scoreDictionnary2probabilityVector(list=Mylist,score_extremes=c(-2,1))
#' @export
scoreDictionnary2probabilityVector <- function(list, score_extremes){
 vector = c()
 vec_names = c()
 for(item in list){
   vector = c(vector, item[[2]])
   vec_names = c(vec_names, item[[1]])
 }
 names(vector) = vec_names
 #check for missing values and set zeros if not present
 #get extremes, then create list and check difference
 missingScores = setdiff(min(as.numeric(names(vector))):max(as.numeric(names(vector))),as.numeric(names(vector)))
 missingValues = numeric(length = length(missingScores))
 names(missingValues) = missingScores
 vector = c(vector, missingValues)
 #order from lowest to highest score
 vector = vector[order(factor(as.numeric(names(vector))))]
 #si la fonction est appele avec des valeurs "limitantes" on coupe ce qui depasse ces valeurs
 if(!missing(score_extremes)){
   vector = subset(vector, as.numeric(names(vector))>=score_extremes[1] & as.numeric(names(vector))<=score_extremes[2])
 }
 return(vector)
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
  
  ExtendedIntegerScore <- MinIntegerScore:MaxIntegerScore
  
  ProbExtendedIntegerScore <- rep(0, times = length(ExtendedIntegerScore))
  
  ProbExtendedIntegerScore[ExtendedIntegerScore %in% IntegerScore] <- ProbRealScore
  names(ProbExtendedIntegerScore) <- ExtendedIntegerScore
  
  list(ExtendedIntegerScore = ExtendedIntegerScore,
       ProbExtendedIntegerScore = ProbExtendedIntegerScore)
}

