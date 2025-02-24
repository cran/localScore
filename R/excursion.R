#' Functions parameters global description 
#' @param theta  alphabet used (vector of character)
#' @param lambda transition probability matrix of size, theta x theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param a score strictly positive
#' @param prob0 Distribution of the first element of the Markov chain. Default to stationary distribution of lambda
#' @name MarkovParameters
NULL
#> NULL

#' @title compute all the descending trajectories of length c from (l-s)>=0 to l<0 with scores include in (Smin:-1), under condition that for every 
#' trajectory (t_1,...,t_c): ((l-s)+t_1+...+t_(c-1))>=0
#' some notations are the same as in the paper:"Exact distribution of excursion heights of the Lindley process in a Markovian model"
#' @param s -1*distance between starting point and arrival point of the trajectory (s<0)
#' @param c value in 1:(-s+1) represents the number of step of the trajectory (same notation as in the paper)
#' @param storage empty vector, used to store the result
#' @param l arrival value of the trajectory (l<0) (same notation as in the paper)
#' @param Smin lowest score possible to use in the trajectory
#' @param allowed allowed values for t_i
#' @return vector containing the trajectories
#' @examples
#' ######The following example should return an empty list
#' s=-4 #highest score=3
#' c=6
#' storage=rep(NA,c)
#' l=-1
#' Smin=-3
#' allowed = seq(-1, Smin, -1)
#' part(s,c,storage,l,Smin,allowed)
#' ######The following example should return this list: (-1,-1,-1)
#' s=-3 #highest score=2
#' c=3
#' storage=rep(NA,c)
#' l=-1
#' Smin=-2
#' allowed = seq(-1, Smin, -1)
#' part(s,c,storage,l,Smin,allowed)
#' ######The following example should return a list of 6 vectors
#' s=-5 #highest score =4
#' c=3
#' storage=rep(NA,c)
#' l=-1
#' Smin=-4
#' allowed = seq(-1, Smin, -1)
#' part(s,c,storage,l,Smin,allowed)
#' @noRd
part <- function(s, c, storage, l, Smin, allowed) {
  return_val = list()
  if ( c == 0 ) {
    if ( s == 0 ) {
      return_val[[length(return_val) + 1]] = storage
    }
  } else if ( c == 1 ) {
    if ( s %in% allowed ) {
      if ( s <= l) {
        storage[[length(storage)]] <- s
        return_val[[length(return_val) + 1]] = storage
      }
    }
  } else if ( (Smin * c) <= s & s <= (-c) ) {
    for (ele in allowed) {
      storage[(length(storage) - c) + 1] <- ele
      part_val = part((s-ele), (c-1), storage, l, Smin, allowed)
      for ( y in part_val) {
        return_val[[length(return_val) + 1]] = y
      }
    }
  }
  return(return_val)
}

#' @title compute all the descending trajectories from j to i with scores include in {Smin;-1}, 
#' under condition that for every trajectory (t_1,...,t_c): (j+t_1+...+t_(c-1))>=0
#' @param l arrival value (should be strictly negative)
#' @param j starting value (should be strictly positive)
#' @param Smin lowest possible score
#' @return list containing trajectories stored as vectors
#' @noRd
#' @examples
#' arrival_value=-4
#' starting_value=5
#' lowest_score=-10
#' traj(l=arrival_value,j=starting_value,Smin=lowest_score)
traj <- function(l, j, Smin) {
  trajectories = c()
  allowed = seq(-1, Smin, -1)
  for (c in 1:(j-l+1)) {
    trajectories <- c(trajectories, part((l-j), c, rep(NA, c), l, Smin, allowed))
  }
  return (trajectories)
}

#' 
#' #' @title Helper function used for the computation of all possible trajectories. Used in the computation of the matrix M
#' #' @param i value between -1 and the lowest score
#' #' @param j value between 1 and the highest score
#' #' @param m lowest score
#' #' @return all possibles trajectories
#' #' @examples
#' #' possible_tuples = lapply(min_val:(-1), function(ell) {
#' #' possible_hill_ending = lapply(1:max_val, function(j){
#' #'   traj(ell, j, min_val)
#' #' })
#' #' names(possible_hill_ending) = c(1:max_val)
#' #' return(possible_hill_ending)
#' #' })
#' traj <- function(i, j, m) {
#'   trajectories = lapply(traj_aux(i,j,m), function(i) c(i))
#'   return(trajectories)
#' }
#' 


#' @title Helper function to turn string back into list
#' @param str containing a letter of theta and a score separated by ":"
#' @return list containing a key (the letter) and a score
#' @examples
#' alpha <- "K:1"
#' str_to_alphabeta(alpha)
#' @noRd
str_to_alphabeta <- function(str) {
  tmp <- unlist(strsplit(str, split = ":", fixed = TRUE))
  return(list(key = tmp[1], score = as.integer(tmp[2])))
}

#' @title cell_value, compute the value of a cell for the matrix P
#' @param alpha letters used as indices of the row
#' @param beta  letters used as indices of the column
#' @param theta  alphabet used (vector of character)
#' @param a threshold score
#' @param lambda transition probability matrix of size, theta x theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @return reel number
#' @examples
#' cell_value(alpha=list(key="L",score=4),
#'            beta=list(key="K",score=3),
#'            theta=c('K','L'),
#'            a=5,
#'            lambda=matrix(c(0.1, 0.9, 0.2, 0.8), ncol = 2, byrow = TRUE),
#'            score_function=c(4,-1))
#' @noRd
cell_value <- function(alpha, beta, theta, a, lambda, score_function) {
  if ( alpha$score == beta$score & ( alpha$score == -1 | alpha$score == as.double(a) ) ) { #XXX Ajouter condition alpha$key == beta$key ?
    return(1)
  } 
  
  index_alpha <- match(alpha$key, theta)
  index_beta  <- match(beta$key , theta)
  
  lambda_a_b <- lambda[index_alpha, index_beta]
  score <- score_function[index_alpha]
  
  if ( !(alpha$score >= 0 & alpha$score <= a-1) ) {
    return(0) 
  }
  
  if ((beta$score == -1) & (alpha$score + score < 0)) {
    return(lambda_a_b)
  } else if ((beta$score == a) & (a <= alpha$score + score) ) {
    return(lambda_a_b)
  } else if ((beta$score == alpha$score + score) & (alpha$score + score >= 0) & (alpha$score + score <= a - 1)) {
    return(lambda_a_b)
  } else {
    return(0)
  }
}


#' @title create_matrix_P, compute the transition probability matrix P of the stopped process Y, same notation as in the 
#' paper:"Exact distribution of excursion heights of the Lindley process in a Markovian model"
#' @param a threshold score
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix, same size and order as theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param combination vector of strings, each string containing a letter of theta and a score separated by ":"
#' @return a transition probability matrix, P
#' @examples
#' create_matrix_P(a=5,
#'                 theta=c("K","L"),
#'                 lambda=matrix(c(0.1, 0.9, 0.2, 0.8), ncol = 2, byrow = TRUE),
#'                 score_function=c(4,-1),
#'                 combination=c("K:-1","L:-1","K:0","L:0","K:1","L:1","K:2","L:2","K:3","L:3","K:4","L:4","K:5","L:5"))
#' @noRd
create_matrix_P <- function(a, theta, lambda, score_function, combination) {
  P <- matrix(nrow = length(combination), ncol = length(combination))
  dimnames(P) <- list(combination, combination)
  i <- 1
  j <- 1
  for (alpha in combination) {
    for (beta in combination) {
      decomp_alpha <- str_to_alphabeta(alpha)
      decomp_beta  <- str_to_alphabeta(beta)
      P[i, j] = cell_value(decomp_alpha, decomp_beta, theta, a, lambda, score_function)
      j <- j + 1
    }
    j <- 1
    i <- i + 1
  }
  P <- P[-(length(combination) - length(theta) + 1):-(length(combination)),]
  P <- P[-1:-length(theta),]
  return(P)
}


#' @title Probability P_alpha(Q(1)>=a) that the height of the first excursion >=a given the sequence begins with the letter alpha given a markovian model on the sequence.
#' @param a threshold score
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix, square, same size and order as theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param alpha first letter of the sequence
#' @return s: reel number between 0 and 1 representing the theoretical probability of reaching a score of a on the first excursion
#' #' @details : The first letter 'alpha' is not included in the score sequence. The conditioning is BEFORE the sequence.
#' Also beware that a sequence beginning with a negative score gives a "flat" excursion, with score 0 are considered.
#' It is also used in \code{\link{proba_theoretical_first_excursion()}}.
#' @examples
#' proba_theoretical_first_excursion_alpha(a=5,
#'                  theta=c("K","L","M"),
#'                  lambda=matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE),
#'                  score_function=c(-2,-1,2),
#'                  alpha="L")
#' @noRd
proba_theoretical_first_excursion_markov_alpha <- function(a, theta, lambda, score_function, alpha) {
  list.theta <- as.list(theta)
  list.a <- as.list(c(-1:a))
  combination <- as.vector(outer(list.theta, list.a, paste, sep = ":"))
  
  P <- create_matrix_P(a, theta, lambda, score_function, combination)
  
  Q_C <- P[,-(length(combination) - length(theta) + 1):-(length(combination))]
  Q_C <- Q_C[,-1:-length(theta)]
  Q_B <- P[,1:length(theta)]                                            
  Q_A <- P[,(length(combination) - length(theta) + 1):length(combination)]
  
  z <- (solve(diag(ncol(Q_C)) - Q_C)) %*% (Q_A %*% rep(1, ncol(Q_A)))
  s <- 0 
  
  index_alpha <- match(alpha, theta)
  
  for (delta in theta) {
    index_delta <- match(delta, theta)
    score <- score_function[index_delta]
    lambda_alpha_delta <- lambda[index_alpha, index_delta]
    
    if (0 <= score & score <= (a-1)) {
      for (mu in theta) {
        index_mu <- match(mu, theta)
        lambda_delta_mu <- lambda[index_delta, index_mu]  
        
        key <- paste(mu, score, sep = ":")
        index_z <- match(key, rownames(z))     
        s <- s + lambda_alpha_delta * lambda_delta_mu * z[index_z]
      }
    }
    if (score >= a) {
      s <- s + lambda_alpha_delta
    }
  }
  return( list( s = s, z = z, p=P) )
}

#' @title generalize_to_all_min, takes a list of trajectories to -1 and return all the trajectories with arrival value from -1 to Smin
#' @description the idea of this function is that the trajectories from s to -2 are the same as those from s to -1 if we remove 1 to the last element of each trajectory,
#' with the condition that the last element of each trajectory have to be >=Smin
#' Useful because it's faster than computing all trajectories with part() function (linear complexity vs exponential)
#' @param list_of_trajectories the list of trajectories from any starting point to -1 with scores included in (Smin:-1), same format as the return of traj() function
#' @param Smin the minimal score, trajectories will be compute with scores included in (Smin:-1)
#' @return a list of lists of vector, containing all the trajectories from the starting point of list_of_trajectories to all the arrival value included in (Smin:-1)
#' generalize_to_all_min[[i]][[j]] is the j-th trajectory from starting point to (-i)
#' @examples
#' ###The following example should return a list of 9 lists of trajectories
#' generalize_to_all_min(traj(-1,4,-10),-10)
#' ###The following example should return a list (of size one) of list (of size one) 
#' with this vector: (-1,-1,-1,-1,-1) generalize_to_all_min(traj(-1,4,-1),-1)
#' @noRd
generalize_to_all_min <- function(list_of_trajectories, Smin)  {
  generalized_list = list()
  i = 0
  for (j in 1:(-Smin)){
    k=1
    generalized_list[[j]] = list()
    generalized_list[[j]][k][[1]] = rep(NA, length(list_of_trajectories[i+1][[1]]))
    for (i in 1:length(list_of_trajectories)){
      if(list_of_trajectories[i][[1]][length(list_of_trajectories[i][[1]])]-(j-1)>=(Smin)){
        generalized_list[[j]][k][[1]] = c(list_of_trajectories[i][[1]][1:length(list_of_trajectories[i][[1]])-1], list_of_trajectories[i][[1]][length(list_of_trajectories[i][[1]])]-(j-1))
        k = k + 1
      }
    }
  }
  return (generalized_list)
}

#' @title Probability \eqn{P(Q(1)\geq a)} that the height of the first excursion is greater or equal to \code{a} given a Markov model on the letters sequence
#' @description Mathematical definition of an excursion of the Lindley process is based on the record times of the partial
#' sum sequence associated to the score sequence (see Karlin and Altschul 1990, Karlin and Dembo 1992) and
#' define the successive times where the partial sums are strictly decreasing. There must be distinguished
#' from the visual excursions of the Lindley sequence. Detailed definitions are given in the documentation.
#' @param a score strictly positive
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix of size, theta x theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param prob0 Distribution of the first letter of the sequence (must be in the same order as theta). If NULL, stationary distribution of lambda is calculated and applied.
#' @return theoretical probability of reaching a score of a on the first excursion given the first letter of the sequence follow the distribution given by parameter 'prob0"
#' @details Beware that a sequence beginning with a negative score gives a "flat" excursion, with score 0, are considered.
#' @seealso \code{\link{proba_theoretical_ith_excursion_markov()}}
#' @examples
#' theta = c("a","b","c","d")
#' lambda = matrix(c(0.1,0.2,0.4,0.3,
#'                   0.3,0.1,0.5,0.1,
#'                   0.2,0.6,0.1,0.1,
#'                   0.4,0.1,0.1,0.4),
#'                 ncol=4, byrow = TRUE, dimnames = list(theta,theta))
#' all(apply(lambda,1,sum)==1) #TRUE for markov transition matrix
#' stationary_dist = stationary_distribution(lambda)
#' score_function = c(a=-3,b=-1,c=1,d=2)
#' sum(score_function*stationary_dist) #Score expectation (should be <0)
#' prob0 = c(.25,.25,.25,.25) #
#' proba_theoretical_first_excursion_markov(3, theta, lambda, score_function, prob0)
#' proba_theoretical_ith_excursion_markov(3, theta, lambda, score_function, i=1, prob0) #for comparison
#' @noRd
proba_theoretical_first_excursion_markov <- function(a, theta, lambda, score_function, prob0=NULL) {
  # RM théoriquement plus utilisée
  if (is.null(prob0)){
    prob0 = stationary_distribution(lambda)
  }

  if (a <= 0) {
    return(1.0)
  }
  
  prob = 0.0
  for (i in 1:length(theta)) {
    if (prob0[i] != 0.0) {
      p <- proba_theoretical_first_excursion_markov_alpha(a, theta, lambda, score_function, theta[i])
      prob  <- prob + prob0[i]*p$s
    }
  }
  return(prob)
}


#' @title Transition matrix of the start of excursions, used for function computing probability of the i-th excursion, same notation as in the 
#' paper:"Exact distribution of excursion heights of the Lindley process in a Markovian model"
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix, square, same size and order as theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param N_iter maximal number of iteration (recommended : 1000)
#' @param epsilon threshold of the difference between matrices of two following iteration, this can reduce the number of iterations
#' @return transition probability matrix, square, same size and order as theta
#' @examples
#' matrix_M(c("K","L","M"), matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), c(-2,-1,2), 1000)
#' @noRd
matrix_M <- function(theta, lambda, score_function, N_iter, epsilon = 1e-16) {
  min_val <- min(score_function)
  max_val <- max(score_function)
  positive_scores <- score_function[score_function > 0]
  possible_tuples <- lapply(1:max_val, function(max) if (max %in% positive_scores) {generalize_to_all_min(traj(-1, max, min_val), min_val)})
  lambda_matrices <- lapply(min_val:max_val, function(j) {
    intermediary_lambda_matrix <- matrix(0, length(theta), length(theta))
    if ( j %in% score_function) {
      indexes = which(score_function %in% j)
      for ( ind in indexes ) {
        intermediary_lambda_matrix[, ind] <- lambda[, ind]
      }
    }
    return( intermediary_lambda_matrix )
  })
  names(lambda_matrices) <- c(min_val:max_val)
  M_iter_aux <- lambda_matrices[as.character(min_val:-1)]
  names(M_iter_aux) <- as.character(min_val:(-1))
  M_iter <- lapply(min_val:(-1), function(i) {
    return(matrix(0, length(theta), length(theta)))
  })
  names(M_iter) <- c(min_val:(-1))
  i = 1
  norme = rep(1, abs(min_val))
  eps = rep(epsilon, abs(min_val))
  while((i < N_iter) &! ((i > 10) & (min(norme < eps)))) {
    for ( ell in min_val:(-1) ) {
      ell_char = as.character(ell)
      zero_char = as.character(0)
      M_iter[[ell_char]] = lambda_matrices[[ell_char]] + lambda_matrices[[zero_char]] %*% M_iter_aux[[ell_char]]
      for ( j in unique(positive_scores) ) {
        j_char = as.character(j)
        for ( tuple in possible_tuples[[j]][[(-ell)]] ) {
          matr_prod = lambda_matrices[[j_char]]
          for ( step in tuple ) {
            matr_prod <- matr_prod %*% M_iter_aux[[as.character(step)]]
          }
          M_iter[[ell_char]] <- M_iter[[ell_char]] + matr_prod
        }
      }
    }
    for ( ell in min_val:(-1) ) {
      ell_char = as.character(ell)
      norme[[-ell]] = max(abs(M_iter_aux[[ell_char]] - M_iter[[ell_char]]))
      M_iter_aux[[ell_char]] = M_iter[[ell_char]]
    }
    i = i + 1
  }
  if (i == N_iter & !min(norme < eps)){
    warning("The computation of the transition probability matrix M has stopped without convergence. It reached the maximal number of iterations and the norm is not < epsilon.")
  }
  M <- matrix(0, length(theta), length(theta))
  for (ell_char in names(M_iter)) {
    M <- M + M_iter[[ell_char]]
  }
  return ( M )
}


#' @title theoretical probability of reaching the threshold score a on the i-th excursion (sequential order) of a markov's score sequence
#' @description
#' This function implements the results of the paper: "Exact distribution of excursion heights of the Lindley process in a Markovian model", this is an
#' exact method. 
#' Scores of score function have to be integers, the expectancy of the score have to be negative, the score function have to contains at least
#' one positive integer, i have to be >0
#' @details
#' Be careful: the computational time is exponential in function of the maximum score of the score function.
#' The computational time also increase with the rise of the threshold score a and with the lowering of the minimum of the score function and have a cubic complexity in function of the size of theta.
#' Lowering epsilon can also decrease the execution time of the function but it can have an impact on the accuracy of the probabilities.
#' @param a threshold score
#' @param theta vector containing the alphabet used
#' @param lambda transition probability matrix of the markov chain, square, same size and order as theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @param i number of wanted excursion (sequential order)
#' @param prob0 probability distribution of the first letter of the sequence (stationary distribution of lambda if the parameter is NULL),
#' used only for computation of proba_q_i_geq_a
#' @param epsilon threshold for the computation of matrix M: to calculate the probabilities, this function need to compute the matrix M (transition probability matrix of the beginning of excursions, as in the paper)
#' with a recurrence, this recurrence stop at a certain number of iterations or when the maximal absolute value of the difference between an iteration and the next one is < epsilon.
#' @return a list containing:\tabular{ll}{
#'     \code{proba_q_i_geq_a} \tab reel number between 0 and 1 representing the theoretical probability that the i-th excursion is greater or equal a \cr
#'     \tab \cr
#'     \code{P_alpha} \tab vector containing the probabilities (of reaching the threshold score on the i-th excursion) depending on the first letter of the sequence (alpha) \cr
#' }
#' @examples
#' ## In this example: the probability to reach a score of a in the third excursion is: 0.04237269
#' ## The conditional probabilities to reach a score of a in the third excursion if the first
#' ## letter of the sequence is K, L or M are respectively: 0.04239004, 0.04247805 and 0.04222251
#' proba_theoretical_ith_excursion_markov(a = 5,
#'                theta = c("K","L","M"),
#'                lambda = matrix(c(0.5, 0.3, 0.2,
#'                                  0.4, 0.2, 0.4,
#'                                  0.4, 0.4, 0.2), 
#'                                ncol = 3, byrow = TRUE),
#'                score_function = c(-2,-1,2),
#'                i = 3,
#'                prob0 = c(0.4444444, 0.2962963, 0.2592593))
#' ### This example implements the numerical application of the paper,
#' ###    the global probability is 0.2095639
#' proba_theoretical_ith_excursion_markov(a = 6,
#'             theta = c("a","b","c","d","e"),
#'             lambda = matrix(c(0.1, 0.7, 0.05, 0.05, 0.1,
#'                               0.3, 0.3, 0.2, 0.15, 0.05,
#'                               0.1, 0.4, 0.15, 0.2, 0.15,
#'                               0.5, 0.05, 0.25, 0.1, 0.1,
#'                               0.25, 0.05, 0.5, 0.15, 0.05),
#'                               ncol = 5,nrow=5, byrow = TRUE),
#'             score_function = c(-3,-2,-1,6,7),
#'             i = 3)
#' @export
proba_theoretical_ith_excursion_markov <- function(a, theta, lambda, score_function, i, prob0 = NULL, epsilon = 1e-16) {
  theta = as.character(theta)
  if(is.null(prob0)){
    prob0 = stationary_distribution(lambda)
  }
  for (x in score_function){
    if(x%%1){
      stop("[Invalid Input] scores of score_function should be integers.")
    }
  }
  if(sum(score_function * stationary_distribution(lambda)) >= 0){
    stop("[Invalid Input] score expectation should be negative.")
  }
  if(all(score_function <= 0)){
    stop("[Invalid Input] score_function should have at least one positive integer.")
  }
  if(length(unique(theta)) != length(theta)){
    stop("[Invalid Input] theta should contains only unique values")
  }
  if (i < 1){
    stop("[Invalid Input] i should be a strictly positive integer.")
  }
  if(!length(theta) == dim(lambda)[1]){
    stop("[Invalid Input] lambda should be the same size as theta.")
  }
  if (a == 0){
    P_alpha = rep(1,length(theta))
    names(P_alpha) = theta
    return(list(proba_q_i_geq_a = 1, P_alpha = P_alpha))
  }
  P_alpha = c()
  for (alpha in theta) {
    P_alpha[alpha] = proba_theoretical_first_excursion_markov_alpha(a, theta, lambda, score_function, alpha)$s
  }
  proba_q_i_geq_a = 0
  for (alpha in theta) {
    index_alpha = which(theta == alpha)
    proba_q_i_geq_a = proba_q_i_geq_a + (prob0[index_alpha] * P_alpha[alpha])
  }
  proba_q_i_geq_a = as.double(proba_q_i_geq_a)
  if ( i == 1 ) { 
    return ( list(proba_q_i_geq_a = proba_q_i_geq_a, P_alpha = P_alpha))
  }
  M = matrix_M(theta, lambda, score_function, N_iter = 1000, epsilon = epsilon)
  if (i > 1) {
    previous_exc = P_alpha
    M_i = M
    while (i > 2){
      M_i = M %*% M_i
      i = i - 1
    }
    P_alpha = M_i %*% P_alpha
    P_alpha = as.vector(P_alpha)
    names(P_alpha) = theta
  }
  proba_q_i_geq_a = 0
  for( alpha in theta ) {
    index_alpha = which(theta == alpha)
    proba_q_i_geq_a = proba_q_i_geq_a + (prob0[index_alpha] * P_alpha[[alpha]])
  }
  return (list(proba_q_i_geq_a = proba_q_i_geq_a, P_alpha = P_alpha))
}

#' @title Probability \eqn{P(Q(1)\geq a)} that the height of the first excursion is greater or equal to \code{a} given a i.i.d. model on the letters sequence
#' @description Mathematical definition of an excursion of the Lindley process is based on the record times of the partial
#' sum sequence associated to the score sequence (see Karlin and Altschul 1990, Karlin and Dembo 1992) and
#' define the successive times where the partial sums are strictly decreasing. There must be distinguished
#' from the visual excursions of the Lindley sequence. The number \code{i} is the number of excursion in sequential order. Detailed definitions are given in the vignette.
#' @param a score strictly positive
#' @param theta vector containing the alphabet used
#' @param theta_distribution distribution vector of theta
#' @param score_function vector containing the scores of each letters of the alphabet (must be in the same order as theta)
#' @return theoretical probability of reaching a score of a on the first excursion supposing an  i.i.d model on the letters sequence
#' @details Beware that a sequence beginning with a negative score gives a "flat" excursion, with score 0 are considered.
#' @examples
#' proba_theoretical_first_excursion_iid(3, c("a","b","c","d"), 
#'                                       c(a=0.1,b=0.2,c=0.4,d=0.3), c(a=-3,b=-1,c=1,d=2))
#' @export
proba_theoretical_first_excursion_iid <- function(a, theta, theta_distribution, score_function){
  #Trick : use markov results with transition matrix with equal lines
  if (a <= 0) {
    return(1.0)
  }
  #Build transition matrix
  l <- length(theta)
  lambda = matrix(ncol = l,nrow = l, byrow = TRUE, dimnames = list(theta, theta))
  for (i in 1:l) {
    lambda[i,] <- theta_distribution
  }
  p <- proba_theoretical_first_excursion_markov_alpha(a, theta, lambda, score_function, theta[1])$s
  return(p)
}

#'@title Probability \eqn{P(Q(i)\geq a)} that the height of the ith excursion (sequential order) is greater or equal to \code{a} given a i.i.d. model on the letters sequence
#'@inherit proba_theoretical_first_excursion_iid  description params 
#'@param i Number of excursion in sequential order
#'@return theoretical probability of reaching a score of a on the ith excursion supposing an  i.i.d model on the letters sequence
#'@details
#'  In the i.i.d., the distribution of the ith excursion is the same as the first excursion. This function is just for convenience, and the result is the same as \code{proba_theoretical_first_excursion_iid}. Beware that a sequence beginning with a negative score gives a "flat" excursion, with score 0 are considered.
#'@examples
#'p1 <- proba_theoretical_ith_excursion_iid(3, c("a","b","c","d"), 
#'                                       c(a=0.1,b=0.2,c=0.4,d=0.3), c(a=-3,b=-1,c=1,d=2), i = 10)
#'p2 <- proba_theoretical_first_excursion_iid(3, c("a","b","c","d"), 
#'                                       c(a=0.1,b=0.2,c=0.4,d=0.3), c(a=-3,b=-1,c=1,d=2))
#'p1 == p2  #TRUE
#'@seealso \code{\link{proba_theoretical_first_excursion_iid}}
#'@export
proba_theoretical_ith_excursion_iid <- function(a, theta, theta_distribution, score_function, i = 1){
  return(proba_theoretical_first_excursion_iid(a, theta, theta_distribution, score_function))
}

