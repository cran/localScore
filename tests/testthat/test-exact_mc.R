context("Tests of exact_mc function (Markov model case)")

##### ERROR HANDLING ############
# Error link to #10215 ticket on forge DGA
test_that("Wrong matrix size error", {
  f <- c(-2,-1,0,1,2)
  min.f <- min(f)
  max.f <- max(f)
  nb_Alpha <- length(f)
  seq <- c(-2, -1, 1 , 1 , 1 , 1 )
  n <- length(seq)
  P <- matrix(ncol = 4, nrow = 4)
  P[1,] <- c(0.3, 0.3, 0.2, 0.2); P[2,] <- c(0.3, 0.3, 0.2, 0.2) ;
  P[3,] <- c(0.2, 0.3, 0.3, 0.2); P[4,] <- c(0.3, 0.2, 0.2, 0.3)
  LS <- 4
  expect_error(exact_mc(LS, P, sequence_length = n, score_values = f),
              "[ERROR exact_mc : Invalid Input] m should be a square matrix of size the length of score_values", fixed = TRUE)
})

############ 
test_that("Non ergodic Markov Chain", {
  #Liée à la demande #10408 sur la forge DGA.
  # Attention le code ci-desous est pour la tracabilité et n'est pas un vrai test (la valeur attendue ici n'a pas été vérifiée)
  min.f <- -2
  max.f <- 2
  f <- min.f:max.f
  nb_Alpha <- max.f - min.f + 1
  seq <- c(-2, -1, 1 , 1 , 1 , 1 )
  n <- length(seq)
  P <- matrix(ncol = nb_Alpha, nrow = nb_Alpha)
  P[1,] <- c(0.3, 0.3, 0.0, 0.2, 0.2)
  P[2,] <- c(0.3, 0.3, 0.0, 0.2, 0.2)
  P[3,] <- c(0.0, 0.0, 1.0, 0.0, 0.0)
  P[4,] <- c(0.2, 0.3, 0.0, 0.3, 0.2)
  P[5,] <- c(0.3, 0.2, 0.0, 0.2, 0.3)
  LS <- 4
  expect_error(exact_mc(LS, P, sequence_length = n, score_values = f), "Markov matrix is not irreductible (many eigenvalues == 1).", fixed = TRUE )
  ## Attention : valeur attendue ci-dessous non vérifiée théoriquement.
  expect_equal(exact_mc(LS, P, sequence_length = n, score_values = f, prob0 = c(1.0,0.0,0.0,0.0,0.0)), 0.30376)
})

test_that("Bad probability vector input", {
  mTransition <- t(matrix(c(0.2, 0.3, 0.5, 0.3, 0.4, 0.3, 0.2, 0.4, 0.4), nrow = 3))
  scoreValues <- -1:1
  expect_error( exact_mc(local_score = 12, m = mTransition, sequence_length = 100, prob0 = -1:1),"[ERROR exact_mc : Invalid Input] prob0 vector should sum to 1",fixed = TRUE)
  expect_error( exact_mc(local_score = 12, m = mTransition, sequence_length = 100, prob0 = c(0,2,-1)),"[ERROR exact_mc : Invalid Input] prob0 vector should contains values between 0 and 1",fixed = TRUE)
})

test_that("prob0 specifed or not", {
  P0 <- matrix(ncol = 4, nrow = 4)
  P0[1,] <- c(0.3, 0.3, 0.2, 0.2)
  P0[2,] <- c(0.3, 0.3, 0.2, 0.2)
  P0[3,] <- c(0.2, 0.3, 0.3, 0.2)
  P0[4,] <- c(0.3, 0.2, 0.2, 0.3)
  seq <- c(-2, -1, 1 , 1 , 1 , 1 )
  n <- length(seq)
  LS <- 4
  mu <- stationary_distribution(P0)
  expect_equal(exact_mc(LS, P0, n, -2:1, mu),0.01976)
  expect_equal(exact_mc(LS, P0, n, -2:1), 0.01976)
})


test_that("Named transition matrix or not",{
  P0 <- matrix(ncol = 4, nrow = 4)
  P0[1,] <- c(0.3, 0.3, 0.2, 0.2)
  P0[2,] <- c(0.3, 0.3, 0.2, 0.2)
  P0[3,] <- c(0.2, 0.3, 0.3, 0.2)
  P0[4,] <- c(0.3, 0.2, 0.2, 0.3)
  seq <- c(-2, -1, 1 , 1 , 1 , 1 )
  n <- length(seq)
  LS <- 4
  mu <- stationary_distribution(P0)
  expect_error(exact_mc(LS, P0, n, prob0 = mu), "[ERROR exact_mc : Invalid Input] Either m is a matrix  with score values as rownames or specify score values via score_values parameter", fixed = TRUE)
  colnames(P0) <- -2:1
  expect_error(exact_mc(LS, P0, n, prob0 = mu), "[ERROR exact_mc : Invalid Input] Either m is a matrix  with score values as rownames or specify score values via score_values parameter", fixed = TRUE)
  rownames(P0) <- -2:1
  expect_equal(exact_mc(LS, P0, n, prob0 = mu), 0.01976)  
})

test_that("Minimal call",{
  P0 <- matrix(ncol = 4, nrow = 4)
  P0[1,] <- c(0.3, 0.3, 0.2, 0.2)
  P0[2,] <- c(0.3, 0.3, 0.2, 0.2)
  P0[3,] <- c(0.2, 0.3, 0.3, 0.2)
  P0[4,] <- c(0.3, 0.2, 0.2, 0.3)
  seq <- c(-2, -1, 1 , 1 , 1 , 1 )
  n <- length(seq)
  LS <- 4
  rownames(P0) <- -2:1
  expect_equal(exact_mc(LS, P0, n), 0.01976)  
})

test_that("Comparison exact_mc with daudin() in cas i.i.d", {
  probs <- c(0.3, 0.1, 0.2, 0.2, 0.1, 0.1)
  P0 <- matrix(ncol = 6, nrow = 6)
  P0[1,] <- probs
  P0[2,] <- probs
  P0[3,] <- probs
  P0[4,] <- probs
  P0[5,] <- probs
  P0[6,] <- probs
  scoremin <- -3
  scoremax <- 2
  rownames(P0) <- scoremin:scoremax
  n <- 100
  LS <- 5
  pEMC <- exact_mc(LS, P0, n)
  pDaudin <- daudin(LS, n , probs, scoremin, scoremax)
  expect_equal(pEMC, pDaudin)
})


POUBELLE <- function() {
  exact_mc3(LS, P0, n)
  exact_mc(P, LS, sequence_length = n, sequence_min = -2, sequence_max = 2)
  exact_mc(P0, LS, sequence_length = n, sequence_min = -2, sequence_max = 1)
  #double exact_mc2(int localScore, NumericMatrix m, NumericVector score_values, long sequence_length, NumericVector prob0){
  exact_mc2(LS, P0, -2:1, n, rep(0.25,4))
  exact_mc3(LS, P0, n,-2:1, rep(0.25,4))
  
  mu <- stationary_distribution(P0)
  exact_mc3(LS, P0, n, -2:1, mu)
  exact_mc3(LS, P0, n, -2:1)
  
  P0.1 <- P0 ; rownames(P0.1) <- -2:1
  exact_mc3(LS, P0.1, n)
  exact_mc3(LS, P0, n) #should give error
  tmp <- sapply(0:10,exact_mc3, P0.1, n)
  sum(-diff(tmp)) # should be 1 if 10 is enough
}