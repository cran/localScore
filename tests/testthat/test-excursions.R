context("Exact distribution of first and ith excursions - theoretical methods - Markov and iid model")
test_that("Results of intermediate function: traj, test when l!=Smin", {
  expected_results = list(
    c(-3, -6),
    c(-4, -5),
    c(-1, -2, -6),
    c(-1, -3, -5),
    c(-2, -1, -6),
    c(-2, -2, -5),
    c(-3, -1, -5),
    c(-1, -1, -1, -6),
    c(-1, -1, -2, -5),
    c(-1, -2, -1, -5),
    c(-2, -1, -1, -5),
    c(-1, -1, -1, -1, -5)
  )
  result <- traj( l = -5, j = 4, Smin = -6)
  expect_equal(result, expected_results)
})

test_that("Results of intermediate function: traj, test when l==Smin", {
  expected_results = list(
    c(-3, -3, -3),
    c(-1, -2, -3, -3),
    c(-1, -3, -2, -3),
    c(-2, -1, -3, -3),
    c(-2, -2, -2, -3),
    c(-2, -3, -1, -3),
    c(-3, -1, -2, -3),
    c(-3, -2, -1, -3),
    c(-1, -1, -1, -3, -3),
    c(-1, -1, -2, -2, -3),
    c(-1, -1, -3, -1, -3),
    c(-1, -2, -1, -2, -3),
    c(-1, -2, -2, -1, -3),
    c(-1, -3, -1, -1, -3),
    c(-2, -1, -1, -2, -3),
    c(-2, -1, -2, -1, -3),
    c(-2, -2, -1, -1, -3),
    c(-3, -1, -1, -1, -3),
    c(-1, -1, -1, -1, -2, -3),
    c(-1, -1, -1, -2, -1, -3),
    c(-1, -1, -2, -1, -1, -3),
    c(-1, -2, -1, -1, -1, -3),
    c(-2, -1, -1, -1, -1, -3),
    c(-1, -1, -1, -1, -1, -1, -3)
  )
  result <- traj(l = -3, j = 6, Smin = -3)
  expect_equal(result, expected_results)
})

test_that("Results of intermediate function: traj, test when there is no trajectory", {
  expected_results = list()
  result <- traj(l = -3, j = 6, Smin = -1)
  expect_equal(result, expected_results)
})

test_that("Results proba theoretical third excursion normal case (Numerical application of the paper:\"Exact distribution of excursion heights of the Lindley
process in a Markovian model\")", {
  #Validated by simulation
  expected_results = list()
  expected_results$proba_q_i_geq_a = c(0.2095639, 0.2095639, 0.2095639, 0.2095639, 0.2095639, 0.2095639, 0.1524696, 0.1196005, 0.1154948, 0.1112926)
  expected_results$P_alpha = list(c(a = 0.2113536, b = 0.2088997, c = 0.2105256, d = 0.2065410, e = 0.2098564),
                                  c(a = 0.1560269, b = 0.1513963, c = 0.1540934, d = 0.1475148, e = 0.1510870),
                                  c(a = 0.1212326, b = 0.1190616, c = 0.1203996, d = 0.1171290, e = 0.1193353),
                                  c(a = 0.1170527, b = 0.1149792, c = 0.1162590, d = 0.1131304, e = 0.1152516),
                                  c(a = 0.1128484, b = 0.1107814, c = 0.1120513, d = 0.1089475, e = 0.1110196)
  )
  
  #a=1 to a=5 should return the same results as a=6, we don't verify them to reduce computational time
  for (a in 6:10){
    result <- proba_theoretical_ith_excursion_markov(a = a,
                                                     theta = c("a","b","c","d","e"),
                                                     lambda = matrix(c(0.1,0.7,0.05,0.05,0.1, 0.3,0.3,0.2,0.15,0.05, 0.1,0.4,0.15,0.2,0.15, 0.5,0.05,0.25,0.1,0.1, 0.25,0.05,0.5,0.15,0.05), ncol = 5,nrow=5, byrow = TRUE),
                                                     score_function = c(-3,-2,-1,6,7),
                                                     i = 3)
    expect_equal(result$proba_q_i_geq_a, expected_results$proba_q_i_geq_a[a], tolerance = 1e-6)
    expect_equal(result$P_alpha, expected_results$P_alpha[a-5][[1]], tolerance = 1e-6)
  }
})

test_that("Results proba theoretical ith excursion, normal case (with an integer in theta)", {
  #Validated by simulation
  expected_result = list()
  expected_result$proba_q_i_geq_a = 0.3083594
  expected_result$P_alpha = c(K = 0.2158516, "1" = 0.3453626)
  result <- proba_theoretical_ith_excursion_markov(a = 5,
                                                   theta = c("K", 1),
                                                   lambda = matrix(c(0.5, 0.5, 0.2, 0.8), ncol = 2, byrow = TRUE),
                                                   score_function = c(-4,1),
                                                   i = 1)
  expect_equal(result$proba_q_i_geq_a, expected_result$proba_q_i_geq_a, tolerance = 1e-6)
  expect_equal(result$P_alpha, expected_result$P_alpha, tolerance = 1e-6)
})

test_that("Results proba theoretical first excursion normal case (2)", {
  result = proba_theoretical_first_excursion_markov(a = 5, 
                                                    theta = c("K","L","M"), 
                                                    lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), 
                                                    score_function = c(-2,-1,2),
                                                    prob0 = c(0.4444444, 0.2962963, 0.2592593) )
  expect_equal(result, 0.04134028)
  #  expect_equal(as.vector(result$z), expected_results$z)
})

test_that("Results proba theoretical second excursion, specific case (transition matrix have only low or high values) ", {
  #Validated by simulation
  expected_result = list()
  expected_result$proba_q_i_geq_a = 0.01818579
  expected_result$P_alpha =  c(K = 0.01998373, L = 0.01800000, M = 0.00396745)
  result <- proba_theoretical_ith_excursion_markov(a = 4,
                                                   theta = c("K", "L", "M"),
                                                   lambda = matrix(c(0.9, 0, 0.1, 0.9, 0.1, 0, 0, 0.8, 0.2), ncol = 3, byrow = TRUE),
                                                   score_function = c(-3,-1,2),
                                                   i = 2)
  expect_equal(result$proba_q_i_geq_a, expected_result$proba_q_i_geq_a, tolerance = 1e-6)
  expect_equal(result$P_alpha, expected_result$P_alpha, tolerance = 1e-6)
})

test_that("Results proba theoretical second excursion, specific case (a = 0) ", {
  #Validated by simulation
  expected_result = list()
  expected_result$proba_q_i_geq_a = 1
  expected_result$P_alpha =  c(K = 1, L = 1, M = 1)
  result <- proba_theoretical_ith_excursion_markov(a = 0,,
                                                   theta = c("K", "L", "M"),
                                                   lambda = matrix(c(0.3, 0.2, 0.5, 0.4, 0.4, 0.2, 0.3, 0.3, 0.4), ncol = 3, byrow = TRUE),
                                                   score_function = c(-3,-1,2),
                                                   i = 2)
  expect_equal(result$proba_q_i_geq_a, expected_result$proba_q_i_geq_a, tolerance = 1e-6)
  expect_equal(result$P_alpha, expected_result$P_alpha, tolerance = 1e-6)
})

test_that("Results proba theoretical second excursion, specific case (0 in the score function) ", {
  #Validated by simulation
  expected_result = list()
  expected_result$proba_q_i_geq_a = 0.2816628
  expected_result$P_alpha = c(K = 0.3691211, L = 0.2132907, M = 0.3155178, N = 0.2312518, O = 0.3066551)
  result <- proba_theoretical_ith_excursion_markov(a = 4,,
                                                   theta = c("K", "L", "M","N","O"),
                                                   lambda = matrix(c(0.1, 0.7, 0.05, 0.05, 0.1,
                                                                     0.3, 0, 0.2, 0.15, 0.35,
                                                                     0.1, 0.4, 0.15, 0.2, 0.15,
                                                                     0.3, 0.05, 0.25, 0.3, 0.1,
                                                                     0.25, 0.45, 0.1, 0.15, 0.05), ncol = 5, nrow = 5, byrow = TRUE),
                                                   score_function = c(-5, -3, 0, 1, 4),
                                                   i = 2)
  expect_equal(result$proba_q_i_geq_a, expected_result$proba_q_i_geq_a, tolerance = 1e-6)
  expect_equal(result$P_alpha, expected_result$P_alpha, tolerance = 1e-6)
})

test_that("Results proba theoretical second excursion, specific case (diagonal of the transition matrix is filled by 0) ", {
  #Validated by simulation
  expected_result = list()
  expected_result$proba_q_i_geq_a = 0.1254411
  expected_result$P_alpha = c(K = 0.1237689, L = 0.1302016, M = 0.1227621)
  result <- proba_theoretical_ith_excursion_markov(a = 5,
                                                   theta = c("K", "L", "M"),
                                                   lambda = matrix(c(0, 0.3, 0.7, 0.4, 0, 0.6, 0.4, 0.6, 0), ncol = 3, byrow = TRUE),
                                                   score_function = c(-2,-1,2),
                                                   i = 2)
  expect_equal(result$proba_q_i_geq_a, expected_result$proba_q_i_geq_a, tolerance = 1e-6)
  expect_equal(result$P_alpha, expected_result$P_alpha, tolerance = 1e-6)
})

test_that("Result that differ from matrix_M", {
  expected_result = matrix(c(1,0,
                           1,0),
                           nrow = 2,
                           byrow = TRUE)

  result <- matrix_M(theta = c("K", "L"),
                            lambda = matrix(c(0.5, 0.5, 0.2, 0.8), ncol = 2, byrow = TRUE),
                            score_function = 2 * c(-4,1),
                            N_iter = 1000,
                            epsilon = 1e-16)
  expect_equal(result, expected_result)
})

test_that("Result where two positive scores are equal", {
  expected_result = matrix(c(0.2192127, 0.72629245, 0.05449489, 0, 0,
                           0.4653842, 0.32935851, 0.20525729, 0, 0,
                           0.3845969, 0.45565097, 0.15975218, 0, 0,
                           0.6615238, 0.08277737, 0.25569883, 0, 0,
                           0.4153842, 0.07935851, 0.50525729, 0, 0),
                         nrow = 5,
                         byrow = TRUE)
  
  result <- matrix_M(theta = c("a","b","c","d","e"),
                     lambda = matrix(c(0.1,0.7,0.05,0.05,0.1, 0.3,0.3,0.2,0.15,0.05, 0.1,0.4,0.15,0.2,0.15, 0.5,0.05,0.25,0.1,0.1, 0.25,0.05,0.5,0.15,0.05), ncol = 5,nrow=5, byrow = TRUE),
                     score_function = c(-10, -2, -1, 6, 6),
                     N_iter = 1000,
                     epsilon = 1e-16)
  expect_equal(result, expected_result, tolerance = 1e-06)
})

test_that("Results proba theoretical fourth excursion, error case (all scores are negative) ", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5,
                                                      theta = c("K", "L", "M"),
                                                      lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE),
                                                      score_function = c(-2, -1, -2),
                                                      i = 4),"[Invalid Input] score_function should have at least one positive integer.",fixed = TRUE)
})

test_that("Results proba theoretical fourth excursion, error case (theta have a non-unique value) ", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5,
                                                      theta = c("K", "L", "L"),
                                                      lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE),
                                                      score_function = c(-2, -1, 2),
                                                      i = 4),"[Invalid Input] theta should contains only unique values",fixed = TRUE)
})

test_that("Results proba theoretical fourth excursion, error case (theta have a different number of values) ", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5,
                                                      theta = c("K", "L"),
                                                      lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE),
                                                      score_function = c(-2, -1, 2),
                                                      i = 4),"[Invalid Input] lambda should be the same size as theta.",fixed = TRUE)
})

test_that("Results proba theoretical second excursion, error case (score expectation is positive)",{
  expect_error(proba_theoretical_ith_excursion_markov(a = 5,
                                                      theta = c("K", "L", "M"),
                                                      lambda = matrix(c(0.5, 0.1, 0.4, 0.4, 0.05, 0.55, 0.3, 0.4, 0.3), ncol = 3, byrow = TRUE),
                                                      score_function = c(2, -1, 2),
                                                      i = 2),"[Invalid Input] score expectation should be negative.",fixed = TRUE)
})

test_that("Results proba theoretical fourth excursion, error case (score are not integers) ", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5,
                                                      theta = c("K","L","M"),
                                                      lambda = matrix(c(0.1, 0.5, 0.4, 0.4, 0.05, 0.55, 0.3, 0.4, 0.3), ncol = 3, byrow = TRUE),
                                                      score_function = c(2.4,-1.8,-2),
                                                      i = 3),"[Invalid Input] scores of score_function should be integers.",fixed = TRUE)
})

test_that("Results proba theoretical fifth excursion, error case (i = 0) ", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5,
                                                      theta = c("K", "L", "M"),
                                                      lambda = matrix(c(0.1, 0.5, 0.4, 0.4, 0.05, 0.55, 0.3, 0.4, 0.3), ncol = 3, byrow = TRUE),
                                                      score_function = c(2,-1,-2),
                                                      i = 0),"[Invalid Input] i should be a strictly positive integer.",fixed = TRUE)
})

test_that("Results proba theoretical fourth excursion, error case (lambda is not a stochastic matrix) ", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5,
                                                      theta = c("K","L","M"),
                                                      lambda = matrix(c(0.1, 0.5, 0.5, 0.4, 0.01, 0.05, 0.3, 0.4, 0.3), ncol = 3, byrow = TRUE),
                                                      score_function = c(-4, -1, 2),
                                                      i = 1),"[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.). Sum of line 1 equal 1.100000",fixed = TRUE)
})
##############

test_that("Results proba theoretical first execution special case only positive scores", {
  result = proba_theoretical_first_excursion_markov(a = 5, 
                                                    theta = c("K","L","M"), 
                                                    lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), 
                                                    score_function = c(2,2,2),
                                                    prob0 = c(0.4444444, 0.2962963, 0.2592593) )
  expected_results = list()
  expected_results$s = 1
  expected_results$z = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
  
  expect_equal(result, expected_results$s)
  #  expect_equal(as.vector(result$z), expected_results$z)
})

test_that("Results proba theoretical first execution special case only negative scores", {
  result = proba_theoretical_first_excursion_markov(a = 5, 
                                                    theta = c("K","L","M"), 
                                                    lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), 
                                                    score_function = c(-2,-2,-2),
                                                    prob0 = c(0.4444444, 0.2962963, 0.2592593) )
  expected_results = list()
  expected_results$s = 0
  expected_results$z = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  expect_equal(result, expected_results$s)
  #  expect_equal(as.vector(result$z), expected_results$z)
})

test_that("Results proba theoretical ith execution normal case", {
  # Note : simulating with 1.000.000 iterations, we found 0.033892. 
  # With this simple score scheme, could it be calculated exactly ?
  result = proba_theoretical_ith_excursion_markov(a = 5, 
                                                  theta = c("K","L","M"), 
                                                  lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), 
                                                  score_function = c(-1,0,1),
                                                  prob0 = c(0.4444444, 0.2962963, 0.2592593), i = 3 )
  expected_results = list()
  expected_results$proba_q_i_geq_a = 0.03392221 # 0.034025 
  #expected_results$proba_q_i_geq_a =  0.0172477 #XXX je ne sais pas laquelle des valeurs est bonne entre l'initial (0.034025) et celle actuelle (0.0172)
  
  expect_equal(result$proba_q_i_geq_a, expected_results$proba_q_i_geq_a, tolerance = 0.000001)
})

test_that("Results proba theoretical ith execution special case only positive scores", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5, 
                                                      theta = c("K","L","M"), 
                                                      lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), 
                                                      score_function = c(2,2,2),
                                                      prob0 = c(0.4444444, 0.2962963, 0.2592593), i = 3 ),
               "[Invalid Input] score expectation should be negative.", fixed = TRUE
  )
})

test_that("Results proba theoretical ith execution special case only negative scores", {
  expect_error(proba_theoretical_ith_excursion_markov(a = 5, 
                                                      theta = c("K","L","M"), 
                                                      lambda = matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE), 
                                                      score_function = c(-2,-2,-2),
                                                      prob0 = c(0.4444444, 0.2962963, 0.2592593), i = 3 ),
               "[Invalid Input] score_function should have at least one positive integer.", fixed = TRUE
  )
})

test_that("Results matrice M", {
  result = matrix_M(c("K","L","M"),
                   matrix(c(0.5, 0.3, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4, 0.2), ncol = 3, byrow = TRUE),
                    c(-2,-1,2),
                    N_iter = 100)
  expected_results = matrix(c(0.6525332, 0.3474668,    0,
                              0.7050664, 0.2949336,    0,
                              0.5525332, 0.4474668,    0), ncol = 3, byrow = TRUE)
  #(theta, a, lambda, score_function, N_iter)
  expect_equal(expected_results, result, tolerance = 1e-4)
})

