context("Test of R function daudin")

test_that("Example from the R documentation of daudin(.)", {
  p <- daudin(local_score = 4, sequence_length = 50, 
              score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), sequence_min = -3, sequence_max = 2)
  expect_equal(p, 0.5700479,tolerance=1e-7)
})

test_that("Case 1 from the issue #7326 on the forge DGA", {
  prob <- c(0.951855990, 0.003892646, 0.000000000, 0.004302711, 0.039948653)
  p <- daudin(local_score=4, sequence_length=39, sequence_min=-2, sequence_max=2, score_probabilities=prob)
  expect_equal(p, 0.05712747)
})

test_that("Full cases from the issue #7326 on the forge DGA", {
  daudin.dat <- structure(list(local_score = c(29L, 15L, 4L, 12L, 12L, 12L, 4L, 
                                               5L, 4L, 14L, 4L, 18L, 4L, 4L, 11L, 5L, 15L, 12L, 4L, 4L, 4L, 
                                               16L, 4L, 4L, 4L, 4L, 4L, 5L, 4L, 4L, 11L, 4L, 4L, 16L, 5L, 5L, 
                                               5L, 5L, 5L, 4L, 4L, 4L, 5L, 15L, 5L, 52L, 4L, 4L, 5L, 4L, 20L, 
                                               5L, 5L, 34L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 12L, 4L, 24L, 
                                               4L, 4L, 4L, 16L, 4L, 12L, 20L, 4L, 26L, 4L, 4L, 4L, 19L, 3L, 
                                               3L, 4L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 3L, 3L, 3L, 15L),
                                debut = c(19L, 
                                          2L, 38L, 6L, 1L, 5L, 17L, 2L, 36L, 2L, 2L, 2L, 25L, 1L, 2L, 3L, 
                                          11L, 11L, 4L, 1L, 1L, 2L, 1L, 4L, 1L, 1L, 24L, 12L, 1L, 1L, 1L, 
                                          10L, 16L, 1L, 1L, 4L, 15L, 1L, 26L, 29L, 11L, 13L, 36L, 1L, 1L, 
                                          1L, 19L, 33L, 3L, 26L, 6L, 1L, 30L, 1L, 35L, 2L, 1L, 5L, 32L, 
                                          2L, 17L, 11L, 16L, 1L, 1L, 1L, 2L, 1L, 4L, 3L, 1L, 5L, 1L, 2L, 
                                          7L, 9L, 3L, 7L, 13L, 31L, 33L, 32L, 9L, 2L, 23L, 9L, 5L, 34L, 
                                          14L, 18L, 30L, 8L, 25L, 1L), 
                                fin = c(36L, 18L, 39L, 15L, 6L, 
                                        10L, 18L, 4L, 37L, 19L, 3L, 12L, 26L, 2L, 7L, 5L, 30L, 23L, 5L, 
                                        2L, 2L, 9L, 2L, 5L, 2L, 2L, 25L, 14L, 2L, 2L, 6L, 11L, 17L, 8L, 
                                        3L, 6L, 17L, 3L, 28L, 30L, 12L, 14L, 38L, 8L, 3L, 39L, 20L, 34L, 
                                        5L, 27L, 19L, 3L, 34L, 22L, 36L, 3L, 2L, 6L, 33L, 3L, 18L, 12L, 
                                        17L, 6L, 2L, 22L, 3L, 2L, 5L, 14L, 2L, 10L, 12L, 3L, 19L, 10L, 
                                        4L, 8L, 35L, 32L, 34L, 33L, 10L, 3L, 24L, 10L, 6L, 35L, 15L, 
                                        19L, 31L, 9L, 26L, 12L), 
                                longeur_sequence = c(39L, 39L, 39L, 
                                                     39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 
                                                     39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 
                                                     39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 
                                                     39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 
                                                     39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 
                                                     39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 
                                                     39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L, 39L
                                ), 
                                pvalue = c(NA, NA, NA, NA, NA, NA, NA, NA, 3.161659, NA, NA, 
                                           NA, 17.29334, NA, NA, NA, NA, NA, 3.161659, NA, NA, NA, NA, NA, 
                                           NA, NA, NA, NA, NA, 79.22113, NA, NA, NA, NA, NA, NA, NA, NA, 
                                           2.505567, NA, 17.29334, NA, NA, NA, NA, 1.578604e+142, NA, NA, 
                                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                                           1.25114, NA, NA, 2631549, NA, NA, 17.29334, NA, NA, NA, NA, NA, 
                                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                                           NA, NA)), 
                           row.names = c(NA, -94L), class = "data.frame")
  #N.B. Those expected p-values are issue from the daudin() function itself in the version 1.0.8. So there is no warranty they are correct. 
  #This case is only to verify the bug #7326 is corrected. Eventually, those expected p-values are to be modified in case of change in daudin function.
  expectedPv = c(1.1433020784519e-19, 4.80416973117856e-10, 0.0571274669482736, 
                 1.59281651368091e-07, 1.59281651368091e-07, 1.59281651368091e-07, 
                 0.0571274669482736, 0.00312004467287479, 0.0571274669482736, 
                 6.46423153947766e-09, 0.0571274669482736, 1.06247092161754e-11, 
                 0.0571274669482736, 0.0571274669482736, 2.59662408482645e-07, 
                 0.00312004467287479, 4.80416973117856e-10, 1.59281651368091e-07, 
                 0.0571274669482736, 0.0571274669482736, 0.0571274669482736, 2.62168444926764e-10, 
                 0.0571274669482736, 0.0571274669482736, 0.0571274669482736, 0.0571274669482736, 
                 0.0571274669482736, 0.00312004467287479, 0.0571274669482736, 
                 0.0571274669482736, 2.59662408482645e-07, 0.0571274669482736, 
                 0.0571274669482736, 2.62168444926764e-10, 0.00312004467287479, 
                 0.00312004467287479, 0.00312004467287479, 0.00312004467287479, 
                 0.00312004467287479, 0.0571274669482736, 0.0571274669482736, 
                 0.0571274669482736, 0.00312004467287479, 4.80416973117856e-10, 
                 0.00312004467287479, 1.70509626232842e-35, 0.0571274669482736, 
                 0.0571274669482736, 0.00312004467287479, 0.0571274669482736, 
                 4.30212967482166e-13, 0.00312004467287479, 0.00312004467287479, 
                 7.41397840765337e-23, 0.0571274669482736, 0.0571274669482736, 
                 0.0571274669482736, 0.0571274669482736, 0.0571274669482736, 0.0571274669482736, 
                 0.0571274669482736, 0.0571274669482736, 0.0571274669482736, 1.59281651368091e-07, 
                 0.0571274669482736, 7.03248132206193e-16, 0.0571274669482736, 
                 0.0571274669482736, 0.0571274669482736, 2.62168444926764e-10, 
                 0.0571274669482736, 1.59281651368091e-07, 4.30212967482166e-13, 
                 0.0571274669482736, 2.83827416760602e-17, 0.0571274669482736, 
                 0.0571274669482736, 0.0571274669482736, 8.72382977649776e-13, 
                 0.0686244605146422, 0.0686244605146422, 0.0571274669482736, 0.0686244605146422, 
                 0.0686244605146422, 0.0686244605146422, 0.0686244605146422, 0.0686244605146422, 
                 0.0686244605146422, 0.0686244605146422, 0.0571274669482736, 0.0686244605146422, 
                 0.0686244605146422, 0.0686244605146422, 4.80416973117856e-10)
  prob <- c(0.951855990, 0.003892646, 0.000000000, 0.004302711, 0.039948653)
  n <- dim(daudin.dat)[1]
  p <- rep(NA,n)
  for (i in 1:n) {
    p[i] <- daudin(local_score=daudin.dat$local_score[i], sequence_length=daudin.dat$longeur_sequence[i], sequence_min=-2, sequence_max=2, score_probabilities=prob)
  }
  # Developper note : the '!!' is to print the value of 'i' in the output (useful to debug failed tests)
  for (i in 1:n) {
    expect_equal(p[!!i],expectedPv[!!i])
  }
})

## Input Error testing
test_that("Entry probability does not sum to 1.0", {
  prob <- c(0.6,0.1,0.1,0.1) # sum to 0.9 (error)
  expect_error(daudin(150, 10000, prob, -2, 1), "[Invalid Input] score_probabilities must sum to 1.0.", fixed = TRUE)
})

############## Link to ticket #1 on forge MIA/INRAE #######################
test_that("Alternative and equivalent parameters specification", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  p <- daudin(150, 10000, prob, sequence_min = -5, sequence_max =  5)
  p0 <- daudin(150, 10000, prob, score_values = -5:5)
  p1 <- daudin(150, 10000, prob, score_values = -5:5, sequence_min = -5, sequence_max =  5)
  
  expect_equal(p, p0)
})

test_that("Alternative and equivalent parameters specification (with 'hole')", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- prob[prob != 0]
  score_values <- which(prob != 0) - 6
  
  p <- daudin(150, 10000, prob, sequence_min = -5, sequence_max =  5)
  p0 <- daudin(150, 10000, prob0, score_values = score_values)
  expect_equal(p, p0)
})

test_that("Alternative and equivalent parameters specification (with leading et trailing 0)", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  probShort <- prob[prob != 0]
  score_values <- which(prob != 0) - 6
  probShort0 <- c(0, probShort)
  score_values0 <- c(-6, score_values)
  probShort00 <- c(probShort0, 0)
  score_values00 <- c(score_values0 , 6)
  
  
  p <- daudin(150, 10000, prob, sequence_min = -5, sequence_max =  5)
  p1 <- daudin(150, 10000, probShort, score_values = score_values)
  p0 <- daudin(150, 10000, probShort0, score_values = score_values0)
  p00 <- daudin(150, 10000, probShort00, score_values = score_values00)
  
  expect_equal(p, p1)
  expect_equal(p, p0)
  expect_equal(p, p00)
})


############## Link to ticket #10382 on forge DGA #######################
test_that("Leading zero probability element does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,prob)
  p <- daudin(150, 10000, prob, -5, 5)
  p0 <- daudin(150, 10000, prob0, -5-1, 5)
  expect_equal(p, p0)
})

test_that("Leading multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,0,0,0,prob)
  p <- daudin(150, 10000, prob, -5, 5)
  p0 <- daudin(150, 10000, prob0, -5-4, 5)
  expect_equal(p, p0)
})

test_that("Trailing zero probability element does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(prob,0)
  p <- daudin(150, 10000, prob, -5, 5)
  p0 <- daudin(150, 10000, prob0, -5, 5+1)
  expect_equal(p, p0)
})

test_that("Trailing multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(prob, 0,0,0,0)
  p <- daudin(150, 10000, prob, -5, 5)
  p0 <- daudin(150, 10000, prob0, -5, 5+4)
  expect_equal(p, p0)
})

test_that("Leading and trailing multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,0,0,0,prob,0,0)
  p <- daudin(150, 10000, prob, -5, 5)
  p0 <- daudin(150, 10000, prob0, -5-4, 5+2)
  expect_equal(p, p0)
})
###################################################################################

test_that("Zero probability elements does not change the results", {
  # En pratique, on multiplie par deux tous les scores ; théoriquement, on doit trouver la même p-value
  prob <- c(0.6, 0.1, 0.1, 0.2)
  scoremin <- -1
  scoremax <- scoremin + length(prob) -1 #2
  scoremin2 <- scoremin*2 #-2
  scoremax2 <- scoremax*2 #4
  prob0 <- c(0.6, 0.0, 0.1, 0.0, 0.1, 0.0, 0.2)
  p <- daudin(80, 10000, prob, scoremin, scoremax)
  p0 <- daudin(80*2, 10000, prob0, scoremin2, scoremax2)
  expect_equal(p, p0)
})
