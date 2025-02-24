context("Tests of diverse tool functions")

test_that("Markov score generation (approximate testing because random)", {
  nState <- 5
  MyTransMat <- matrix(c(0.3,0.1,0.1,0.1,0.4,
                         0.2,0.2,0.1,0.2,0.3, 
                         0.3,0.4,0.1,0.1,0.1, 
                         0.3,0.3,0.1,0.0,0.3,
                         0.1,0.1,0.2,0.3,0.3), ncol = nState, byrow = TRUE)
  sequence <- transmatrix2sequence(matrix = MyTransMat,length = 100000, score = c(-2,-1,0,2,3))
  count <- table(sequence[-length(sequence)],sequence[-1])
  count <- as.matrix(addmargins(count))
  relFreq <- count[1:5,1:5]/count[1:5,6]
  expect_true(all(abs(relFreq - MyTransMat) < 0.01))
})

test_that("Markov score generation (good output type)", {
  B <-  matrix(c(0.2, 0.8, 0.4, 0.6), nrow = 2, byrow = TRUE)
  s1 <- transmatrix2sequence(B, length = 10)
  s2 <- transmatrix2sequence(B, length = 10, score = c(-2,1))
  s3 <- transmatrix2sequence(B, length = 10, score = c("A","B"))
  s4 <- transmatrix2sequence(B, length = 10, score = c(-2.5,1))
  expect_true(is.integer(s1))
  expect_true(is.integer(s2))
  expect_false(is.integer(s3))
  expect_false(is.integer(s4))
  expect_true(is.numeric(s4))
})


test_that("Results of CharSequence2ScoreSequence", {
  expected_results = c(2, 4, -1, 5, -1, -1, -1, 3, 0, 3, 4, 4, 2, 2, 4, -1, 5, -1, -1, 4, 4, 3, 5, 0, 4, -4, -4, 5, -5, 4, 5)
  data("Seq31")
  data("Seq219")
  data("HydroScore")
  result <- CharSequence2ScoreSequence(Seq31, HydroScore)
  expect_equal(result, expected_results)
  
  expected_results = c(2, -1, 0, 4, -1, 0, -2, -2, 2, -5, -5, 0, -2, 3, -2, 4, 2, 4, 4, 4, 4, 3, 4, 4, 0, -2, -5, 4, 4, 4, 2, 5, -1, 3, -3,
                       4, -2, 5, -4, -1, -5, -4, 3, 4, -5, -4, -4, 5, -3, -4, -4, 4, 4, 4, -1, 0, 2, -1, -4, 5, -1, -4, -4, -1, 0, 0, 2, 0, 0,
                       4, -5, -1, -3, 4, -4, 5, -1, -4, -1, 2, 0, -3, 5, 4, -1, -1, -4, -4, -4, 2, -1, -4, 0, -4, 3, 2, 3, -1, -1, -4, -4, -1,
                       -4, 2, 3, -4, 4, 3, 3, -4, -1, -4, 0, -1, 0, -5, 5, -2, -4, -4, 4, 4, 5, 4, -4, 2, -4, -3, 0, 4, -4, 2, -4, -4, -1, -4,
                       -4, 5, 2, -4, 4, -4, -4, 4, -4, -2, 4, -4, 4, -4, 4, -5, -5, 4, -4, -4, 4, -1, -4, -1, 5, 4, -4, -4, 3, 2, -1, 2, -4, -4,
                       -5, -4, -4, -4, 2, -5, -4, -1, -4, -4, -1, -1, -4, -1, -5, 4, 4, -1, 3, -1, 5, 3, -1, 2, 3, 3, 4, 5, 0, 4, 2, -1, -1, -4,
                       4, 3, -1, 4, -5, -5, 3, 3, -4, 2, -4, -4, 4, 5, -4)
  result <- CharSequence2ScoreSequence(Seq219, HydroScore)
  expect_equal(result, expected_results)
  
  expected_results = c(2, -1, 0, 4, -1, 0, -2, -2, 2, -5, -5, 0, -2, 3, -2, 4, 2, 4, 4, 4, 4, 3, 4, 4, 0, -2, -5, 4, 4, 4, 2, 5, -1, 3, -3,
                       4, -2, 5, -4, -1, -5, -4, 3, 4, -5, -4, -4, 5, -3, -4, -4, 4, 4, 4, -1, 0, 2, -1, -4, 5, -1, -4, -4, -1, 0, 0, 2, 0, 0,
                       4, -5, -1, -3, 4, -4, 5, -1, -4, -1, 2, 0, -3, 5, 4, -1, -1, -4, -4, -4, 2, -1, -4, 0, -4, 3, 2, 3, -1, -1, -4, -4, -1,
                       -4, 2, 3, -4, 4, 3, 3, -4, -1, -4, 0, -1, 0, -5, 5, -2, -4, -4, 4, 4, 5, 4, -4, 2, -4, -3, 0, 4, -4, 2, -4, -4, -1, -4,
                       -4, 5, 2, -4, 4, -4, -4, 4, -4, -2, 4, -4, 4, -4, 4, -5, -5, 4, -4, -4, 4, -1, -4, -1, 5, 4, -4, -4, 3, 2, -1, 2, -4, -4,
                       -5, -4, -4, -4, 2, -5, -4, -1, -4, -4, -1, -1, -4, -1, -5, 4, 4, -1, 3, -1, 5, 3, -1, 2, 3, 3, 4, 5, 0, 4, 2, -1, -1, -4,
                       4, 3, -1, 4, -5, -5, 3, 3, -4, 2, -4, -4, 4, 5, -4)
  result <- CharSequence2ScoreSequence(Seq219, HydroScore)
  expect_equal(result, expected_results)
})


