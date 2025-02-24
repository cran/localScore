context("Test of R function mcc")

test_that("Example 1 from the R documentation of mcc(.)", {
  p <- mcc(40, 100, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
  expect_equal(p, 0.003537593)
})

test_that("Example 2 from the R documentation of mcc(.)", {
  p <- mcc(40, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
  expect_equal(p, 0.292068)
})

## Check of input error
test_that("Entry probability does not sum to 1.0", {
  prob <- c(0.6,0.1,0.1,0.1) # sum to 0.9 (error)
  expect_error(mcc(15, 1000, prob, -2, 1), "[Invalid Input] score_probabilities must sum to 1.0.", fixed = TRUE)
})

############## Link to ticket #10382 on forge DGA #######################
test_that("Leading zero probability element does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,prob)
  p <- mcc(150, 10000, prob, -5, 5)
  p0 <- mcc(150, 10000, prob0, -5-1, 5)
  expect_equal(p, p0)
})

test_that("Leading multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,0,0,0,prob)
  p <- mcc(150, 10000, prob, -5, 5)
  p0 <- mcc(150, 10000, prob0, -5-4, 5)
  expect_equal(p, p0)
})

test_that("Trailing zero probability element does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(prob,0)
  p <- mcc(150, 10000, prob, -5, 5)
  p0 <- mcc(150, 10000, prob0, -5, 5+1)
  expect_equal(p, p0)
})

test_that("Trailing multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(prob, 0,0,0,0)
  p <- mcc(150, 10000, prob, -5, 5)
  p0 <- mcc(150, 10000, prob0, -5, 5+4)
  expect_equal(p, p0)
})

test_that("Leading and trailing multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,0,0,0,prob,0,0)
  p <- mcc(150, 10000, prob, -5, 5)
  p0 <- mcc(150, 10000, prob0, -5-4, 5+2)
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
  p <- mcc(80, 10000, prob, scoremin, scoremax)
  p0 <- mcc(80*2, 10000, prob0, scoremin2, scoremax2)
  expect_equal(p, p0)
})


# Below : this check is running fine except on M1Mac configuration (polynomial root finding issue which throw a controlled error)
# test_that("Example from the a pig dataset (unit score scheme)", {
#   score1.prob.ext = c(`-2` = 0.0681561425858712, `-1` = 0.759474076388722, `0` = 0.121671514073504, 
#                       `1` = 0.0310891041511429, `2` = 0.0125039059683196, `3` = 0.00563901641708531,
#                       `4` = 0.0010912674566738, `5` = 0.000242770954017739, `6` = 7.93212027978752e-05,
#                       `7` = 3.60550921808523e-05, `8` = 1.44220368723409e-05, `9` = 2.40367281205682e-06)
#   p <- mcc(248, 416030, score1.prob.ext, -2, 9)
#   expect_equal(p, 0.0)
# })
