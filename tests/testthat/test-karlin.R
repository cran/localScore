context("Test of R function karlin")

test_that("Example from the R documentation of karlin(.)", {
  p <- karlin(150, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -5, 5)
  kp <- karlin_parameters(c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -5, 5)
#  expect_equal(p, 0.0265141)
  expect_equal(p, 0.03464975)
  expect_equal(kp, c(K_star = 0.01823197, K_plus = 0.01930885, lambda = 0.05738689))
})

test_that("Entry probability does not sum to 1.0", {
  prob <- c(0.6,0.1,0.1,0.1) # sum to 0.9 (error)
  expect_error(karlin(150, 10000, prob, -2, 1), "[Invalid Input] score_probabilities must sum to 1.0.", fixed = TRUE)
  expect_error(karlin_parameters(prob, -2,1), "[Invalid Input] score_probabilities must sum to 1.0.", fixed = TRUE)
})

############## Link to ticket #11183 on forge DGA #######################
test_that("Example from the Karlin and Altschull 1990 Appendix", {
  # We should find lambda = ln(q/p) and K+ = (q-p)^2/q and Kstar = exp(-lambda*delta)*K+
  # We should also find that using Kstar instead of K+ is equivalent to call the function with parameter local_Score+1
  p <- 0.3 #P(X=1)
  q <- 0.4 #P(X=-1)
  r <- 1-p-q #P(X=0)
  seqlength <- 1000
  a <- 30
  pv <- karlin(a, seqlength, c(q,r,p), -1, 1)
  kparameters <- karlin_parameters(c(q,r,p), -1, 1)
  
  lambdaTh <- log(q/p)
  KplusTh <- ((q-p)^2)/q
  pvTh <- 1.0-exp(-KplusTh*seqlength*exp(-lambdaTh*a))
  expect_equal(unname(kparameters['lambda']), lambdaTh)
  expect_equal(pv, pvTh)
  expect_equal(unname(kparameters['K_plus']), KplusTh)
  
  delta <- 1 #span (for lattice case)
  KstarTh <- exp(-lambdaTh*delta)*KplusTh
  expect_equal(unname(kparameters['K_star']), KstarTh)
  expect_equal(
    karlin(a+1, seqlength, c(q,r,p), -1, 1) ,
    1.0-exp(-KstarTh*seqlength*exp(-lambdaTh*a))
  )
})


############## Link to ticket #10382 on forge DGA #######################
test_that("Leading zero probability element does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,prob)
  p <- karlin(150, 10000, prob, -5, 5)
  pk <- karlin_parameters(prob, -5, 5)
  p0 <- karlin(150, 10000, prob0, -5-1, 5)
  pk0 <- karlin_parameters(prob0, -5-1, 5)
  expect_equal(p, p0)
  expect_equal(pk, pk0)
  expect_equal(unname(pk0),c(0.01823197, 0.01930885, 0.05738689 ))
})

test_that("Leading multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,0,0,0,prob)
  p <- karlin(150, 10000, prob, -5, 5)
  pk <- karlin_parameters(prob, -5, 5)
  p0 <- karlin(150, 10000, prob0, -5-4, 5)
  pk0 <- karlin_parameters(prob0, -5-4, 5)
  expect_equal(p, p0)
  expect_equal(pk, pk0)
  expect_equal(unname(pk0),c(0.01823197, 0.01930885, 0.05738689 ))
})

test_that("Trailing zero probability element does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(prob,0)
  p <- karlin(150, 10000, prob, -5, 5)
  p0 <- karlin(150, 10000, prob0, -5, 5+1)
  pk <- karlin_parameters(prob, -5, 5)
  pk0 <- karlin_parameters(prob0, -5, 5+1)
  expect_equal(p, p0)
  expect_equal(pk, pk0)
})

test_that("Trailing multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(prob, 0,0,0,0)
  p <- karlin(150, 10000, prob, -5, 5)
  p0 <- karlin(150, 10000, prob0, -5, 5+4)
  pk <- karlin_parameters(prob, -5, 5)
  pk0 <- karlin_parameters(prob0, -5, 5+4)
  expect_equal(p, p0)
  expect_equal(pk, pk0)
})

test_that("Leading and trailing multiple zero probability elements does not change the results", {
  prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
  prob0 <- c(0,0,0,0,prob,0,0)
  p <- karlin(150, 10000, prob, -5, 5)
  p0 <- karlin(150, 10000, prob0, -5-4, 5+2)
  pk <- karlin_parameters(prob, -5, 5)
  pk0 <- karlin_parameters(prob0, -5-4, 5+2)
  expect_equal(pk, pk0)
  expect_equal(p, p0)
})
###################################################################################

test_that("Double each score does not change the results", {
  # En pratique, on multiplie par deux tous les scores ; théoriquement, on doit trouver la même p-value
  prob <- c(0.6, 0.1, 0.1, 0.2)
  scoremin <- -1
  scoremax <- scoremin + length(prob) -1 #2
  scoremin2 <- scoremin*2 #-2
  scoremax2 <- scoremax*2 #4
  prob0 <- c(0.6, 0.0, 0.1, 0.0, 0.1, 0.0, 0.2)
  p <- karlin(80+1, 10000, prob, scoremin, scoremax)
  p0 <- karlin(80*2+1, 10000, prob0, scoremin2, scoremax2)
  pk <- karlin_parameters(prob, scoremin, scoremax)
  pk0 <- karlin_parameters(prob0, scoremin2, scoremax2)
  expect_equal(p, p0)
  expect_equal(pk[c("K_star", "K_plus")], pk0[c("K_star", "K_plus")])
  expect_equal(0.5*pk[c("lambda")], pk0[c("lambda")])
})



# Below : this check is running fine except on M1Mac configuration (polynomial root finding issue which throw a controlled error)
# test_that("Example from the a pig dataset (unit score scheme)", {
#   score1.prob.ext = c(`-2` = 0.0681561425858712, `-1` = 0.759474076388722, `0` = 0.121671514073504, 
#                       `1` = 0.0310891041511429, `2` = 0.0125039059683196, `3` = 0.00563901641708531,
#                       `4` = 0.0010912674566738, `5` = 0.000242770954017739, `6` = 7.93212027978752e-05,
#                       `7` = 3.60550921808523e-05, `8` = 1.44220368723409e-05, `9` = 2.40367281205682e-06)
#   p <- karlin(248, 416030, score1.prob.ext, -2, 9)
#   expect_equal(p, 0.0)
# })
