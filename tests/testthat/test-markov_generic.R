context("Test of generic functions around Markov models")

test_that("Stationary distribution", {
  
  # Following is a list of transition matrix and associated stationnary distribution with "no particularities"
  p_all = list(
    list(p = matrix(c(0.6, 0.1, 0.3,
                  0.1, 0.7, 0.2,
                  0.2, 0.2, 0.6), byrow=TRUE, ncol=3),
         mu = c(0.2758621, 0.3448276, 0.3793103)),
    list(p = matrix(c(0.0, 0.7, 0.3,
                  0.1, 0.7, 0.2,
                  0.2, 0.2, 0.6), byrow=TRUE, ncol=3),
         mu= c(0.1230769, 0.5230769, 0.3538462)),
    list(p = matrix(c(0.6, 0.1, 0.3,
                      0.1, 0.7, 0.2,
                      0.5, 0.5, 0.0), byrow=TRUE, ncol=3),
         mu = c(0.3571429, 0.4464286, 0.1964286)),
    list(p = matrix(c(0.6, 0.4, 0.0,
                      0.1, 0.7, 0.2,
                      0.0, 0.4, 0.6), byrow=TRUE, ncol=3),
         mu = c(0.1428571, 0.5714286, 0.2857143)),
    list(p = matrix(c(0.6, 0.4, 0.0,
                      0.3, 0.0, 0.7,
                      0.0, 0.5, 0.5), byrow=TRUE, ncol=3),
         mu = c(0.2380952, 0.3174603, 0.4444444)),
    list(p = matrix(c(0.0, 0.4, 0.6,
                      0.3, 0.0, 0.7,
                      0.5, 0.5, 0.0), byrow=TRUE, ncol=3),
         mu = c(0.2914798, 0.3139013, 0.3946188)),
    # next matrix : state '1' is not reacheable in finite time
    list(p = matrix(c(0.0, 0.0, 1.0,
                      0.0, 0.7, 0.3,
                      0.0, 0.5, 0.5), byrow=TRUE, ncol=3),
         mu = c(0.000, 0.625, 0.375))
  )
  
  #Next function calculate stationnary distribution in R code (without any checking and assuming a lot of unchecked hypothesis). 
  #It is better to not use it in test, but use the results hard-coded.
  cal_vp_R = function(ptrans) {
    dist_stat = eigen(t(ptrans), symmetric = FALSE)$vectors[,1]
    return(dist_stat/sum(dist_stat))
  }
  
  # Developper note : the '!!' is to print the value of 'i' in the output (useful to debug failed tests)
  for (i in 1:length(p_all)) {
    # Check that stationary distribution sums to 1.
    expect_equal(sum(stationary_distribution(p_all[[!!i]]$p)),1)
    # Check that stationary distribution is all positive valued (null included)
    expect_true(all(stationary_distribution(p_all[[!!i]]$p)>=0))
    # Check correct value of stationary distribution
    expect_equal(stationary_distribution(p_all[[!!i]]$p),p_all[[!!i]]$mu, tolerance = 1e-6, scale = 1)
  }

    # Check error on non stochastic matrix
  p1 = matrix(c(5.2, 1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0), byrow=TRUE, ncol=3)
  expect_error(stationary_distribution(p1),"[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.). Sum of line 1 equal 6.200000", fixed = TRUE)
  
  # Check error/warning on non irreductible Markov matrix (XXX TODO, or not)
  # p2 : cyclic
  p2 = matrix(c(0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                1.0, 0.0, 0.0), byrow=TRUE, ncol=3)
  expect_equal(stationary_distribution(p2),c(0.3333333, 0.3333333, 0.3333333), tolerance = 1e-6, scale = 1)
  # p3 : absorbing state '2'
  p3 = matrix(c(0.2, 0.5, 0.3,
                0.0, 1.0, 0.0,
                0.1, 0.8, 0.1), byrow=TRUE, ncol=3)
  expect_equal(stationary_distribution(p3),c(0.0, 1.0, 0.0), tolerance = 1e-6, scale = 1)
})
