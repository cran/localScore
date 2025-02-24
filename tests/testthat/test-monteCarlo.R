test_that("Parameter checks", {
  #link to bug #11238 (forge DGA)
  mySeq <- sample(-7:6, replace = TRUE, size = 100000) 
  expect_error(karlinMonteCarlo(local_score = 66, sequence_length = 100000,  
                                FUN = function(x, simulated_sequence_length) {
                                  return(sample(x = x, 
                                                size = simulated_sequence_length,
                                                replace = TRUE))
                                }, 
                                x = mySeq, 
                                simulated_sequence_length = 1000,
                                numSim = 10) ,
               "[Invalid Input] Parameters of FUN should differ from parameters of karlinMonteCarlo function", fixed = TRUE)
  expect_error(monteCarlo(local_score = 50.5,
                          FUN = function(numSim) {return(sample(x = x, 
                                                         size = numSim, replace = TRUE))},
                                                         x = mySeq, numSim = 1000) ,
               "[Invalid Input] Parameters of FUN should differ from parameters of monteCarlo function", fixed = TRUE)
})
