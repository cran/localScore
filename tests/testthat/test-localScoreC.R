context("Score local detection and Record times")
#library(localScore)

test_that("localScoreC and recordTimes results checks", {
  ##This check 2 things :
  #  1. output of localScoreC is expected (all field)
  #  2. output of recordTimes is expected AND is equal to the output of localScoreC()$RecordTime[-1]
  #Helper function to build object return by localScoreC
  built_expected_results  <-  function(localscore, suboptimalSeg, RecordTime) {
    names(localscore) <- c("value", "begin", "end")
    suboptimalSeg <- as.data.frame(suboptimalSeg) # un truc rapide pour ne pas récrire tous les tests (qui utilisent matrix)
    colnames(suboptimalSeg) <- c("value", "begin", "end")
    resLS <- list(localScore = localscore, suboptimalSegmentScores = suboptimalSeg, RecordTime = RecordTime)
    resRecordFn <- RecordTime[-1]
    return(list(resLS = resLS, resRecordFn = resRecordFn))
  }
  # Helper function to call localScoreC() AND recordTimes()
  built_call <- function(score_sequence) {
    resLS <- localScoreC(score_sequence, suppressWarnings = TRUE)
    resRecordFn <- recordTimes(score_sequence)
    return(list(resLS = resLS, resRecordFn = resRecordFn))
  }
 
  # Test with all null values as input
  expect_equal(built_call(c(0,0,0,0,0)), built_expected_results(localscore = c(0,0,0), 
                                                                suboptimalSeg =  matrix(c(0,0,0), ncol = 3, byrow = TRUE),
                                                                RecordTime = 0))
  
  # Test with all null values as input
  expect_equal(built_call(c(0,0,0,0,0)), built_expected_results(localscore = c(0,0,0), 
                                                                suboptimalSeg =  matrix(c(0,0,0), ncol = 3, byrow = TRUE),
                                                                RecordTime = 0))
  # Test with all negative values as input
  expect_equal(built_call(c(-1,-1,-1,-1,-1)), built_expected_results(localscore = c(0,0,0), 
                                                                     suboptimalSeg =  matrix(c(0,0,0), ncol=3, byrow=TRUE),
                                                                     RecordTime = 0:5))
  # Test with all positive values as input
  expect_equal(built_call(c(1,2,3,4,5)), built_expected_results(localscore = c(15,1,5), 
                                                                suboptimalSeg =  matrix(c(15,1,5), ncol=3, byrow=TRUE),
                                                                RecordTime = 0))
  # Test local score in the middle and sequence beginning with a mountain
  expect_equal(built_call(c(1,-2,1,3,-5)), built_expected_results(localscore = c(4,3,4), 
                                                                  suboptimalSeg =  matrix(c(1,1,1,
                                                                                            4,3,4), ncol=3, byrow=TRUE),
                                                                  RecordTime = c(0,2,5)))
  # Test local score in the middle and sequence ending with a mountain
  expect_equal(built_call(c(-2,1,1,-5,1)), built_expected_results(localscore = c(2,2,3), 
                                                                  suboptimalSeg =  matrix(c(2,2,3,
                                                                                            1,5,5), ncol=3, byrow=TRUE),
                                                                  RecordTime = c(0,1,4)))
  #built_call(c(3,-2,1,1,-5)) (problème pour suboptimal : soit une seule montagne de 1 à 4 soit 2 montagnes [1-1] [3-4])
  # Pour le dire autrement (?) : la fonction localScoreC ne renvoie pas le "même" résultat sur cette séquence, et sa séquence prise dans l'autre sens.
  
  expect_equal(built_call(c(-1,3,-2,1,1,-5)), built_expected_results(localscore = c(3,2,2), 
                                                                     suboptimalSeg =  matrix(c(3,2,2), ncol=3, byrow=TRUE),
                                                                     RecordTime = c(0, 1, 6)))
  # lindley finishing in 0 (not below 0)
  expect_equal(built_call(c(-1,3,-2,1,1,-3)), built_expected_results(localscore = c(3,2,2), 
                                                                     suboptimalSeg =  matrix(c(3,2,2), ncol=3, byrow=TRUE),
                                                                     RecordTime = c(0, 1)))
  
  expect_equal(built_call(c(1,-2,1,1,-5)), built_expected_results(localscore = c(2,3,4), 
                                                                  suboptimalSeg =  matrix(c(1,1,1,
                                                                                            2,3,4), ncol=3, byrow=TRUE),
                                                                  RecordTime = c(0,2,5)))
  # Test two mountains in the middle
  expect_equal(built_call(c(-2,1,1,-5,1,-4)), built_expected_results(localscore = c(2,2,3), 
                                                                     suboptimalSeg =  matrix(c(2,2,3,
                                                                                               1,5,5), ncol=3, byrow=TRUE),
                                                                     RecordTime = c(0, 1, 4, 6)))
})

test_that("localScoreC_double and recordTimes results checks", {
  ##This check 2 things :
  #  1. output of localScoreC_double is expected (all field)
  #  2. output of recordTimes is expected AND is equal to the output of localScoreC_double()$RecordTime[-1]
  #Helper function to build object return by localScoreC_double

  built_expected_results_double  <-  function(localscore, suboptimalSeg, RecordTime) {
    names(localscore) <- c("value", "begin", "end")
    colnames(suboptimalSeg) <- c("value", "begin", "end")
    resLS <- list(localScore = localscore, suboptimalSegmentScores = suboptimalSeg, RecordTime = RecordTime)
    resRecordFn <- RecordTime[-1]
    return(list(resLS = resLS, resRecordFn = resRecordFn))
  }
  # Helper function to call localScoreC() AND recordTimes()
  built_call_double <- function(score_sequence) {
    resLS <- localScoreC_double(score_sequence, suppressWarnings = TRUE)
    resRecordFn <- recordTimes(score_sequence)
    return(list(resLS = resLS, resRecordFn = resRecordFn))
  }
  
  # Test with all null values as input
  expect_equal(built_call_double(0.1*c(0,0,0,0,0)), built_expected_results_double(localscore = c(0.0,0.0,0.0), 
                                                                suboptimalSeg =  data.frame(0,0,0),
                                                                RecordTime = 0))
  
  # Test with all negative values as input
  expect_equal(built_call_double(0.1*c(-1,-1,-1,-1,-1)), built_expected_results_double(localscore = c(0.0,0.0,0.0),  
                                                                                suboptimalSeg =  data.frame(0,0,0),
                                                                                RecordTime = 0:5))
  # Test with all positive values as input
  expect_equal(built_call_double(0.1*c(1,2,3,4,5)), built_expected_results_double(localscore = c(1.5,1,5), 
                                                                suboptimalSeg =  data.frame(1.5,1,5),
                                                                RecordTime = 0))
  # Test local score in the middle and sequence beginning with a mountain
  expect_equal(built_call_double(0.1*c(1,-2,1,3,-5)), built_expected_results_double(localscore = c(0.4,3,4), 
                                                                  suboptimalSeg =  data.frame(
                                                                    value = c(0.1,0.4),
                                                                    begin = c(1,3),
                                                                    end = c(1,4)),
                                                                  RecordTime = c(0,2,5)))
  # Test local score in the middle and sequence ending with a mountain
  expect_equal(built_call_double(0.1*c(-2,1,1,-5,1)), built_expected_results_double(localscore = c(0.2,2,3), 
                                                                  suboptimalSeg =  data.frame(
                                                                    value = c(0.2,0.1),
                                                                    begin = c(2,5),
                                                                    end = c(3,5)),
                                                                  RecordTime = c(0,1,4)))
  #built_call_double(0.1*c(3,-2,1,1,-5)) (problème pour suboptimal : soit une seule montagne de 1 à 4 soit 2 montagnes [1-1] [3-4])
  # Pour le dire autrement (?) : la fonction localScoreC ne renvoie pas le "même" résultat sur cette séquence, et sa séquence prise dans l'autre sens.
  
  expect_equal(built_call_double(0.1*c(-1,3,-2,1,1,-5)), built_expected_results_double(localscore = c(0.3,2,2), 
                                                                     suboptimalSeg =  data.frame(
                                                                       value = 0.3,
                                                                       begin = 2,
                                                                       end = 2),
                                                                     RecordTime = c(0, 1, 6)))
  # lindley finishing in 0 (not below 0)
  expect_equal(built_call_double(0.1*c(-1,3,-2,1,1,-3)), built_expected_results_double(localscore = c(0.3,2,2), 
                                                                                suboptimalSeg =  data.frame(
                                                                                  value = 0.3,
                                                                                  begin = 2,
                                                                                  end = 2),
                                                                                RecordTime = c(0, 1)))
  
  expect_equal(built_call_double(0.1*c(1,-2,1,1,-5)), built_expected_results_double(localscore = c(0.2,3,4), 
                                                                             suboptimalSeg =  data.frame(
                                                                               value = c(0.1,0.2),
                                                                               begin = c(1,3),
                                                                               end = c(1,4)),
                                                                             RecordTime = c(0,2,5)))
  # Test two mountains in the middle
  expect_equal(built_call_double(0.1*c(-2,1,1,-5,1,-4)), built_expected_results_double(localscore = c(0.2,2,3), 
                                                                                suboptimalSeg =  data.frame(
                                                                                  value = c(0.2,0.1),
                                                                                  begin = c(2,5),
                                                                                  end = c(3,5)),
                                                                     RecordTime = c(0, 1, 4, 6)))
})


test_that("Errors and warnings checks", {
  # Test warning when no local score found in input sequence
  expect_warning(localScoreC(c(-1,-1,-1,-1,-1), suppressWarnings=FALSE),"No local score found")
  # Test warning when postive mean of input sequence
  expect_warning(localScoreC(c(-1,10,-1,-1,-1), suppressWarnings=FALSE),"The mean of this sequence is greater than 0. The sequence may be trivial")
  # Test warning when no local score found in input sequence
  expect_warning(localScoreC_double(c(-0.1,0.-1,0.-1,0.-1,0.-1), suppressWarnings=FALSE),"No local score found")
  # Test warning when postive mean of input sequence
  expect_warning(localScoreC_double(c(-0.1,1.0,-0.1,-0.1,-0.1), suppressWarnings=FALSE),"The mean of this sequence is greater than 0. The sequence may be trivial")
  })


test_that("Type output", {
  expect_false(is.integer(localScoreC(c(-2,1,1,-5,1,-4))$localScore))  #attention à ce comportement
  expect_true(is.integer(localScoreC(as.integer(c(-2,1,1,-5,1,-4)))$localScore))  #attention à ce comportement
  expect_false(is.integer(localScoreC(0.5 + c(-2,1,1,-5,1,-4))$localScore))
  expect_true(is.numeric(localScoreC(0.5 + c(-2,1,1,-5,1,-4))$localScore))
})