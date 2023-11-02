context("Score local detection")
library(localScore)

test_that("Results checks", {
  #Helper function to build object return by localScoreC
  built_expected_results = function(localscore, suboptimalSeg, RecordTime) {
    names(localscore) <- c("value", "begin", "end")
    colnames(suboptimalSeg) <- c("value", "begin", "end")
    list(localScore=localscore, suboptimalSegmentScores=suboptimalSeg, RecordTime=RecordTime)
  }
  # Test with all null values as input
  expect_equal(localScoreC(c(0,0,0,0,0), supressWarnings=TRUE), built_expected_results(localscore = c(0,0,0), 
                                                                                       suboptimalSeg =  matrix(c(0,0,0), ncol=3, byrow=TRUE),
                                                                                       RecordTime = integer(0)))
  # Test with all negative values as input
  expect_equal(localScoreC(c(-1,-1,-1,-1,-1), supressWarnings=TRUE), built_expected_results(localscore = c(0,0,0), 
                                                                                            suboptimalSeg =  matrix(c(0,0,0), ncol=3, byrow=TRUE),
                                                                                            RecordTime = integer(0)))
  # Test with all positive values as input
  expect_equal(localScoreC(c(1,2,3,4,5), supressWarnings=TRUE), built_expected_results(localscore = c(15,1,5), 
                                                                                       suboptimalSeg =  matrix(c(15,1,5), ncol=3, byrow=TRUE),
                                                                                       RecordTime = 1))
  # Test local score in the middle and sequence beginning with a mountain
  expect_equal(localScoreC(c(1,-2,1,3,-5), supressWarnings=TRUE), built_expected_results(localscore = c(4,3,4), 
                                                                                         suboptimalSeg =  matrix(c(1,1,1,
                                                                                                                   4,3,4), ncol=3, byrow=TRUE),
                                                                                         RecordTime = c(1,3)))
  # Test local score in the middle and sequence ending with a mountain
  expect_equal(localScoreC(c(-2,1,1,-5,1), supressWarnings=TRUE), built_expected_results(localscore = c(2,2,3), 
                                                                                         suboptimalSeg =  matrix(c(2,2,3,
                                                                                                                   1,5,5), ncol=3, byrow=TRUE),
                                                                                         RecordTime = c(2,5)))
  #localScoreC(c(3,-2,1,1,-5), supressWarnings=TRUE) (problème : soit une seule montagne de 2 à 5 soit 2 montagnes [2-2] [4-5] ce qui n'est pas le cas)
  # Pour le dire autrement (?) : la fonction localScoreC ne renvoie pas le "même" résultat sur cette séquence, et sa séquence prise dans l'autre sens.
  expect_equal(localScoreC(c(-1,3,-2,1,1,-5), supressWarnings=TRUE), built_expected_results(localscore = c(3,2,2), 
                                                                                            suboptimalSeg =  matrix(c(3,2,2), ncol=3, byrow=TRUE),
                                                                                           RecordTime = 2))
  # Same problem that above, but different behaviour
  expect_equal(localScoreC(c(1,-2,1,1,-5), supressWarnings=TRUE), built_expected_results(localscore = c(2,3,4), 
                                                                                         suboptimalSeg =  matrix(c(1,1,1,
                                                                                                                   2,3,4), ncol=3, byrow=TRUE),
                                                                                         RecordTime = c(1,3)))
  # Test two mountains in the middle
  expect_equal(localScoreC(c(-2,1,1,-5,1,-4), supressWarnings=TRUE), built_expected_results(localscore = c(2,2,3), 
                                                                                            suboptimalSeg =  matrix(c(2,2,3,
                                                                                                                      1,5,5), ncol=3, byrow=TRUE),
                                                                                            RecordTime = c(2,5)))
})

test_that("Errors and warnings checks", {
  # Test warning when no local score found in input sequence
  expect_warning(localScoreC(c(-1,-1,-1,-1,-1), supressWarnings=FALSE),"No local score found")
  # Test warning when postive mean of input sequence
  expect_warning(localScoreC(c(-1,10,-1,-1,-1), supressWarnings=FALSE),"The mean of this sequence is greater than 0. The sequence may be trivial")
})
