#include <iostream>
#include <vector>
#include <math.h>       //round
#include <cmath>        //fpclassify() and FP_ZERO

#include <Rcpp.h>

using namespace Rcpp;

//template <typename T>
template <int RTYPE>
List localScoreC_T(const Vector<RTYPE>& v, bool suppressWarnings = false) {
  Vector<RTYPE> localScoresValues ;  // Les Excursions (valeur)
  IntegerVector localScoresBegin ;   // Les Excursions (début)
  IntegerVector localScoresEnd ;     // Les Excursions (fin)
  int numberOfLocalScores = 0;       // Number of Excursions
  int localScoreStart = 0 ;          // Excursion en cours (début)
  Vector<RTYPE> globalMaxScoreValue(1) ; // Score Local (valeur)
  // Rcpp::traits::storage_type<RTYPE>::type globalMaxScoreValue ; // A tester
  int globalMaxScoreBegin ;          // Score Local (Debut)
  int globalMaxScoreEnd ;            // Score Local (Fin)
  globalMaxScoreValue[0] = 0;
  globalMaxScoreBegin = 0;
  globalMaxScoreEnd = 0;
  Vector<RTYPE> localMaxScoreValue(1) ;    // Excursion en cours (valeur)
  int localMaxScoreEnd ;                   // Excursion en cours (fin)
  localMaxScoreValue[0] = 0;
  localMaxScoreEnd = 0; 
  
  std::vector<int> tempsRecords;     // Temps de Record au sens Karlin et Dembo (1990)
  Vector<RTYPE> sum(1) ;
  Vector<RTYPE> next_sum(1);
  sum[0] = 0 ;
  next_sum[0] = 0 ;

  float average = 0.0 ;
  for (unsigned int i = 0 ; i < v.size() ; i++)
    average += v[i] ;
  if (average >= 0.0 && !suppressWarnings)
     warning("The mean of this sequence is greater than 0. The sequence may be trivial");
  
  tempsRecords.push_back(0); // By convention, first record time is zero
  for (unsigned int index = 0; index<v.size(); index++){  // 0 est le premier composant de séquence
    next_sum[0] = sum[0]+v[index];
    if(next_sum[0] < 0) { // Nouveau temps de record au sens Karlin et Dembo
      tempsRecords.push_back(index + 1);//R is not zero indexed
      sum[0] = 0;
    } else if (fpclassify(next_sum[0]) == FP_ZERO){  // compare to zero (also for double)
      sum[0] = 0;
    } else { // alors on va être dans une "vraie" excursion (ie strictly positive)
      if (fpclassify(sum[0]) == FP_ZERO){  // on sera au tout début d'une vraie excursion et on vient de finir la montagne précédente
        if(localMaxScoreValue[0] > 0){ //new temps de record au sens indice de début d'une "vraie" excursion
          localScoresValues.push_back(localMaxScoreValue[0]); 
          localScoresBegin.push_back(localScoreStart); 
          //          localScores(numberOfLocalScores, 1) = tempsRecords[tempsRecords.size()-1]; 
          localScoresEnd.push_back(localMaxScoreEnd); 
          numberOfLocalScores++;
          if(localMaxScoreValue[0]>globalMaxScoreValue[0]){ // mise à jour du score local global (Hn)
            globalMaxScoreValue[0] = localMaxScoreValue[0];
            globalMaxScoreBegin = localScoreStart;
            //            globalMaxScore[1] = tempsRecords[tempsRecords.size()-1];
            globalMaxScoreEnd = localMaxScoreEnd;
          }
          localMaxScoreValue[0] = 0;
          localMaxScoreEnd = 0;
        } // fin if new temps de record
        //tempsRecords.push_back(index + 1);//R is not zero indexed
        localScoreStart = index + 1; 
      }// fin if sum==0 on vient de finir une vraie excursion et on va en commencer une autre
      sum[0] += v[index]; // on actualise la valeur du processus de Lindley
      if(sum[0]>localMaxScoreValue[0]){ // une nouvelle valeur de score de la montagne en cours score vient d'être réalisé
        localMaxScoreValue[0] = sum[0];
        localMaxScoreEnd = index + 1 ;//R is not zero indexed on a indexé la fin de la réalization du score en cours
      }
    }// fin else vraie montagne
  } // fin de la lecture de la séquence
  // on étudie la dernière montagne non complète (fin de la séquence)
  if(localMaxScoreValue[0] > 0){
    //adding last local score if there's any
    localScoresValues.push_back(localMaxScoreValue[0]); 
    localScoresBegin.push_back(localScoreStart); 
    //    localScores(numberOfLocalScores, 1) = tempsRecords[tempsRecords.size()-1]; 
    localScoresEnd.push_back(localMaxScoreEnd);     
    numberOfLocalScores++;
    if(localMaxScoreValue[0]>globalMaxScoreValue[0]){
      globalMaxScoreValue[0] = localMaxScoreValue[0];
      globalMaxScoreBegin = localScoreStart;
      //      globalMaxScore[1] = tempsRecords[tempsRecords.size()-1];
      globalMaxScoreEnd = localMaxScoreEnd;
    }
  }// fin if localMaxScore[0] > 0)
  
  int fixNoScoreMatrixSize = (numberOfLocalScores==0) ? 1 : 0; // needed to create finalScores Matrix in case of no localScore found
//  NumericMatrix finalScores(numberOfLocalScores+fixNoScoreMatrixSize,3);
  Vector<RTYPE> glMxScore = Vector<RTYPE>::create(Named("value") = globalMaxScoreValue[0],
                                                  Named("begin") = globalMaxScoreBegin,
                                                  Named("end") = globalMaxScoreEnd);

  if(numberOfLocalScores==0){ // pas de vraie excrusion 
    if(!suppressWarnings)
      warning("No local score found");
    localScoresValues.push_back(0); 
    localScoresBegin.push_back(0); 
    localScoresEnd.push_back(0);     
  } // fin pas de vraie excrusion 
  // else { // il y a au moins une vraie excursion 
  //   finalScores = localScores( Range(0,numberOfLocalScores-1), Range(0,2));
  // } // fin il y a au moins 1 vraie excursion
  //colnames(finalScores) = CharacterVector::create("value", "begin", "end");
  DataFrame finalScores = DataFrame::create(Named("value") = localScoresValues,
                                            Named("begin") = localScoresBegin,
                                            Named("end") = localScoresEnd) ;
  return Rcpp::List::create(
    Named("localScore") = glMxScore,
    Named("suboptimalSegmentScores") = finalScores,
    Named("RecordTime") = tempsRecords
  );
} 

//' @title Local score
//' @description Calculates the local score for a sequence of scores, the sub-optimal segments, and the associated record times. The local score is the maximal sum of values contained in a segment among all possible segments of the sequence. In other word, it generalizes a sliding window approach, considering of all possible windows size.
//' @param v a sequence of numerical values as vector (integer or double).
//' @param suppressWarnings (optional) if warnings should not be displayed
//' @return A list containing:
//' \item{localScore}{the local score value and the begin and end index of the segment realizing this optimal score;}
//' \item{suboptimalSegmentScores}{An array containing sub-optimal local scores, that is all the local maxima of the Lindley rocess (non negative excursion) and their begin and end index;}
//' \item{RecordTime}{The record times of the Lindley process as defined in Karlin and Dembo (1990).}
//' @details The \code{localScoreC} function is implemented in a templated C function. Be aware that the type of the output (\code{integer} or \code{double}) depends on the type of the input. The function \code{localScoreC_double} \code{localScoreC_int} explicitly use the corresponding type (with an eventual conversion in case of integer). Warning: in R, \code{typeof(c(1,3,4,10)) == "double"}. You can set a type of a vector with \code{mode()} or \code{as.integer()} functions for example. \cr
//' \code{localScoreC_int} is just a call to \code{as.integer()} before calling \code{localScoreC}. \code{localScore_double} is just a call to \code{localScoreC}, and as such is deprecated.
//' @seealso \code{\link{lindley}}
//' @examples 
//' localScoreC(c(1.2,-2.1,3.5,1.7,-1.1,2.3))
//' # one segment realizing the local score value
//' seq.OneSegment <- c(1,-2,3,1,-1,2)
//' localScoreC(seq.OneSegment) 
//' seq.TwoSegments <- c(1,-2,3,1,2,-2,-2,-1,1,-2,3,1,2,-1,-2,-2,-1,1)
//' # two segments realizing the local score value
//' localScoreC(seq.TwoSegments) 
//' # only the first realization
//' localScoreC(seq.TwoSegments)$localScore 
//' # all the realization of the local together with the suboptimal ones
//' localScoreC(seq.TwoSegments)$suboptimalSegmentScores 
//' # for small sequences, you can also use lindley() function to check if 
//' # several segments achieve the local Score
//' lindley(seq.TwoSegments) 
//' plot(1:length(seq.TwoSegments),lindley(seq.TwoSegments),type='b')
//' seq.TwoSegments.InSameExcursion <- c(1,-2,3,2,-1,0,1,-2,-2,-4,1)
//' localScoreC(seq.TwoSegments.InSameExcursion)
//' # lindley() shows two realizations in the same excursion (no 0 value between the two LS values)
//' lindley(seq.TwoSegments.InSameExcursion) 
//' # same beginning index but two possible ending indexes
//' # only one excursion realizes the local score even in there is two possible lengths of segment
//' localScoreC(seq.TwoSegments.InSameExcursion)$suboptimalSegmentScores 
//' plot(1:length(seq.TwoSegments.InSameExcursion),lindley(seq.TwoSegments.InSameExcursion),type='b')
//' # Technical note about type correspondance
//' seq.OneSegment <- c(1,-2,3,1,-1,2)
//' seq.OneSegmentI <- as.integer(seq.OneSegment)
//' typeof(seq.OneSegment)  # "double" (beware)
//' typeof(seq.OneSegmentI) # "integer"
//' LS1 <- localScoreC(seq.OneSegment, suppressWarnings = TRUE)
//' LS1I <- localScoreC(seq.OneSegmentI, suppressWarnings = TRUE)
//' typeof(LS1$localScore)  # "double"
//' typeof(LS1I$localScore) # "integer"
//' typeof(LS1$suboptimalSegmentScores$value)  # "double"
//' typeof(LS1I$suboptimalSegmentScores$value) # "integer"
//' # Force to use integer values (trunk values if necessary)
//' seq2 <- seq.OneSegment + 0.5
//' localScoreC(seq2, suppressWarnings = TRUE)
//' localScoreC_int(seq.OneSegment, suppressWarnings = TRUE)
//' @export
// [[Rcpp::export]] 
SEXP localScoreC(SEXP v, bool suppressWarnings = false) {
  // see for inspiration : https://gallery.rcpp.org/articles/rcpp-return-macros/
  switch (TYPEOF(v)) {
  case INTSXP: {
    return localScoreC_T(as<IntegerVector>(v),suppressWarnings) ;
  }
  case REALSXP: {
    return localScoreC_T(as<NumericVector>(v),suppressWarnings) ;
  }
  default: {
    warning(
      "localScoreC() : Invalid input SEXPTYPE %d (%s).\n",
      TYPEOF(v), type2name(v)
    );
    return R_NilValue;
  }
  }
}
