#include <iostream>
#include <vector>
#include <math.h>       //round

#include <Rcpp.h>

using namespace Rcpp;

//' @description Calculates the local score for a sequence of integer scores. Only provides the
//' first occurrence of the local score. Use function suboptimalSegment() or Lindley() to obtain the others localizations of the different realizations of the local score.
//' @title Local score
//' @return A structure containing: the local score value and the begin and end index of the segment realizing this optimal score ; all the local maxima of the Lindley process (non negative excursion) and their begin and ens index ; the record times of the Lindley process but only the ones corresponding to the begin index of non negative excursions
//' @param v : a sequence of integer values as vector.
//' @param supressWarnings : if warnings should not be displayed
//' @examples 
//' seq.OneSegment=c(1,-2,3,1,-1,2)
//' # one segment realizing the local score value
//' localScoreC(seq.OneSegment) 
//' seq.TwoSegments=c(1,-2,3,1,2,-2,-2,-1,1,-2,3,1,2,-1,-2,-2,-1,1)
//' # two segments realizing the local score value
//' localScoreC(seq.TwoSegments) 
//' # only the first realization
//' localScoreC(seq.TwoSegments)$localScore 
//' # all the realization of the local together with the suboptimal ones
//' localScoreC(seq.TwoSegments)$suboptimalSegmentScores 
//' # for small sequences, you can also use lindley() fonction to check if 
//' # several segments achieve the local Score
//' lindley(seq.TwoSegments) 
//' plot(1:length(seq.TwoSegments),lindley(seq.TwoSegments),type='b')
//' seq.TwoSegments.InSameExcursion=c(1,-2,3,2,-1,0,1,-2,-2,-4,1)
//' localScoreC(seq.TwoSegments.InSameExcursion)
//' # lindley() shows two realizations in the same excursion (no 0 value between the two LS values)
//' lindley(seq.TwoSegments.InSameExcursion) 
//' # same beginning index but two possible ending indexes
//' # only one excursion realizes the local score even in there is two possible length of segment
//' localScoreC(seq.TwoSegments.InSameExcursion)$suboptimalSegmentScores 
//' plot(1:length(seq.TwoSegments.InSameExcursion),lindley(seq.TwoSegments.InSameExcursion),type='b')
//' @export
// [[Rcpp::export]] 
List localScoreC(std::vector<int> v, bool supressWarnings = false) {
  NumericMatrix localScores(std::round(v.size()/2), 3);
  int numberOfLocalScores = 0;
  std::vector<int> globalMaxScore;
  globalMaxScore.push_back(0);
  globalMaxScore.push_back(0);
  globalMaxScore.push_back(0);
  std::vector<int> localMaxScore;
  localMaxScore.push_back(0);
  localMaxScore.push_back(0); 
  std::vector<int> tempsRecords;
  int sum = 0;
  float average = accumulate(v.begin(), v.end(), 0.0) / v.size();
  if (average >= 0.0 && !supressWarnings)
    warning("The mean of this sequence is greater than 0. The sequence may be trivial");

  
  for (unsigned int index = 0; index<v.size(); index++){  // 0 est le premier composant de séquence
    if(sum + v[index] <= 0){
      sum = 0;
    } else { // alors on va être dans une "vraie" excursion (ie strictt positive)
      if (sum==0){  // on sera au tout début d'une vraie excursion et on vient de finir la montagen précédente
        if(localMaxScore[0] > 0){ //new temps de record au sens indice de début d'une "vraie" excursion
          localScores(numberOfLocalScores, 0) = localMaxScore[0]; 
          localScores(numberOfLocalScores, 1) = tempsRecords[tempsRecords.size()-1]; 
          localScores(numberOfLocalScores, 2) = localMaxScore[1]; 
          numberOfLocalScores++;
          if(localMaxScore[0]>globalMaxScore[0]){ // mise à jour du score local global (Hn)
            globalMaxScore[0] = localMaxScore[0];
            globalMaxScore[1] = tempsRecords[tempsRecords.size()-1];
            globalMaxScore[2] = localMaxScore[1];
          }
          localMaxScore[0] = 0;
          localMaxScore[1] = 0;
        } // fin if new temps de record
        tempsRecords.push_back(index + 1);//R is not zero indexed
      }// fin if sum==0 on vient de finir une vraie excursion et on va en commencer une autre
      sum += v[index]; // on actualise la valeur du processus de Lindley
      if(sum>localMaxScore[0]){ // une nouvelle valeur de score de la montagne en cours score vient d'être réalisé
        localMaxScore[0] = sum;
        localMaxScore[1] = index + 1 ;//R is not zero indexed on a indexé la fin de la réalization du score en cours
      }
    }// fin else vraie montagne
  } // fin de la lecture de la séquence
  // on étudie la dernière montagne non complète (fin de la séquence)
  if(localMaxScore[0] > 0){
    //adding last local score if there's any
    localScores(numberOfLocalScores, 0) = localMaxScore[0]; 
    localScores(numberOfLocalScores, 1) = tempsRecords[tempsRecords.size()-1]; 
    localScores(numberOfLocalScores, 2) = localMaxScore[1];     
    numberOfLocalScores++;
    if(localMaxScore[0]>globalMaxScore[0]){
      globalMaxScore[0] = localMaxScore[0];
      globalMaxScore[1] = tempsRecords[tempsRecords.size()-1];
      globalMaxScore[2] = localMaxScore[1];
    }
  }// fin if localMaxScore[0] > 0)

  NumericMatrix finalScores(1,3);
  NumericVector glMxScore = NumericVector::create(Named("value") = globalMaxScore[0], Named("begin") = globalMaxScore[1], Named("end") = globalMaxScore[2]);
  if(numberOfLocalScores==0){ // pas de vraie excrusion 
    if(!supressWarnings)
       warning("No local score found");
    finalScores(0,0) = 0;
    finalScores(0,1) = 0;
    finalScores(0,2) = 0;
/*    colnames(finalScores) = CharacterVector::create("value", "[", "]");
    return Rcpp::List::create(
      Named("localScore") = globalMaxScore,
      Named("suboptimalSegmentScores") = finalScores,
      Named("RecordTime") = tempsRecords
    );*/
  } // fin pas de vraie excrusion 
  else { // il y a au moins une vraie excursion 
      finalScores = localScores( Range(0,numberOfLocalScores-1), Range(0,2));
      /*NumericVector glMxScore = NumericVector::create(Named("value") = globalMaxScore[0], Named("begin") = globalMaxScore[1], Named("end") = globalMaxScore[2]);
      return Rcpp::List::create(
        Named("localScore") = glMxScore,
        Named("suboptimalSegmentScores") = finalScores,
        Named("RecordTime") = tempsRecords
      );*/
  } // fin il y a au moins 1 vraie excursion
  colnames(finalScores) = CharacterVector::create("value", "begin", "end");
  return Rcpp::List::create(
    Named("localScore") = glMxScore,
    Named("suboptimalSegmentScores") = finalScores,
    Named("RecordTime") = tempsRecords
  );
} // boucle sur la lecture de la sequence

//' @description Calculates the local score for a sequence of doubles. Only provides the
//' first occurrence. Use function suboptimalSegment() or Lindley() to obtain the others localizations of the different realizations of the local score.
//' @title Local score for sequences of floating values
//' @return A structure containing: the local score value and the begin and end index of the segment realizing this optimal score ; all the local maxima of the Lindley process (non negative excursion) and their begin and ens index ; the record times of the Lindley process but only the ones corresponding to the begin index of non negative excursions 
//' @param v A sequence of values as vector.
//' @param supressWarnings if warnings should be displayed
//' @examples 
//' localScoreC_double(c(1.2,-2.1,3.5,1.7,-1.1,2.3))
//' seq.TwoSegments=c(1.2,-2.1,3.5,1.7,2,-2,-2,-3.5,1,3.5,1.7,1,-2,-2)
//' # two segments realizing the local score value
//' localScoreC(seq.TwoSegments) 
//' # only the first realization
//' localScoreC(seq.TwoSegments)$localScore 
//' # all the realization of the local together with the suboptimal ones
//' localScoreC(seq.TwoSegments)$suboptimalSegmentScores 
//' # for small sequences, you can also use lindley() fonction to check if 
//' # several segments achieve the local score
//' lindley(seq.TwoSegments) 
//' plot(1:length(seq.TwoSegments),lindley(seq.TwoSegments),type='b')
//' seq.TwoSegments.InSameExcursion=c(1,-2,3,2,-1,0,1,-2,-2)
//' localScoreC(seq.TwoSegments.InSameExcursion)
//' # lindley() shows two realizations in the same excursion (no 0 value between the two LS values)
//' lindley(seq.TwoSegments.InSameExcursion) 
//' plot(1:length(seq.TwoSegments.InSameExcursion),lindley(seq.TwoSegments.InSameExcursion),type='b')
//' # same beginning index but two possible ending indexes
//' # only one excursion realizes the local score even in there is two possible length of segment
//' localScoreC(seq.TwoSegments.InSameExcursion)$suboptimalSegmentScores 
//' @export
// [[Rcpp::export]] 
List localScoreC_double(std::vector<double> v, bool supressWarnings = false) {
  NumericMatrix localScores(round(v.size()/2), 3);
  int numberOfLocalScores = 0;
  std::vector<double> globalMaxScore;
  globalMaxScore.push_back(0);
  globalMaxScore.push_back(0);
  globalMaxScore.push_back(0);
  std::vector<double> localMaxScore;
  localMaxScore.push_back(0);
  localMaxScore.push_back(0); 
  std::vector<int> tempsRecords;
  double sum = 0;
  float average = accumulate(v.begin(), v.end(), 0.0) / v.size();
  if (average >= 0.0 && !supressWarnings)
    warning("The mean of this sequence is greater than 0. The sequence may be trivial");
  
  
  for (unsigned int index = 0; index<v.size(); index++){
    if(sum + v[index] <= 0){
      sum = 0.0;
    } else {
      if (sum==0.0){  //new temps de record
        if(localMaxScore[0] > 0.0){
          localScores(numberOfLocalScores, 0) = localMaxScore[0]; 
          localScores(numberOfLocalScores, 1) = tempsRecords[tempsRecords.size()-1]; 
          localScores(numberOfLocalScores, 2) = localMaxScore[1]; 
          numberOfLocalScores++;
          if(localMaxScore[0]>globalMaxScore[0]){
            globalMaxScore[0] = localMaxScore[0];
            globalMaxScore[1] = tempsRecords[tempsRecords.size()-1];
            globalMaxScore[2] = localMaxScore[1];
          }
          localMaxScore[0] = 0.0;
          localMaxScore[1] = 0.0;
        }
        tempsRecords.push_back(index + 1);//R is not zero indexed
      }
      sum += v[index];
      if(sum>localMaxScore[0]){
        localMaxScore[0] = sum;
        localMaxScore[1] = index + 1 ;//R is not zero indexed
      }
    }
  }
  if(localMaxScore[0] > 0.0){
    //adding last local score if there's any
    localScores(numberOfLocalScores, 0) = localMaxScore[0]; 
    localScores(numberOfLocalScores, 1) = tempsRecords[tempsRecords.size()-1]; 
    localScores(numberOfLocalScores, 2) = localMaxScore[1];     
    numberOfLocalScores++;
    if(localMaxScore[0]>globalMaxScore[0]){
      globalMaxScore[0] = localMaxScore[0];
      globalMaxScore[1] = tempsRecords[tempsRecords.size()-1];
      globalMaxScore[2] = localMaxScore[1];
    }
  }
  
  NumericMatrix finalScores;
  if(numberOfLocalScores==0){
    if(!supressWarnings)
      warning("No local score found");
    return Rcpp::List::create(
      Named("maxScore") = globalMaxScore,
      Named("stopping times") = tempsRecords
    );
  } else {
    finalScores = localScores( Range(0,numberOfLocalScores-1), Range(0,2));
    colnames(finalScores) = CharacterVector::create("value", "begin", "end");
    NumericVector glMxScore = NumericVector::create(Named("value") = globalMaxScore[0], Named("begin") = globalMaxScore[1], Named("end") = globalMaxScore[2]);
    return Rcpp::List::create(
      Named("localScore") = glMxScore,
      Named("suboptimalSegmentScores") = finalScores,
      Named("RecordTime") = tempsRecords
    );
  }
}
