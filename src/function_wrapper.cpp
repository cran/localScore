#include <Rcpp.h>
#include <cstdlib>
#include <string>
#include <iostream>


#include "pValueMethods.h"
//#include "malloc.h"

#include "Eigen/Dense"

using namespace Rcpp;
using Eigen::MatrixXd;

//' @description Calculates the exact p-value in the identically and independantly distributed of a given local score, a sequence length that 'must not be too large' and for a given score distribution
//' @details Small in this context depends heavily on your machine. On a 3,7GHZ machine this means for daudin(1000, 5000, c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1), -2, 3)
//' an execution time of ~2 seconds. This is due to the calculation method using matrix exponentation which becomes very fast very slow. The size of the matrix of the exponentiation is equal to a+1 with a the local score value. The matrix must be put at the power n, with n the sequence length.
//' Moreover, it is known that the local score value is expected to be in mean of order log(n).
//' @title Daudin [p-value] [iid]
//' @return A double representing the probability of a local score as high as the one given as argument
//' @param localScore the observed local score
//' @param sequence_length length of the sequence
//' @param score_probabilities the probabilities for each score from lowest to greatest
//' @param sequence_min minimum score
//' @param sequence_max maximum score
//' @examples 
//' daudin(localScore = 4, sequence_length = 50, 
//' score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), sequence_min = -3, sequence_max = 2)
//' @export
// [[Rcpp::export]]
double daudin(int localScore, int sequence_length, NumericVector score_probabilities, int sequence_min, int sequence_max){
  if(localScore<0)
    stop("[Invalid Input] local score must be positive.");
  if(sequence_length<= 0 )
    stop("[Invalid Input] sequence length must be positive.");
  if(score_probabilities.size()!= sequence_max - sequence_min + 1)
    stop("[Invalid Input] score probability distribution must contain as much elements as the range from sequence_min to sequence_max.");
  if(sequence_max<=0)
    stop("[Invalid Input] sequence_max must be positive.");
  if(sequence_min>=0)
    stop("[Invalid Input] sequence_min must be negative.");
  return calcul_daudin(localScore, sequence_length, as<std::vector<double>>(score_probabilities), sequence_min, sequence_max);
}

//' @description Calculates an approximated p-value of a given local score value and a long sequence length in the identically and independantly distributed model for the sequence. See also mcc() function for another approximated method in the i.i.d. model 
//' @details This method works the better the longer the sequence is. Important note : the calculus of the parameter of the distribution uses
//' the resolution of a polynome which is a function of the score distribution, of order max(score)-min(score). There exists only empirical methods to solve a polynome of order greater that 5
//' with no warranty of reliable solution.
//' The found roots are checked internally to the function and an error message is throw in case of inconsistent. In such case, you could try to change your score scheme (in case of discretization)
//' or use the function \code{\link{karlinMonteCarlo}} .
//' @title Karlin [p-value] [iid]
//' @return A double representing the probability of a localScore as high as the one given as argument
//' @param localScore the observed local score
//' @param sequence_length length of the sequence (at least several hundreds)
//' @param score_probabilities the probabilities for each unique score from lowest to greatest
//' @param sequence_min minimum score
//' @param sequence_max maximum score
//' @examples 
//' karlin(150, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -5, 5)
//' @export
// [[Rcpp::export]]
double karlin(int localScore, int sequence_length, NumericVector score_probabilities, int sequence_min, int sequence_max){
  if(localScore<0)
    stop("[Invalid Input] local score must be positive.");
  if(sequence_length<= 0 )
    stop("[Invalid Input] sequence length must be positive.");
  if(score_probabilities.size()!= sequence_max - sequence_min + 1)
    stop("[Invalid Input] score probability distribution must contain as much elements as the range from sequence_min to sequence_max.");
  if(sequence_max<=0)
    stop("[Invalid Input] sequence_max must be positive.");
  if(sequence_min>=0)
    stop("[Invalid Input] sequence_min must be negative.");
  double esp = 0.0 ;
  for (int i=sequence_min ; i<=sequence_max ; i++) {
    esp += i*score_probabilities[i-sequence_min] ;
  }
  if(esp>=0.0)
    stop("[Invalid Input] Score expectation must be strictly negative.");
  double p = calcul_karlin(localScore, as<std::vector<double>>(score_probabilities), sequence_max, -sequence_min, (long)sequence_length);
  if (fabs(p+1.0)<1e-10) //p==-1.0 in case of error in calcul_karlin (polynomial roots problem)
    stop("karlin() function cannot be used in your case due to numerical instability (polynomial roots solver). Check the documentation of 'karlin()' for details.\n You could try to change your scoring discretisation step or use karlinMonteCarlo()");
  if (fabs(p+2.0)<1e-10) //p==-2.0 in case of error in calcul_karlin (polynomial roots problem)
    stop("ERROR karlin: u and/or v are not compatible with the size of 'distribution'");
  return p;
}

//' @description Calculates an approximated p-value for a given local score value and a medium to long sequence length in the identically and independantly distributed model
//' @details This methods is actually an improved method of Karlin and produces more precise results. It should be privileged whenever possible. \cr
//' As with karlin, the method works the better the longer the sequence. Important note : the calculus of the parameter of the distribution uses
//' the resolution of a polynome which is a function of the score distribution, of order max(score)-min(score). There exists only empirical methods to solve a polynome of order greater that 5
//' with no warranty of reliable solution.
//' The found roots are checked internally to the function and an error message is throw in case of inconsistency. In such case, you could try to change your score scheme (in case of discretization)
//' or use the function \code{\link{karlinMonteCarlo}} .
//' @title MCC [p-value] [iid]
//' @return A double representing the probability of a local score as high as the one given as argument
//' @param localScore the observed local score
//' @param sequence_length length of the sequence (up to one hundred)
//' @param score_probabilities the probabilities for each unique score from lowest to greatest
//' @param sequence_min minimum score
//' @param sequence_max maximum score
//' @examples 
//' mcc(40, 100, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
//' mcc(40, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
//' @export
// [[Rcpp::export]]
double mcc(int localScore, int sequence_length, NumericVector score_probabilities, int sequence_min, int sequence_max){
  if(localScore<0)
    stop("[Invalid Input] local score must be strictly positive.");
  if(sequence_length<= 0 )
    stop("[Invalid Input] sequence length must be strictly positive.");
  if(score_probabilities.size()!= sequence_max - sequence_min + 1)
    stop("[Invalid Input] score probability distribution must contain as much elements as the range from sequence_min to sequence_max.");
  if(sequence_max<=0)
    stop("[Invalid Input] sequence_max must be strictly positive.");
  if(sequence_min>=0)
    stop("[Invalid Input] sequence_min must be strictly negative.");
  double esp = 0.0 ;
  for (int i=sequence_min ; i<=sequence_max ; i++) {
    esp += i*score_probabilities[i-sequence_min] ;
  }
  if(esp>=0.0)
    stop("[Invalid Input] Score expectation must be strictly negative.");
  double p = calcul_mcc(localScore, as<std::vector<double>>(score_probabilities), sequence_max, -sequence_min, (long)sequence_length);
  if (fabs(p+1.0)<1e-10) //p==-1.0 in case of error in calcul_karlin
    stop("mcc() function cannot be used in your case. Check the documentation of 'mcc()' for details.\n You could try to change your scoring discretisation step or use karlinMonteCarlo()");
  return p;
}

//' @description Calculates the distribution of the maximum of the partial sum process for a given value in the identically and independantly distributed model
//' @details Implement the formula (4) of the article Mercier, S., Cellier, D., & Charlot, D. (2003). An improved approximation for assessing the statistical significance of molecular sequence features. Journal of Applied Probability, 40(2), 427-441. doi:10.1239/jap/1053003554 \cr
//' Important note : the calculus of the parameter of the distribution uses
//' the resolution of a polynome which is a function of the score distribution, of order max(score)-min(score). There exists only empirical methods to solve a polynome of order greater that 5
//' with no warranty of reliable solution.
//' The found roots are checked internally to the function and an error message is throw in case of inconsistency. 
//' @title Maximum of the partial sum [probability] [iid]
//' @return A double representing the probability of the maximum of the partial sum process equal to k
//' @param k value at which calculates the probability
//' @param score_probabilities the probabilities for each unique score from lowest to greatest
//' @param sequence_min minimum score
//' @param sequence_max maximum score
//' @examples 
//' maxPartialSumd(10, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
//' @export
// [[Rcpp::export]]
double maxPartialSumd(int k, NumericVector score_probabilities, int sequence_min, int sequence_max){
     if(k<0)
     stop("[Invalid Input] local score must be strictly positive.");
   if(score_probabilities.size()!= sequence_max - sequence_min + 1)
     stop("[Invalid Input] score probability distribution must contain as much elements as the range from sequence_min to sequence_max.");
   if(sequence_max<=0)
     stop("[Invalid Input] sequence_max must be strictly positive.");
   if(sequence_min>=0)
     stop("[Invalid Input] sequence_min must be strictly negative.");
   double esp = 0.0 ;
   for (int i=sequence_min ; i<=sequence_max ; i++) {
     esp += i*score_probabilities[i-sequence_min] ;
   }
   if(esp>=0.0)
     stop("[Invalid Input] Score expectation must be strictly negative.");
   double p = calcul_probMaxPartialSum(k, as<std::vector<double>>(score_probabilities), sequence_max, -sequence_min);
   if (fabs(p+1.0)<1e-10) //p==-1.0 in case of error in calcul_probMaxPartialSum
     stop("probMaxPartialSum() function cannot be used in your case. Check the documentation of 'probMaxPartialSum()' for details.\n You could try to change your scoring discretisation step or use karlinMonteCarlo()");
   return p;
 }


//' @description Calculates stationary distribution of markov transition matrix by use of eigenvectors of length 1
//' @title Stationary distribution [Markov chains]
//' @return A vector with the probabilities
//' @param m Transition Matrix [matrix object]
//' @examples 
//' B = t(matrix (c(0.2, 0.8, 0.4, 0.6), nrow = 2))
//' stationary_distribution(B)
//' @export
// [[Rcpp::export]]
NumericVector stationary_distribution(NumericMatrix m){
  // std::accumulate(x.begin(),x.end(), 0.0)
  float rowsum ;
  int i, j ;
  for(i = 0; i< m.nrow(); i++){
    NumericMatrix::Row row_i = m(i,_) ;
    rowsum = 0.0 ;
    for ( j = 0; j<row_i.size();j++) {
      rowsum += row_i[j] ;
    }
    if (rowsum != 1.0 ) {
      String the_message = "[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.). Sum of line " + std::to_string(i+1) + " equal " + std::to_string(rowsum) ;
      stop(the_message) ;
    }
  }

  MatrixXd mm(m.nrow(), m.ncol());
  for(int i = 0; i< m.nrow(); i++){
    for(int j = 0; j< m.ncol(); j++)
      mm(i, j) = m(i, j);
  }
  
  std::vector<Eigen::VectorXcd> v = stationary_distribution(mm); // XXX revoir ce bordel avec les types pour le simplifier
  NumericVector results(v[0].size());
  if(v.size()> 1){
    stop("Markov matrix is not irreductible (many eigenvalues == 1).");
  }
  if(v.size() == 0){
    stop("no eigenvector found.");
  } else {
    for(int i = 0; i < v[0].size(); i++)
      results[i] = v[0](i).real();
    return results;
  }
}

// XXX Pour les problèmes de warning CRAN concernant la doc, voir l'astuce dans la deuxième réponse : https://stackoverflow.com/questions/26263479/function-param-default-value-stdvector-initialization-with-rcpp-and-c11
// Voir aussi : https://github.com/RcppCore/Rcpp/issues/190 et https://github.com/Rcpp11/attributes/issues/36
//' @description Calculates the exact p-value for short numerical Markov chains. Memory usage and time computation can be too large for a high local score value and high score range (see details).
//' @title Exact method for p-value [Markov chains]
//' @return A double representing the probability of a localScore as high as the one given as argument
//' @param localScore Integer local score for which the p-value should be calculated
//' @param m Transition matrix [matrix object]. Optionnaly, rownames can be corresponding score values. m should be a transition matrix of an ergodic Markov chain.
//' @param sequence_length Length of the sequence
//' @param score_values A integer vector of sequence score values (optional). If not set, the rownames of m are used if they are numeric and set.
//' @param prob0 Vector of probability distribution of the first score of the sequence (optional). If not set, the stationnary distribution of m is used.
//' @details This method computation needs to allocate a square matrix of size localScore^(range(score_values)). This matrix is then exponentiated to sequence_length.
//' @examples 
//' mTransition <- t(matrix(c(0.2, 0.3, 0.5, 0.3, 0.4, 0.3, 0.2, 0.4, 0.4), nrow = 3))
//' scoreValues <- -1:1
//' initialProb <- stationary_distribution(mTransition)
//' exact_mc(localScore = 12, m = mTransition, sequence_length = 100, 
//'         score_values = scoreValues, prob0 = initialProb)
//' exact_mc(localScore = 150, m = mTransition, sequence_length = 1000, 
//'          score_values = scoreValues, prob0 = initialProb)
//' rownames(mTransition) <- scoreValues
//' exact_mc(localScore = 12, m = mTransition, sequence_length = 100, prob0 = initialProb)
//' # Minimal specification
//' exact_mc(localScore = 12, m = mTransition, sequence_length = 100)
//' @export
// [[Rcpp::export]]
double exact_mc(int localScore, NumericMatrix m, int sequence_length, Nullable<NumericVector> score_values = R_NilValue, Nullable<NumericVector> prob0 = R_NilValue){
//  double exact_mc(int localScore, NumericMatrix m, int sequence_length, NumericVector score_values = Rcpp::NumericVector::create(), NumericVector prob0 = Rcpp::NumericVector::create()){
    
   double eps = 1e-12 ;
  NumericVector score_values_ ;
  NumericVector prob0_ ;
  if(score_values.isUsable()) {
    score_values_=score_values;
  } else {
    score_values_ = NumericVector::create();
  }
  if(prob0.isUsable()) {
    prob0_=prob0;
  } else {
    prob0_ = NumericVector::create();;
  }
  
  if(localScore<0) 
     stop("[Invalid Input] local score must be positive.");
   if(sequence_length<=0)
     stop("[Invalid Input] sequence length must be positive.");
   if (m.nrow()!= m.ncol())
     stop("[ERROR exact_mc : Invalid Input] m should be a square matrix");
   if ((score_values_.size() != 0) && (m.nrow()!= score_values_.size() || m.ncol()!= score_values_.size())) 
     stop("[ERROR exact_mc : Invalid Input] m should be a square matrix of size the length of score_values");
   if ((prob0_.size() != 0) && (m.nrow()!= prob0_.size() || m.ncol()!= prob0_.size())) 
     stop("[ERROR exact_mc : Invalid Input] prob0 size should be equal to the number of rows of m");
   if ((prob0_.size() != 0) && (fabs(sum(prob0_) - 1.0)> eps))
     stop("[ERROR exact_mc : Invalid Input] prob0 vector should sum to 1");
   for(NumericVector::iterator i = prob0_.begin(); i != prob0_.end(); ++i) {
     //if ((fabs(*i) > eps) || (fabs(*i-1.0) > eps)) {
       if ((*i < -eps) || (*i > 1.0+eps)) {
         stop("[ERROR exact_mc : Invalid Input] prob0 vector should contains values between 0 and 1");
     }
   }
   
   // Deal with default parameters of R call
   if (score_values_.size() == 0) { // utilisation de la valeur par défaut
     if ((m.attr("dimnames")==R_NilValue) || (rownames(m)==R_NilValue)){
       stop("[ERROR exact_mc : Invalid Input] Either m is a matrix  with score values as rownames or specify score values via score_values parameter");
     } else { 
       CharacterVector ch = rownames(m);
       //     if (Rf_isNull(ch)) {
       if (ch.size() != m.ncol()) {
         stop("[ERROR exact_mc : Invalid Input] Either m is a matrix  with score values as rownames or specify score values via score_values parameter");
       } else { 
         // Initialisation of score values
         NumericVector out(ch.size());
         std::transform(ch.begin(), ch.end(), out.begin(), std::atoi);
         score_values_ = out;
         // XXX TODO gestion erreur ou les rownames ne sont pas numeric
       }
     }
   }
   if(max(score_values_)<=0)
     stop("[Invalid Input] sequence_max must be positive.");
   if(min(score_values_)>=0)
     stop("[Invalid Input] sequence_min must be negative.");
   
   if(prob0_.size() == 0) { // utilisation de la valeur par défaut (donné par stationary_distribution(m))
     NumericVector stationnary = stationary_distribution(m);
     prob0_ = stationnary ;
   }
   
   MatrixXd mm(m.nrow(), m.ncol());
   Eigen::VectorXi scoreValues(m.nrow());
   Eigen::VectorXd probInitial(m.nrow());
   for(int i = 0; i< m.nrow(); i++){
     scoreValues(i) = score_values_(i) ;
     probInitial(i) = prob0_(i) ;
     for(int j = 0; j< m.ncol(); j++)
       mm(i, j) = m(i, j);
   }
   return mh_markov(localScore, mm, scoreValues, sequence_length, probInitial);
 }
