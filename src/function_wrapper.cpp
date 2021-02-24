#include <Rcpp.h>
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
//' @details This method works the better the longer the sequence is. 
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
  return calcul_karlin(localScore, as<std::vector<double>>(score_probabilities), sequence_max, -sequence_min, (long)sequence_length);
  }

//' @description Calculates an approximated p-value for a given local score value and a medium to long sequence length in the identically and independantly distributed model
//' @details This methods is actually an improved method of Karlin and produces more precise results. It should be privileged whenever possible. \cr
//' As with karlin, the method works the better the longer the sequence.
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
  return calcul_mcc(localScore, as<std::vector<double>>(score_probabilities), sequence_max, -sequence_min, (long)sequence_length);
}

//' @description Calculates stationary distribution of markov transitition matrix by use of eigenvectors of length 1
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
      //stop(std::strcat("[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.)",to_string(rowsum))) ;
      //stop(strcat("[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.)", rowsum);
    }
  }

  MatrixXd mm(m.nrow(), m.ncol());
  for(int i = 0; i< m.nrow(); i++){
    for(int j = 0; j< m.ncol(); j++)
      mm(i, j) = m(i, j);
  }
  
  // check that m is a stochastic matrix (row sum == 1 for each row)
  for(int i = 0; i< m.nrow(); i++){
    for(int j = 0; j< m.ncol(); j++)
      mm(i, j) = m(i, j);
  }
  std::vector<Eigen::VectorXcd> v = stationary_distribution(mm);
  NumericVector results(v[0].size());
  if(v.size()> 1){
    stop("Markov matrix is not irreductible (many eigenvalues == 1).");
    return results; //useless in case of 'stop' above
  }
  if(v.size() == 0){
    stop("no eigenvector found.");
    return results;  //useless in case of 'stop' above
  } else {
    for(int i = 0; i < v[0].size(); i++)
      results[i] = v[0](i).real();
    return results;
  }
}

//' @description Calculates the exact p-value for short numerical Markov chains. Time computation can be too large for a sequence length of several thousands, specially for a data set.
//' @title Exact method for p-value [Markov chains]
//' @return A double representing the probability of a localScore as high as the one given as argument
//' @param m Transition matrix [matrix object]
//' @param sequence_length length of the sequence
//' @param localScore score for which the p-value should be calculated
//' @param sequence_min minimum score
//' @param sequence_max maximum score
//' @examples 
//' matrix = t(matrix(c(0.2, 0.3, 0.5, 0.3, 0.4, 0.3, 0.2, 0.4, 0.4), nrow = 3))
//' exact_mc(localScore = 12, m = matrix, sequence_length = 100, sequence_min = -1, sequence_max = 1)
//' exact_mc(localScore = 150, m = matrix, sequence_length = 1000, sequence_min = -1, sequence_max = 1)
//' @export
// [[Rcpp::export]]
double exact_mc(NumericMatrix m, int localScore, long sequence_length,  int sequence_min, int sequence_max){
  if(localScore<0)
    stop("[Invalid Input] local score must be positive.");
  if(sequence_length<= 0 )
    stop("[Invalid Input] sequence length must be positive.");
  if(sequence_max<=0)
    stop("[Invalid Input] sequence_max must be positive.");
  if(sequence_min>=0)
    stop("[Invalid Input] sequence_min must be negative.");
  MatrixXd mm(m.nrow(), m.ncol());
  for(int i = 0; i< m.nrow(); i++){
    for(int j = 0; j< m.ncol(); j++)
      mm(i, j) = m(i, j);
  }
  return mh_markov(localScore, mm, sequence_length, sequence_min, sequence_max);
}
