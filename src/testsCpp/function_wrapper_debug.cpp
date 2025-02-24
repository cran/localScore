#include <Rcpp.h>
#include "pValueMethods.h"
//#include "malloc.h"

#include "Eigen/Dense"

using namespace Rcpp;
using Eigen::MatrixXd;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


//vector <double> calcul_karlin_parameters(vector<double> distribution, int u, int v) {
// [[Rcpp::export]]
std::vector <double> icalcul_karlin_parameters(NumericVector distribution, int u, int v) {
  return calcul_karlin_parameters(as<std::vector<double>>(distribution), u, v) ;
}
// 
// karlin(segment_scores$value[i], length(chro$score1),score_probabilities = score1.prob.ext, sequence_min = scoremin, sequence_max = scoremax)      

// [[Rcpp::export]]
std::vector <double> icalcul_poly (int u, int v, NumericVector probabilities){
  return calcul_poly (u, v, as<std::vector<double>>(probabilities));
}

// [[Rcpp::export]]
std::vector<std::complex<double>> ieq_bairstow(NumericVector polynom, double eps = 1e-15){
  std::vector<double> polycoef = Rcpp::as<std::vector<double>>(polynom) ;
  return eq_bairstow(polycoef, eps) ;
}
// [[Rcpp::export]]
std::vector<std::complex<double>> ieq_bairstow_Rpolyroot(NumericVector polynom, double eps = 1e-15){
  std::vector<double> polycoef = Rcpp::as<std::vector<double>>(polynom) ;
  return eq_bairstow_Rpolyroot(polycoef, eps) ;
}

// [[Rcpp::export]]
std::vector<std::complex<double>> ieq_bairstow_eigen(NumericVector polynom, double eps = 1e-15){
  std::vector<double> polycoef = Rcpp::as<std::vector<double>>(polynom) ;
  return eq_bairstow_eigen(polycoef, eps) ;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
scoremax = 13 ; scoremin = -2
## note : use 'dput(myRvec)' to get this format output
score1.prob.ext = c(`-2` = 0.0498654613882978, `-1` = 0.735683157172755, `0` = 0.14835546332284, 
                    `1` = 0.043958424930972, `2` = 0.0142114982149452, `3` = 0.00559469583721707, 
                    `4` = 0.00154059900459014, `5` = 0.000411529871089147, `6` = 0.00020259932115158, 
                    `7` = 9.07476125991453e-05, `8` = 4.71324809623468e-05, `9` = 2.18075658183993e-05, 
                    `10` = 9.84857811153515e-06, `11` = 5.62775892087723e-06, `12` = 7.03469865109653e-07, 
                    `13` = 7.03469865109653e-07)  #chro 1

score1.prob.ext = c(`-2` = 0.0681561425858712, `-1` = 0.759474076388722, `0` = 0.121671514073504, 
  `1` = 0.0310891041511429, `2` = 0.0125039059683196, `3` = 0.00563901641708531, 
  `4` = 0.0010912674566738, `5` = 0.000242770954017739, `6` = 7.93212027978752e-05, 
  `7` = 3.60550921808523e-05, `8` = 1.44220368723409e-05, `9` = 2.40367281205682e-06) #chro 18 "NC_010460.4"
polynom = icalcul_poly(scoremax, -scoremin, score1.prob.ext)
length(polynom)
roots.bairstow = ieq_bairstow(polynom)

#### test root eigen
roots.eigen = ieq_bairstow_eigen(polynom)
#karlin.params = icalcul_karlin_parameters(score1.prob.ext, scoremax, -scoremin)

#poly5 = c(2.0, 1.0, 0.0, 0.0, -2.0, -1.0)
#polynom = poly5

polynomrev= rev(polynom)
x=roots.eigen[length(roots.eigen)]
#x=roots.polyroot[15]
#x=-15.8353342
res = polynomrev[1]
print(length(polynomrev))
for (i in 2:length(polynomrev)) {
  res = res + polynomrev[i]*x^(i-1)
}
print(res)  # doit être égale à 0

polynomrev= rev(polynom)
res = rep(NA,length(roots.eigen))
for (j in 1:length(roots.eigen)) {
  x=roots.eigen[j]
  res[j] = polynomrev[1]
  for (i in 2:length(polynomrev)) {
    res[j] = res[j] + polynomrev[i]*x^(i-1)
  }
}
print(Mod(res))


### Test ieq_bairstow_Rpolyroot
roots.polyr = ieq_bairstow_Rpolyroot(polynom)

polynomrev= rev(polynom)
x=roots.polyr[length(roots.polyr)]
#x=roots.polyroot[15]
#x=-15.8353342
res = polynomrev[1]
print(length(polynomrev))
for (i in 2:length(polynomrev)) {
  res = res + polynomrev[i]*x^(i-1)
}
print(res)  # doit être égale à 0

polynomrev= rev(polynom)
res = rep(NA,length(roots.polyr))
for (j in 1:length(roots.polyr)) {
  x=roots.polyr[j]
  res[j] = polynomrev[1]
  for (i in 2:length(polynomrev)) {
    res[j] = res[j] + polynomrev[i]*x^(i-1)
  }
}
print(Mod(res))
### XXX pour écrire automatique =ment une fonction.
fct = as.character(polynomrev[1])
for (i in 2:length(polynomrev)) {
  fct = paste(fct, " + ", as.character(polynomrev[i]),"*x^",i-1,sep="")
}
fct
tmp=parse(text=fct)

### Test native R polyroot
roots.polyroot = polyroot(rev(polynom))
polynomrev= rev(polynom)
x=roots.polyroot[length(roots.polyroot)]
res = polynomrev[1]
for (i in 2:length(polynomrev)) {
  res = res + polynomrev[i]*x^(i-1)
}
print(res)  # doit être égale à 0
polynomrev= rev(polynom)
res = rep(NA,length(roots.polyroot))
for (j in 1:length(roots.polyroot)) {
  x=roots.polyroot[j]
  res[j] = polynomrev[1]
  for (i in 2:length(polynomrev)) {
    res[j] = res[j] + polynomrev[i]*x^(i-1)
  }
}
print(Mod(res))

# Degree > 2:  2x^5+x^4-2x-1 = 0   Expected solution : -1, 1, -1/2, i, -i
poly5 = c(2.0, 1.0, 0.0, 0.0, -2.0, -1.0)
#roots.test = polyroot(rev(poly5r))
roots.test = ieq_bairstow_Rpolyroot(poly5)
f = function(x) {2*x^5+x^4-2*x-1}
print(roots.test) # Expected solution : -1, 1, -1/2, i, -i
print(f(roots.test)) # Expected solution : -1, 1, -1/2, i, -i
sum(f(roots.test))

prob = c(0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08)
u = 4 ; v = 5
icalcul_poly(u, v, prob) # attended_result = c(0.1, 0.05, 0.02, 0.03, 0.2, -0.7, 0.04, 0.06, 0.12, 0.08)
*/
