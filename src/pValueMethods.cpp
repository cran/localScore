//C++11 Standard

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#include "pValueMethods.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/unsupported/MatrixFunctions"

using namespace std;
using namespace Eigen;

bool sortByModule_desc(const std::complex<double> &one, const std::complex<double> &two) { return std::norm(one) > std::norm(two); }

bool sortByModule_asc(const std::complex<double> &one, const std::complex<double> &two) { return std::norm(one) < std::norm(two); }

bool sortByType(const std::complex<double> &one, const std::complex<double> &two) { return (one.imag() == 0.0 && two.imag() != 0.0) || ((one.imag() == 0.0 && two.imag() == 0.0) && std::norm(one) > std::norm(two)); }


//INPUT:    array probabilities for all x in scores from v to u
//OUTPUT:   array with polynomial coefficients belonging to the power of x corresponding to its index.
//The underlying formula is:
//
// P(x) = sum from i=1 to u(p_i*x^(u-i) + (p_0 -1)*x^u + sum from j = 1 to v (q_j * x^(u+j)
//
// ref (3) "An improved approximation for assessing the satistical significance of molecular sequence features",
// Applied Probability Trust 18 april 2003 by Mercier S, Cellier D, Charlot D University of Toulouse / Rouen
vector<double> calcul_poly (int u, int v, vector<double> probabilities){
  vector<double> polynom(u+v+1);                     //polynom: coefficents of polynom, descending order (x^n, x^(n-1),...x^0)
  polynom = probabilities;
  polynom[v] = probabilities[v]-1;
  return polynom;
}

//calculates the probability for each member of the score [-v..u] in sequence based on number of occurrences.
vector<double> calcul_probabilities (vector<int> sequence, int u, int v){
  vector<double> probabilities(u+v+1);
  int s ;
  //count occurences
  for(s=0;s<(int)sequence.size();s++)
    probabilities[v+sequence[s]]++;
  //devide by number of elements
  for(s=0;s<u+v;s++)
    probabilities[s]/=sequence.size();
  return probabilities;
}

//calculates the real or imaginary roots of the trinome x^2+px+q and returns the results as a couple of complex numbers
vector<complex<double>> eq_trinome(double p,double q)
{
  vector<complex<double>> roots;
  complex<double> x1 = complex<double>(0.0,0.0);
  complex<double> x2 = complex<double>(0.0,0.0);
  double delta,pr,pi;
  delta=p*p-4*q;
  pr=-p/2;
  pi=sqrt((double)abs(delta))/2;
  if(delta>=0)
  {
    x1.real(pr+pi);
    x2.real(pr-pi);
    x1.imag(0.0);
    x2.imag(0.0);
  }
  else
  {
    x1.real(pr);
    x2.real(pr);
    x1.imag(pi);
    x2.imag(-pi);
  }
  roots.push_back(x1);
  roots.push_back(x2);
  return roots;
}

//calculates the roots of a polynom with real coefficients a and degree n with the method of Bairstow and returns the results as complex numbers.
//copied and adapted from http://www.polytech-lille.fr/cours-algos-calcul-scientifique/progrceq.html
//n degree of polynom, a polynom coefficients IN DESCENDING ORDER, eps level of precision
vector<complex<double>> eq_bairstow(vector<double> polynom, double eps = 1e-15)
{
    int n = polynom.size()-1;             //degree of polynom is the number of coefficients -1
    unsigned int niterMax = 10000 ;
    unsigned int niter ;
    vector<complex<double>> roots;
    double * b = new double[n+1];                               //replaced
    double * c = new double[n+1];
    int i,j,k;
    double p,q,d,dp,dq,q0,p0;
    bool restart ;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> unif_random(0.0, 1.0);
    
    p0 = 1.0;
    q0 = -1.0;
    k=1;
    p=p0;
    q=q0;
    b[0]=polynom[0];
    c[0]=polynom[0];
    if (n>2) // polyn�me d'entree de degr� au moins > � 2
    {
        do
        {
            niter = 0;
            do
            {
                b[1] = polynom[1] - p * b[0];
                c[1] = b[1] - p * c[0];
                for (i = 2; i <= n; i++)
                    b[i] = polynom[i] - p * b[i - 1] - q * b[i - 2];
                for (j = 2; j <= n; j++)
                    c[j] = b[j] - p * c[j - 1] - q * c[j - 2];
                d = c[n - 2] * c[n - 2] + (b[n - 1] - c[n - 1]) * c[n - 3]; // XXX bug dans le cas n == 2
                dp = (b[n - 1] * c[n - 2] - b[n] * c[n - 3]) / d;           // XXX bug dans le cas n == 2
                dq = (b[n] * c[n - 2] + b[n - 1] * (b[n - 1] - c[n - 1])) / d;
                p += dp;
                q += dq;
                niter++ ;
                restart = (niter == niterMax) ;
                if (restart)  // non convergence after niterMax iteration, so restart with different start point
                {
//                    p = p0 + rand() ;  //new start point
//                    q = q0 - rand() ;  //new start point
                    p = p0 + unif_random(gen) ;  //new start point
                    q = q0 - unif_random(gen) ;  //new start point
                    dp = 10.0 * p ;  // to ensure the while condition is met
                    dq = 10.0 * q ; // to ensure the while condition is met
                    niter = 0 ;
                    restart = false ;
                }
            }
            while (((abs(dp) + abs(dq)) / (abs(p) + abs(q)) > eps) && (niter <= niterMax));

            //            printf("%de facteur quadratique:  p = %.5e   q = %.5e\n", (k + 1) / 2, p, q);
            //            printf("Reste: ");
            //            for (i = 0; i <= n - 2; i++)
            //                printf(" b%d = %.5e ", i, b[i]);
            //            printf("\n");
            vector<complex<double>> temp = eq_trinome(p, q);
            roots.insert(roots.end(), temp.begin(), temp.end());
/*            if (niter >= niterMax)
            {
                std::cerr << "ERROR : Non convergence of polynomial resolution by Bairstow algorithm\n" ;
                std::cerr << "Number of found root: " << roots.size() << "\n" ;
                for(std::vector<complex<double>>::iterator it = roots.begin(); it != roots.end(); ++it)
                {
                    std::cerr << *it << "\n" ;
                }
                exit(2) ;
            } */
            k += 2;
            n -= 2;
            for (i = 1; i <= n; i++)
                polynom[i] = b[i];
        }
        while (n > 2);
    }
    if (n == 2)
    {
        b[0]=polynom[0];
        b[1]=polynom[1];
        b[2]=polynom[2];

        p = b[1] / b[0];
        q = b[2] / b[0];
        vector<complex<double>> temp = eq_trinome(p, q);
        roots.insert(roots.end(), temp.begin(), temp.end());
    }
    else
    {
        b[0]=polynom[0];
        b[1]=polynom[1];
        roots.push_back(complex<double>(-b[1] / b[0], 0.0));
    }
//  } else {                //poly de degré = 1
//    // de la  forme ax+b=0 => x=-b/polynom
//    roots.push_back(complex<double>(-polynom[1] / polynom[0], 0.0));
//  }
    delete [] b;
    delete [] c;
    return roots;
}

// Verification of roots validity ()
// XXX TODO : retourner un entier valant 0 si OK et un entier en fonction de l'erreur.
bool verif_roots(vector <double> polynome, vector <complex<double>> roots, int u, int v, double eps)
{
  // 1. Number of roots
  if ((polynome.size()-1) != roots.size()){
    // Warning : racine multiple ou erreur
  }
  for(complex<double> root: roots){
    // calcul du polynome pris à la valeur 'root'
    // N.B. la boucle sur le vecteur polynome se fait de la fin vers le début (avec 'rbegin' et 'rend')
    // N.B.2 la boucle commence à l'avant dernier operande car le dernier est la constante.
    double cste = polynome.back();
    complex<double> value_poly(cste,0.0) ; // valeur de la constante du polynome
    complex<double> rootsPowI(1,0.0);
    for (auto it = polynome.rbegin()+1; it != polynome.rend(); ++it){
      rootsPowI = rootsPowI*root; // puissance i-1 de 'root'
      value_poly += *it * rootsPowI;
    }
    if (norm(value_poly)>eps){
      // root n'est pas une racine du polynome
      // les 3 lignes sont mise en commentaires (ne pas utiliser std::cerr avec R). A remplacer par un code d'erreur.
      // cerr << "root: " << root << "  ";
      // cerr << "value_poly: " << value_poly << "  ";
      // cerr << "norm(value_poly): " << norm(value_poly) << endl;
      return(false);
    }
  }
  return(true);
}

//splits the given array into two vectors: 1 contains all complex numbers with norm>1, and 2 all others.
Roots separate (vector<complex<double>> allroots){
  Roots results;
  for(complex<double> root: allroots){
    if(norm(root)<0.9999999){
      results.mod_smaller_one.push_back(root);
    } else {
      results.mod_greaterequal_one.push_back(root);
    }
  }
  return results;
}

RootTypes separateByType (vector<complex<double>> allroots){
  RootTypes results;
  for(complex<double> root: allroots){
    if(root.imag()==0.0){
      results.real_roots.push_back(root.real());
    } else {
      results.complex_roots.push_back(root);
    }
  }
  return results;
}

//puts real roots in the first half of the vector and complex roots in the second half
void orderbytype(vector<complex<double>> &roots){
  //keep order of modules....currently not
  sort(roots.begin(), roots.end(), sortByType);
}

//calculates probabilities for the first negative partial sum
vector<double> calcul_TabSmoins_cas_general_complexe(vector<complex<double>> roots_greater_one, int v){
  //-1 trier  racines de facon croissante
  sort(roots_greater_one.begin(), roots_greater_one.end(), sortByModule_asc);
  //0 declare variables
  MatrixXcd m(v,v);
  vector<double> TabSMoins;
  //1 Create Vandermonde Matrix for all roots modulus greater 1
  for(int s = 0; s<v; s++){
    for(int t = 0; t<v; t++)
      m(s, t) = pow(roots_greater_one[s], t+1);
  }

  //2 Decomposition of matrix
  Eigen::FullPivLU<Eigen::MatrixXcd> lu(m);

  Eigen::VectorXcd id = Eigen::VectorXcd(v);
  for(int s = 0; s < v; s++){
    complex<double> one = complex<double>(1.0, 0.0);
    id(s) = one;
  }
  //3 Resolution of system Ax = Id: Create Id, solve
  VectorXcd solutions = lu.solve(id);
  //4 add real part to vector and return result
  for(int s = 0; s<solutions.rows();s++)
    TabSMoins.push_back(solutions(s).real());
  return TabSMoins;
}

vector <double> calcul_deltaI_Complexe(RootTypes roots, vector<double> probabilities, int u, int v){
  vector <complex<double>> pkMoins1;
  complex <double> coeffC;
  double coeffR ;
  vector<double> resultats;
  unsigned int s, t;
  unsigned int idx;
  //0 Calculate size of the square matrix m
  int lSize; //Number of lines of m matrix
  int cSize; //Number of columns of m matrix
  cSize = roots.real_roots.size()+2*roots.complex_roots.size();
  lSize = roots.real_roots.size()+2*roots.complex_roots.size();
  
  // 1  Create Matrix m of size mSize*mSize
  Eigen::MatrixXd m(lSize,cSize);
  // 1.1 second row
  // real roots
  for(t = 0; t < roots.real_roots.size(); t++) {
    coeffR = p_1(probabilities, roots.real_roots[t], u, v);
    m(0, t) = coeffR;
  }
  //complex roots
  for (s = 0; s<roots.complex_roots.size() ; s++){
    coeffC = p_1(probabilities, roots.complex_roots[s], u, v);
    pkMoins1.push_back(coeffC);
    idx = roots.real_roots.size() + 2*s ;
    m(0, idx) = coeffC.real() ;
    m(0, idx+1) = coeffC.imag() ;
  }

  //1.2 row 3 to u-2
  // //1.2 row 3 to u-1
  //for(int k = 1; k < u; k++){ //XXX changement David (voir demande #6733 sur la forge DGA)
  for(int k = 1; k < u-1 ; k++){ //XXX changement David (voir demande #6733 sur la forge DGA)
    //real roots
    for(t = 0; t < roots.real_roots.size(); t++) {
      coeffR = roots.real_roots[t]*m(k-1, t);            //P_k(x) = x*P_(k-1)(x)+P[X=k]
      coeffR += probabilities[v+k+1];   //+1 for code confirmity with urbano code (XXX revoir +1 bizarre)
      m(k, t) = coeffR;
    }
    //complex roots (2 columns in m for each complex roots)
    for (s = 0; s<roots.complex_roots.size() ; s++){
      coeffC = roots.complex_roots[s]*pkMoins1[s];       //P_k(x) = x*P_(k-1)(x)+P[X=k] 
      coeffC.real(coeffC.real()+probabilities[v+k+1]);   //+1 for code confirmity with urbano code (XXX revoir +1 bizarre)
      pkMoins1[s] = coeffC;
      idx = roots.real_roots.size() + 2*s ;
      m(k, idx) = coeffC.real();
      m(k, idx+1) = coeffC.imag();     // XXX la partie imaginaire ne subit pas l'ajout de la proba P(X=k) comme la partie réelle ?
    }
  }
  // 1.1 last row (first row in Urbano rapport p.23)
  // real roots
  for(t = 0; t < roots.real_roots.size(); t++) {
    m(u-1, t) = 1. / (1. - roots.real_roots[t]);
  }
  //complex roots
  for (s = 0; s<roots.complex_roots.size() ; s++){
    coeffC = complex<double>(1.0, 0.0)/(complex<double>(1.0, 0.0) - roots.complex_roots[s]);
    idx = roots.real_roots.size() + 2*s ;
    m(u-1, idx) = coeffC.real();
    m(u-1, idx+1) = coeffC.imag();
  }
  
  // Supplementary rows to add constraints for complex root : delta_jR == delta_jI
  unsigned int idxLine;
  for (s = 0; s<roots.complex_roots.size() ; s++){
    idxLine = u + s; 
    for (t = 0; t<cSize; t++){ // initialisation de la ligne avec des 0
      m(idxLine,t) = 0.0 ;
    }
    idx = roots.real_roots.size() + 2*s ;
    m(idxLine, idx) = 1.0 ;
    m(idxLine, idx+1) = -1.0 ;
  }  
  
  // 2  Création du vecteur b pour la résolution de m*x = b
  Eigen::VectorXd b = Eigen::VectorXd(lSize);

  // Attention ci-dessous : la première ligne de la matrice M du rapport urbano est ici placée en dernier.
  for(s = 0; s<lSize; s++)
    b(s) = 0.0;
  b(u-1) = 1.0;

  // 3  Decomposition LU
  Eigen::FullPivLU<Eigen::MatrixXd> lu(m);
  // 4  Résolution système
  VectorXd solutions = lu.solve(b);
  // 5  Convertion en vecteur
  for(int s = 0; s < solutions.rows(); s++)
    resultats.push_back(solutions[s]);

  return resultats;
}


//Rapport de stage 2002 Anne-Benedicte Urbano, "Conception d'un logiciel à but biologique", page 24
complex<double> p_1(vector<double> distribution, complex<double> x, int u, int v){
  complex<double> result;
  result = x*(distribution[v]-1);                                         //+ (P[X=0]-1)*x
  result.real(result.real()+ distribution[v+1]);                         //P[X=1]
  for(int s=-v;s<=-1;s++)
    result=result+(distribution[s+v]*pow(x, -s+1));     //sum from s=1 to u of (P[X = -s]*x^(s+1)
  return(result);
}
double p_1(vector<double> distribution, double x, int u, int v){
  double result;
  result = x*(distribution[v]-1);                                         //+ (P[X=0]-1)*x
  result = result + distribution[v+1];                         //P[X=1]
  for(int s=-v;s<=-1;s++)
    result=result+(distribution[s+v]*pow(x, -s+1));     //sum from s=1 to u of (P[X = -s]*x^(s+1)
  return(result);
}

// Return P(M = k) = \sum_i delta_i*R_i^k : Distribution du maximum des sommes partiels.
// Formule (4) page 5 du papier MCC (2003)
double calcul_probMaxPartialSum(int k, vector<double> distribution, int u, int v){
  //1  Calcul polynome
  vector <double> polynome = calcul_poly(u, v, distribution);
  //2  Calcul racines
  vector <complex<double>> roots = eq_bairstow(polynome);
  // Roots validity verification (last parameter is the limit value of norm of the calculus of the polynome for each root)
  if (!verif_roots(polynome, roots, u, v, 1e-10)){
    // cerr << "ERROR probMaxPartialSum (cpp): 'roots' are not valid for the 'polynome'" <<endl;
    return -1.0;
  }
  //3  Separation des racines avec un module < 1
  Roots roots_separated = separate(roots);
  if(roots_separated.mod_smaller_one.size() != u){
    // cerr << "ERROR probMaxPartialSum (cpp): roots_separated.mod_smaller_one.size() != u" <<endl;
    return -1.0;
  }
  if(roots_separated.mod_greaterequal_one.size() != v){
    // cerr << "ERROR probMaxPartialSum (cpp): roots_separated.mod_greaterequal_one.size() != v" <<endl;
    return -1.0;
  }
  //4  Triage des racines avec un module < 1 par ordre decroissant
  sort(roots_separated.mod_smaller_one.begin(), roots_separated.mod_smaller_one.end(), sortByModule_desc);
  //5  Put real roots into first part of vector and complex ones into second part
  orderbytype(roots_separated.mod_smaller_one);
  //5.1 Put real and complex roots in different vectors
  RootTypes roots_mod_less_one = separateByType(roots_separated.mod_smaller_one);
  //5.2 Remove doubles from complex part
  roots_mod_less_one.complex_roots.erase( unique( roots_mod_less_one.complex_roots.begin(), roots_mod_less_one.complex_roots.end() ), roots_mod_less_one.complex_roots.end() );
  //6  Calcul deltaI
  //vector <double> deltaI = calcul_deltaI_Complexe(roots_mod_less_one, distribution, u, v);
  vector <double> deltaI = calcul_deltaI_Complexe(roots_mod_less_one, distribution, u, v);
  
  //7.1 calcul avec les racines réelles
  double somme1 = 0.0;
  double Ri = 0.0;
  for (unsigned int s = 0 ; s < roots_mod_less_one.real_roots.size(); s++) {
    Ri = roots_mod_less_one.real_roots[s];
    somme1 += deltaI[s]*std::pow(Ri,k) ;
  }

  // //7.2 calcul avec les racines complexes
  double somme2 = 0.0 ;
  complex<double> Cj ;
  complex<double> Aux ;

  for (unsigned int s=0 ; s<roots_mod_less_one.complex_roots.size() ; s++){
    Cj = roots_mod_less_one.complex_roots[s];
    Aux = std::pow(Cj,k) ;
    somme2+= Aux.real()*deltaI[roots_mod_less_one.real_roots.size()+s];
    somme2+= Aux.imag()*deltaI[roots_mod_less_one.real_roots.size()+s];
  }
  
  return(somme1+somme2);
}

// Return P(. >=localScore) with Karlin-Altschul approximation pvalue = 1.0-exp(-K_star*sequence_length*exp(-lambda*(localScore-1)))
// localscore : value at calculation of the probability
// u : maximum (individual) score
// v : -minimum (individual) score
// distribution : vector of probability of the (individual) scores (from -v to u, length u+v+1)
// sequence_length : length of the sequence of scores.
// Code error return (< 0.0) :
//    -1.0 : problème with polynomial roots
//    -2.0 : distribution vector size and u, v are not compatible
double calcul_karlin(int localScore, vector<double> distribution, int u, int v, long sequence_length){
  // Some checks of error
  if (distribution.size() != (u+v+1)) {
    //cerr << "ERROR calcul_karlin (cpp): u and/or v are not compatible with the size of 'distribution'" <<endl;
    return -2.0;
  }
  //    Algorithme de Karlin
  // 0   Trivial cases localScore <=0 and =1
  if (localScore<=0)
    return 1.0;
  if (localScore==1)
    return 1.0 ;
  
  // 1   Calcul Polynome
  vector <double> polynome = calcul_poly(u, v, distribution);

  // 2   Calcul Racines
  vector <complex<double>> roots = eq_bairstow(polynome);
  // Roots validity verification (last parameter is the limit value of norm of the calculus of the polynome for each root)
  if (!verif_roots(polynome, roots, u, v, 1e-10)){
    // cerr ne doit pas être utilisé avec R. Mise en commentaires des lignes ci-dessous.
    // cerr << "ERROR calcul_karlin (cpp): 'roots' are not valid for the 'polynome'" <<endl;
    // cerr << "karlin() function cannot be used in your case. Check the documentation of 'karlin()' for details" <<endl;
    // cerr << "You could try to change your scoring discretisation step or use karlinMonteCarlo()" << endl;
    // for (auto root: roots) {
    //   cerr << root << endl ;
    // }
    return -1.0;
  }

  // 3 Separation des racines avec un module <1
  Roots roots_separated = separate(roots);
  if(roots_separated.mod_smaller_one.size() != u){
    // cerr << "ERROR calcul_karlin (cpp): roots_separated.mod_smaller_one.size() != u" <<endl;
    return -1.0;
  }
  if(roots_separated.mod_greaterequal_one.size() != v){
    // cerr << "ERROR calcul_karlin (cpp): roots_separated.mod_greaterequal_one.size() != v" <<endl;
    return -1.0;
  }
  
  // 4   Triage des Racines avec un module <1 par module decroissante
  sort(roots_separated.mod_smaller_one.begin(), roots_separated.mod_smaller_one.end(), sortByModule_desc);

  // 5   Put real roots into first part of vector and complex ones into second part
  orderbytype(roots_separated.mod_smaller_one);

  // 6   Calcul Lambda
  double lambda = log(1.0/roots_separated.mod_smaller_one[0].real());
  // 7   Substitution x pour formule de Karlin (log(n)/lambda+x)
  // double x = localScore- floor((log(sequence_length)/lambda))-1; 
  // 8   Calcul Esperance de probabilités: p(x)*{-v..u}
  double E = real(-v);
  for(unsigned int s = 0; s< distribution.size(); s++)
    E+=distribution[s]*real(s);
  // 9   Calcul probabilités S-
  vector<double> probaSmoins = calcul_TabSmoins_cas_general_complexe(roots_separated.mod_greaterequal_one, v);
  // 10  Calcul Esperance Probabilités S-
  double ESmoins = 0.0;
  for(unsigned int s = 0; s< probaSmoins.size(); s++)
    ESmoins+=probaSmoins[s]*(-real(s)-1.0);
  // 11   Calcul mu
  double mu = ESmoins/E;
  // 12   Calcul K*
  //  12.1 Calcul Esperance numerateur E[e^(lambda*S-)]
  double E1 = 0.0;
  for(unsigned int s = 0; s< probaSmoins.size(); s++)
    E1+=probaSmoins[s]*exp(lambda*(-real(s)-1));  // XXX le -1 me semble de trop.
  //   12.2 Calcul numerateur
  double numerateur = pow(1-E1,2);
  //   12.3 Calcul Esperance denominateur
  //        D'abord le calcul de E[X*e^(lambda*X)]
  double E2 = 0.0;
  for(unsigned int s = 0; s< distribution.size(); s++)
    {
    double x = double(s)-double(v) ;
    E2+=distribution[s]*x*exp(lambda*x);
  }
  //   12.4 Calcul denominateur
    double denominateur = pow(mu, 2)*E2*(exp(lambda)-1.0);
  //   12.5 K*
  double K_star = numerateur/denominateur;
  // 13   Calcul p-Value
  // localScore-1 needed because we want return P(. >=localScore) = 1-P(.<localScore) = 1-P(.<=localScore-1) [because localScore is integer]
  return (1.0-exp(-K_star*sequence_length*exp(-lambda*(localScore-1)))); 
  //  return (1.0-exp(-K_star*exp(-lambda*x)));
}

//INPUT:    array probabilities for all x in scores de v to u
//OUTPUT:   array with polynomial coefficients belonging to the power of x corresponding to its index.
//The underlying formula is:
//
// P(x) = sum from i=1 to u(p_i*x^(u-i) + (p_0 -1)*x^u + sum from j = 1 to v (q_j * x^(u+j)
//
// ref (3) "An improved Approximation for assessing the satistical significance of molecular sequence features",
// Applied Probability Trust 18 april 2003 by Mercier S, Cellier D, Charlot D University of Toulouse / Rouen
double calcul_mcc(int localScore, vector<double> distribution, int u, int v, long sequence_length){
  unsigned int s;
  
  //0   Trivial cases localScore <=0 and =1
  if (localScore<=0)
    return 1.0;
  // if (localScore==1) NON, pas trivial
  //   return 1.0 ;
  
  //0.1 localScore shift
  // On veut P(H_n>=localScore) = 1 - P(H_n<= localScore-1) car distribution discrète
  localScore = localScore - 1;

    //1  Calcul polynome
  vector <double> polynome = calcul_poly(u, v, distribution);
  //2  Calcul racines
  vector <complex<double>> roots = eq_bairstow(polynome);
  // Roots validity verification (last parameter is the limit value of norm of the calculus of the polynome for each root)
  if (!verif_roots(polynome, roots, u, v, 1e-10)){
    // cerr << "ERROR calcul_mcc (cpp): 'roots' are not valid for the 'polynome'" <<endl;
    // cerr << "mcc() function cannot be used in your case. Check the documentation of 'mcc()' for details" <<endl;
    // cerr << "You could try to change your scoring discretisation step or use karlinMonteCarlo()" << endl;
    // for (auto root: roots) {
    //   cerr << root << endl ;
    // }
    return -1.0;
  }
  //3  Separation des racines avec un module < 1
  Roots roots_separated = separate(roots);
  if(roots_separated.mod_smaller_one.size() != u){
    // cerr << "ERROR calcul_mcc (cpp): roots_separated.mod_smaller_one.size() != u" <<endl;
    return -1.0;
  }
  if(roots_separated.mod_greaterequal_one.size() != v){
    // cerr << "ERROR calcul_mcc (cpp): roots_separated.mod_greaterequal_one.size() != v" <<endl;
    return -1.0;
  }
  //4  Triage des racines avec un module < 1 par ordre decroissant
  sort(roots_separated.mod_smaller_one.begin(), roots_separated.mod_smaller_one.end(), sortByModule_desc);
  //5  Put real roots into first part of vector and complex ones into second part
  orderbytype(roots_separated.mod_smaller_one);
  //5.1 Put real and complex roots in different vectors
  RootTypes roots_mod_less_one = separateByType(roots_separated.mod_smaller_one);
  //5.2 Remove doubles from complex part
  roots_mod_less_one.complex_roots.erase( unique( roots_mod_less_one.complex_roots.begin(), roots_mod_less_one.complex_roots.end() ), roots_mod_less_one.complex_roots.end() );
  //6  Calcul deltaI
  vector <double> deltaI = calcul_deltaI_Complexe(roots_mod_less_one, distribution, u, v);
  
  //7  Calcul S-
  vector<double> probaSmoins = calcul_TabSmoins_cas_general_complexe(roots_separated.mod_greaterequal_one, v);
  //8  Calcul E[X]
  double E = real(-v);
  for(unsigned int s = 0; s< distribution.size(); s++)
    E+=distribution[s]*real(s);
  
  //9  Calcul E[S-]
  double ESmoins = 0.0;
  for(unsigned int s = 0; s< probaSmoins.size(); s++)
    ESmoins+=probaSmoins[s]*(-real(s)-1.0);

  //10  Calcul mu
  double mu = ESmoins/E;
  //11 Calcul lambda
  double R = roots_separated.mod_smaller_one[0].real();
  double lambda = -log(R);
  //12 Substitution x
  //double x = localScore - floor((log(sequence_length)/lambda));  //inutile car directement mis dans la formule (qui se simplifie)
  //13 Calcul Somme 1 > Ref p. 12 Rapport Anne Benedicte Urbano 2002
  double somme1 = 0.0;
  {
    double ERi = 0.0;
    double Ri = 0.0;
    double coeffi2 = 0.0;
    double coeffi3 = 0.0;
    for (s = 0; s < roots_mod_less_one.real_roots.size(); s++) {
      Ri = roots_mod_less_one.real_roots[s];
      ERi = 0.0;
      for (int t = 0; t < v; t++)
        ERi += probaSmoins[t] * pow(Ri, t + 1);
      coeffi2 = Ri * (1 - ERi) / (1 - Ri);
      coeffi3 = pow(Ri, localScore); // test conformément au mail de Sabine du 18/09/2023 à 15H30
      //      coeffi3 = pow(Ri, localScore+1);
      somme1 += deltaI[s]*coeffi2*coeffi3; 
    }
  }
  
  //14 Calcul Somme 2
  double somme2 = 0.0;
  {
    complex<double> coeffj2, coeffj3;
    unsigned int idxDeltaI ;
    complex<double> Cj, ECj, Aux;
    for(s=0; s<roots_mod_less_one.complex_roots.size(); s++){
      Cj = roots_mod_less_one.complex_roots[s];
      ECj.real(0.0);
      ECj.imag(0.0);
      for(int t = 0;t<v; t++)
        ECj = ECj + (pow(Cj, t+1) * probaSmoins[t]);
      coeffj2 = Cj * (complex<double>(1.0,0.0) - ECj)/(complex<double>(1.0,0.0)-Cj);
      coeffj3 = pow(Cj, localScore); // test conformément au mail de Sabine du 18/09 à 15H30
      Aux = coeffj2*coeffj3;
      idxDeltaI = roots_mod_less_one.real_roots.size() + 2*s ;
      somme2+= Aux.real()*deltaI[idxDeltaI];
      somme2+= Aux.imag()*deltaI[idxDeltaI+1]; // ajout de '-1' par rapport au code initial v1.0.8
    }
  }

  //15 Calcul MCC
  double MCC = pow(1.0-somme1-somme2, (sequence_length/mu)+1.0);
  return (1.0-MCC);
}


MatrixXd ind(MatrixXd input, int power){
  MatrixXd result(input.cols(), input.rows());
  result.setIdentity();
  if(power == 0)
    return result;
  if( power == 1)
    return input;
  while( power > 1){
    if(power%2 == 0){
      input = input*input;
      power = power / 2;
    } else {
      result = input*result;
      input = input*input;
      power = (power - 1)/2;
    }
  }
  return input*result;
}

vector<VectorXcd> stationary_distribution( MatrixXd transitionMatrix){
  vector<double> result;
  double eps = 1e-12;
  
  //1 Calcul vecteurs propres droits et valeurs propres
  EigenSolver<MatrixXd> es(transitionMatrix);
  //2 Determination du nombre des vecteurs propres avec longeur 1: si 1, bien, si 2 qu'est-ce qu'on fait? Si zero?
  vector<int> indices;
  for (int i = 0; i<es.eigenvalues().rows();i++){
    if(std::fabs(es.eigenvalues()[i].real()-1.0)<eps)
      indices.push_back(i);
  }
  //3 Calcul vecteurs propres gauches: inverse des vecteurs droits
  MatrixXcd inv = es.eigenvectors().inverse();

  //4 Calcul des distributions stationnaires
  vector<VectorXcd> distributions;
  for(unsigned int i = 0; i < indices.size(); i++){
    VectorXcd distribution = inv.row(indices[i]);
    distribution = distribution / distribution.sum();
    distributions.push_back(distribution);
  }

  return distributions;
}

double calcul_daudin(int a, int n, vector <double> score_probabilities, int smin, int smax) {
  MatrixXd Pi = ind(creation_pi_new(score_probabilities,a, smin, smax), n);
  return Pi(0, a);
}

MatrixXd creation_pi_new(vector <double> score_probabilities, int a, int smin, int smax){
  MatrixXd pi(a+1, a+1);
  //Fill matrix pi
  for(int h = 0; h < a; h++){                                 //row by row
    pi(h, 0) = f(-h, score_probabilities, smin, smax);      //first column
    for(int l = 1; l < a; l++)                            //per row of matrix
      pi(h, l) = p(l-h, score_probabilities, smin, smax);
    pi(h , a) = 1 - f(a-h-1, score_probabilities, smin, smax);
  }
  //last row
  for(int l = 0; l < a; l++)                                //per row of matrix
    pi(a, l) = 0.0;
  pi(a, a) = 1.0;
  return pi;
}

//Daudin page 7
double f(int k, vector<double> probabilities, int smin, int smax){
  if (k<smin)
    return 0.0;
  else if (k>smax)
    return 1.0;
  else {
    int index = k-smin;
    double result = 0.0;
    for(int t = 0; t <= index; t++)
      result += probabilities[t];
    return result;
  }
}
//Daudin page 7
double p(int k, vector<double> probabilities, int smin, int smax){
  if (k<smin || k>smax)
    return 0.0;
  else {
    int index = k-smin;
    return probabilities[index];
  }
}

double mh_markov(int localscore, MatrixXd transitionMatrix, VectorXi scoreValues, int sequence_length, VectorXd probInitial){
  
  int i, l, c; //used for iterating
  
  int s_min = scoreValues.minCoeff(); 
  int nbScores = scoreValues.size();

  // Créer la table de hashage (score local, score) -> n° ligne de la grande matrice P
  std::map<std::pair<int,int>, long> Pindex;
  long idx = 0;
  for (int sl=0 ; sl<=localscore ; sl++){
    //for(auto score : scoreValues) {
    for(i=0 ; i<nbScores ; i++) {
        Pindex[std::make_pair(sl,scoreValues(i))] = idx ;
      idx++ ;
    }
  }

  // table de hashage pour la matrice de transition (score) -> n° ligne ou colonne de la matrice TransitionMatrix
  std::map<int, int> itransitionMatrix ;
  for(i=0 ; i<nbScores ; i++) {
      itransitionMatrix[scoreValues(i)] =i;
  }
  
  //2 Create Matrix from transition matrix where m = transitionMatrix - x
  long matrix_length = Pindex.size();
  MatrixXd P(matrix_length, matrix_length) ;
  P.setZero() ;
  //2.1 fill Matrix
  long idxLine;
  long idxColumn;
  int slc, scorec, sll, scorel;
  for (slc=0 ; slc<=localscore ; slc++){
    for(c=0 ; c<nbScores ; c++) {
      scorec = scoreValues(c);
      idxColumn = Pindex[std::make_pair(slc,scorec)] ;
      for (sll=0 ; sll<=localscore ; sll++){
        for(l=0 ; l<nbScores ; l++) {
          scorel = scoreValues(l);
          idxLine = Pindex[std::make_pair(sll,scorel)] ;
          if(sll == localscore && slc != localscore)
            P(idxLine, idxColumn) = 0.0;
          else if(sll == localscore && slc == localscore)
            P(idxLine, idxColumn) = transitionMatrix(itransitionMatrix[scorel], itransitionMatrix[scorec]); 
          else if(slc == 0 && sll+scorel <= 0) 
            P(idxLine, idxColumn) = transitionMatrix(itransitionMatrix[scorel], itransitionMatrix[scorec]);
          else if (sll+scorel>=1 && sll+scorel<localscore && slc == sll+scorel)
            P(idxLine, idxColumn) = transitionMatrix(itransitionMatrix[scorel],itransitionMatrix[scorec]);
          else if (sll+scorel>=localscore && slc == localscore)
            P(idxLine, idxColumn) = transitionMatrix(itransitionMatrix[scorel], itransitionMatrix[scorec]);
          else
            P(idxLine, idxColumn) = 0.0;
        }
      }
    }
  }

  //3 Calcul p-valeur
  //3.1 Matrix ^ length of sequence...
  //  MatrixXd Pn = ind(P, sequence_length);
  MatrixXd Pn(matrix_length, matrix_length) ;
  Pn = P.pow(sequence_length);

  //3.2 sum
  double p_value = 0.0;
  int index_Pn_0_u;
  int index_Pn_a_v;
  double lambda_u;
  for(l=0 ; l<nbScores ; l++) {
    scorel = scoreValues(l);
    lambda_u = probInitial(itransitionMatrix[scorel]);
    index_Pn_0_u = Pindex[std::make_pair(0,scorel)] ;
    for(c=0 ; c<nbScores ; c++) {
      scorec = scoreValues(c);
      index_Pn_a_v = Pindex[std::make_pair(localscore,scorec)] ;
      p_value += lambda_u*Pn(index_Pn_0_u, index_Pn_a_v);
    }
  }
  return p_value;
}
