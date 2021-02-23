#ifndef PVALUEMETHODS_H
#define PVALUEMETHODS_H

//#include "Eigen/src/Core/util/DisableStupidWarnings.h"
#include "Eigen/Dense"
//#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

typedef struct {
  std::vector<std::complex<double>> mod_smaller_one;
  std::vector<std::complex<double>> mod_greaterequal_one;
} Roots;

typedef struct {
  std::vector<double> real_roots;
  std::vector<std::complex<double>> complex_roots;
} RootTypes;

typedef struct {
  long t1;
  long t2;
} Tuple;

std::vector<double> calcul_poly (int u, int v, std::vector<double> probabilities);

//calculates the probability for each member of the score [-v..u] in sequence based on number of occurrences.
std::vector<double> calcul_probabilities (std::vector<int> sequence, int u, int v);

//calculates the reel or imaginary roots of the trinome x^2+px+q and saves them in the passed complex pointer fields
std::vector<std::complex<double> > eq_trinome(double p,double q);

//calculates real and complex roots of polynom with coefficients a and degre n with level of precision e
std::vector<std::complex<double> > eq_bairstow(std::vector<double> polynom, double eps);

void printArray (int length, double* array);

void testing ();

Roots separate (std::vector<std::complex<double> >);

std::vector<double> calcul_TabSmoins_cas_general_complexe(std::vector<std::complex<double> > roots_greater_one, int v);

double calcul_karlin(int localScore, std::vector<double> distribution, int u, int v, long sequence_length);

double calcul_mcc(int localScore, std::vector<double> distribution, int u, int v, long sequence_length);

std::complex<double> p_1(std::vector<double> distribution, std::complex<double> x, int u, int v);

double p_1(std::vector<double> distribution, double x, int u, int v);

std::vector <double> calcul_deltaI_Complexe(std::vector<std::complex<double> > complexRoots, std::vector<double> probabilities, int u, int v);

Eigen::MatrixXd ind(Eigen::MatrixXd input, int power);

double calcul_daudin(int a, int n, std::vector <double> score_probabilities, int smin, int smax);

Eigen::MatrixXd creation_pi(std::vector <double> score_probabilities,int a, int smin, int smax);

double f(int k, std::vector<double> probabilities, int smin, int smax);

double p(int k, std::vector<double> probabilities, int smin, int smax);

std::vector<Eigen::VectorXcd> stationary_distribution( Eigen::MatrixXd transitionMatrix);

Eigen::MatrixXd creation_pi_new(std::vector <double> score_probabilities, int a, int smin, int smax);

//maps tuple to integer index for matrices with different bases: ex: tuple (i,j) with (i,j) in [-1,1] X [3,5] yields an index from
//0 to 8 with 0 -> (-1, 3), 1 -> (-1, 4).... 8 -> (1,5)
int tuple2index(Tuple i, Tuple intervall_t1, Tuple intervall_t2);

Tuple index2tuple(int i, Tuple intervall_t1, Tuple intervall_t2);

double mh_markov(int localscore, Eigen::MatrixXd transitionMatrix, long sequence_length,  int s_min, int s_max);

int card(Tuple interval);
//deprecated:
/*
void indienne(double *MatE, int tE, int nE, double *SortieExp);
double calcul_daudinOLD(int a, int n, std::vector <double> score_probabilities, int smin, int smax);
void creationPI(double *distC, int mC, int tailleC, double *PIC,int a, int smin, int smax);
void decompose(int ka, int pe, int *alpha);
void mult(double *A, int nbLigA, int nbColA,double *B, int nbLigB, int nbColB, double *Res);
void carre(double *MatC, int taye, double *Calcul);
*/
#endif
