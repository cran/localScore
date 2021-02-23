#define BOOST_TEST_MODULE First_TestSuite
#include "../Eigen/src/Core/util/DisableStupidWarnings.h"
#include <boost/test/included/unit_test.hpp>
#include "../pValueMethods.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
//#include <src/karlinfuncs.h>
//#include <src/pValueMethods.h>
#include "../Eigen/Dense"

using namespace std;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE( auxiliary_programs_tests )

    BOOST_AUTO_TEST_CASE( polynome ){
        int u = 4;
        int v = 5;
        vector <double> probabilities = {0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08};
        vector<double> poly = calcul_poly(u, v, probabilities);
        // P(x) = sum from i=1 to u(p_i*x^(u-i) + (p_0 -1)*x^u + sum from j = 1 to v (q_j * x^(u+j)
        //
        //ref (3) "An improved Approximation for assessing the satistical significance of molecular sequence features",
        // Applied Probability Trust 18 april 2003 by Mercier S, Cellier D, Charlot D University of Toulouse / Rouen
        vector<double> attended_result = {0.1, 0.05, 0.02, 0.03, 0.2, -0.7, 0.04, 0.06, 0.12, 0.08};
        BOOST_CHECK_EQUAL_COLLECTIONS(poly.begin(), poly.end(),
                                      attended_result.begin(), attended_result.end());
    }
    BOOST_AUTO_TEST_CASE( indienne ){
        MatrixXd mat(2, 2);
        mat << 1, 2, 3, 4;
        MatrixXd result;
        result = ind(mat, 5);
    BOOST_TEST_MESSAGE("Matrix Result: "<<result);
    }
    BOOST_AUTO_TEST_CASE( stationarydistribution ){
        Matrix3d m;
        m << 0.6, 0.1, 0.3,
                0.1, 0.7, 0.2,
                0.2,0.2,0.6;
        BOOST_TEST_MESSAGE("Stationary Distribution of \n"<<m<<": \n"<<stationary_distribution(m)[0]);
    }
    BOOST_AUTO_TEST_CASE( indexconversion_index_to_tuple ){
            Tuple intervall1 = {-2, 3};
            Tuple intervall2 = {6, 8};
            BOOST_TEST_MESSAGE("Interval1: -2, 3, Interval2: 6, 8");
            for(long i = 0; i< (intervall1.t2 - intervall1.t1+1)*(intervall2.t2 - intervall2.t1+1); i++) {
                Tuple tup = index2tuple(i,intervall1, intervall2);
                BOOST_TEST_MESSAGE(i << " -> (" <<tup.t1<<", "<<tup.t2<<")");
            }
    }
    BOOST_AUTO_TEST_CASE( indexconversion_tuple_to_index ){
        Tuple intervall1 = {-4, -1};
        Tuple intervall2 = {3, 5};
        BOOST_TEST_MESSAGE("Interval1: -4, -1, Interval2: 3, 5");
        for(long i = intervall1.t1; i<= intervall1.t2; i++) {
            for(long j = intervall2.t1; j<= intervall2.t2;j++) {
                Tuple tup = {i,j};
                BOOST_TEST_MESSAGE("(" <<tup.t1<<", "<<tup.t2<<") -> "<< tuple2index(tup, intervall1, intervall2));
            }
        }
    }
    BOOST_AUTO_TEST_CASE( roots ) {
        //Degree 1
        vector<double> polynom = {1.0, 1.0};
        vector<complex<double>> attended_result = {complex<double>(-1.0, 0.0)};
        vector<complex<double>> test_result = eq_bairstow(polynom, 1e-15);
        BOOST_CHECK_EQUAL(test_result.size(), attended_result.size());
        for (int i = 0; i < test_result.size(); i++) {
            BOOST_REQUIRE_CLOSE(attended_result[i].real(), test_result[i].real(), 1e-15);
        }
        //Degree 2
        vector<double> quadratic = {1.0, 0.0, -1.0};
        vector<complex<double>> attended_result2 = {complex<double>(1.0, 0.0), complex<double>(-1.0, 0.0)};
        vector<complex<double>> test_result2 = eq_bairstow(quadratic, 1e-15);
        BOOST_CHECK_EQUAL(test_result2.size(), attended_result2.size());
        for (int i = 0; i < test_result.size(); i++) {
            BOOST_REQUIRE_CLOSE(attended_result2[i].real(), test_result2[i].real(), 1e-15);
        }

        //Complex Roots
        //Degree 2 x^2-2x+1
        vector<double> quadratic2 = {1.0, -2.0, 1.0};
        vector<complex<double>> attended_result_quadratic2 = {complex<double>(1.0, 0.0), complex<double>(1.0, 0.0)};
        vector<complex<double>> test_result_quadratic2 = eq_bairstow(quadratic2, 1e-15);
        BOOST_CHECK_EQUAL(test_result_quadratic2.size(), attended_result_quadratic2.size());
        for (int i = 0; i < test_result.size(); i++) {
            BOOST_REQUIRE_CLOSE(attended_result_quadratic2[i].real(), test_result_quadratic2[i].real(), 1e-15);
            BOOST_REQUIRE_CLOSE(attended_result_quadratic2[i].imag(), test_result_quadratic2[i].imag(), 1e-15);
        }

        //Degree > 2:  2x^5+x^4-2x-1 = 0
    // Expected solution : -1, 1, -1/2, i, -i
    vector<double> poly5 = {2.0, 1.0, 0.0, 0.0, -2.0, -1.0};
    vector<complex<double>> attended_result_poly5 = {complex<double>(-0.5, 0.0), complex<double>(-1.0, 0.0),
                                                    complex<double>(0.0, 1.0), complex<double>(0.0, -1.0), complex<double>(1.0, 0.0)};
    vector<complex<double>> test_result_poly5 = eq_bairstow(poly5, 1e-15);
    BOOST_CHECK_EQUAL(test_result_poly5.size(), attended_result_poly5.size());
    for (int i = 0; i < test_result_poly5.size(); i++)
    {
        BOOST_REQUIRE_SMALL(attended_result_poly5[i].real()-test_result_poly5[i].real(), 1e-15); //Change of Boost Macro method because of the null value
        BOOST_REQUIRE_SMALL(attended_result_poly5[i].imag()-test_result_poly5[i].imag(), 1e-15);
    }

//        //Degree > 2: 5x^6+0.1x^5-2x^4+x^3-30x^2+0.1x+2
//        int u = 4;
//        int v = 5;
//        vector<double> probabilities = {0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08};
//        vector<double> poly = calcul_poly(u, v, probabilities);
//        vector<complex<double>> roots = eq_bairstow(poly, 1e-15);
////        for (int i = 0; i < roots.size(); i++)
////            BOOST_TEST_MESSAGE(roots[i]);
    }
    BOOST_AUTO_TEST_CASE( creation_pii ){
        int u = 2;
        int v = 3;
        int length = 10 ;
        vector <double> probabilities = {0.2, 0.3, 0.1, 0.2, 0.1, 0.1};
        //MatrixXd m = creation_pi(probabilities, 5, -v, u);
        MatrixXd n = creation_pi_new(probabilities, 5, -v, u);
        //BOOST_TEST_MESSAGE("pi old: \n"<<m);
        BOOST_TEST_MESSAGE("pi new: \n"<<n);
        int s = 4;
        int t = 1;
        vector <double> probabilities2 = {0.2, 0.3, 0.1, 0.2, 0.1, 0.1};
        //MatrixXd o = creation_pi(probabilities2, 12, -s, t);
        MatrixXd p = creation_pi_new(probabilities2, 12, -s, t);
        //BOOST_TEST_MESSAGE("pi old: \n"<<o);
        BOOST_TEST_MESSAGE("pi new: \n"<<p);
        int y = 1;
        int x = 4;
        vector <double> probabilities3 = {0.2, 0.3, 0.1, 0.2, 0.1, 0.1};
        //MatrixXd q = creation_pi(probabilities3, 3, -x, y);
        MatrixXd r = creation_pi_new(probabilities3, 3, -x, y);
        //BOOST_TEST_MESSAGE("pi old: \n"<<q);
        BOOST_TEST_MESSAGE("pi new: \n"<<r);

    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( mcc_test )
    BOOST_AUTO_TEST_CASE( mccDeltaI ){
        int u = 5;
        int v = 5;
        int length = 10000;
        vector <double> probabilities = { 0.077670 ,0.315534 ,0.082524 ,0.000000 ,0.082524 ,0.000000 ,0.000000 ,0.077670 ,0.019417 ,0.320388 ,0.024273 };
        double p_karlin = calcul_karlin(150, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("Karlin p-value: "<<p_karlin);
        double p_mcc = calcul_mcc(150, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("mcc p-value: "<<p_mcc);
        double p_daudin = calcul_daudin(150, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value: "<<p_daudin);
    }
    // BOOST_AUTO_TEST_CASE( mccErrors ){
    //   // Test negative score expectation return error
    //   int u = 3;
    //   int v = 2;
    //   int length = 10000;
    //   vector <double> probabilities = { 0.90 ,0.02 ,0.02 ,0.02 ,0.02 ,0.02 };
    //   double p_mcc = calcul_mcc(150, probabilities, u, v, length);
    //   //BOOST_CHECK_EXCEPTION() ;
    //     BOOST_CHECK_THROW()
    // }
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( p_values_score_range )
    BOOST_AUTO_TEST_CASE( karlin_bug6732 ){
        int u = 3;
        int v = 2;
        int length = 2000;
        int localScore = 6;
        vector <double> probabilities = {0.2555,0.3935,0.2445,0.0915,0.0135,0.0015};
        double p_karlin = calcul_karlin(localScore, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("Karlin p-value: "<<p_karlin);
        double p_mcc = calcul_mcc(localScore, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("mcc p-value: "<<p_mcc);
        double p_daudin = calcul_daudin(localScore, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value: "<<p_daudin);
    }
    BOOST_AUTO_TEST_CASE( karlin_MCC_Daudin_2_7 ){
        int u = 3;
        int v = 6;
        int length = 100000;
        int localScore = 30;
        vector <double> probabilities = {0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08};
        double p_karlin = calcul_karlin(localScore, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("Karlin p-value: "<<p_karlin);
        double p_mcc = calcul_mcc(localScore, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("mcc p-value: "<<p_mcc);
        double p_daudin = calcul_daudin(localScore, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value: "<<p_daudin);
    }
    BOOST_AUTO_TEST_CASE ( karlin_MCC_Daudin_3_6){
        int u = 6;
        int v = 3;
        int length = 1000;
        int localScore = 300;
        vector <double> probabilities = {0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08}; // XXX Attention : espérance positive
        double p_karlin = calcul_karlin(localScore, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("Karlin p-value: "<<p_karlin);
        double p_mcc = calcul_mcc(localScore, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("mcc p-value: "<<p_mcc);
        double p_daudin = calcul_daudin(localScore, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value: "<<p_daudin);
    }
    BOOST_AUTO_TEST_CASE ( karlin_MCC_Daudin_4_5 ){
        int u = 5;
        int v = 4;
        int length = 100000;
        vector <double> probabilities = {0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08}; // XXX Attention : espérance positive
        double p_karlin = calcul_karlin(150, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("Karlin p-value: "<<p_karlin);
        double p_mcc = calcul_mcc(150, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("mcc p-value: "<<p_mcc);
        double p_daudin = calcul_daudin(150, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value: "<<p_daudin);
    }
    BOOST_AUTO_TEST_CASE ( karlin_MCC_Daudin_6_3 ){
        int u = 3;
        int v = 6;
        int scoreLocal = 40;
        int length = 100000;
        vector <double> probabilities = {0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08};
        double p_karlin = calcul_karlin(scoreLocal, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("Karlin p-value: "<< p_karlin);
        double p_mcc = calcul_mcc(scoreLocal, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("mcc p-value: "<< p_mcc);
        double p_daudin = calcul_daudin(scoreLocal, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value: "<<p_daudin);
    }
    BOOST_AUTO_TEST_CASE ( mcc_2_4 ){
        int u = 3;
        int v = 2;
        int length = 2000 ;
        int scoreLocal = 5;
        vector <double> probabilities = {0.2670, 0.3910, 0.2430, 0.0875, 0.0105, 0.0010};
        double p_mcc = calcul_mcc(scoreLocal, probabilities, u, v, length);
        BOOST_TEST_MESSAGE("mcc p-value: "<< p_mcc);
        double p_daudin = calcul_daudin(scoreLocal, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value: "<<p_daudin);
    }
BOOST_AUTO_TEST_CASE ( mcc_2_3 ){
  int u = 3;
  int v = 2;
  int length = 2000 ;
  int scoreLocal = 3;
  vector <double> probabilities = {0.2555,0.3935,0.2445,0.0915,0.0135,0.0015};
  double p_mcc = calcul_karlin(scoreLocal, probabilities, u, v, length);
  BOOST_TEST_MESSAGE("Karlin p-value: "<< p_mcc);
}
BOOST_AUTO_TEST_SUITE_END()


/*
BOOST_AUTO_TEST_SUITE( daudin_new )
    BOOST_AUTO_TEST_CASE( daudin_new_1 ){
        int u = 4;
        int v = 5;
        int length = 10000;
        vector <double> probabilities = {0.1, 0.05, 0.02, 0.03, 0.2, 0.3, 0.04, 0.06, 0.12, 0.08};
        double res = calcul_daudin(150, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value new: "<<res);
        double res2 = calcul_daudinOLD(150, length, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value old: "<<res2);
    }
    BOOST_AUTO_TEST_CASE( daudin_new_2 ){
        int u = 2;
        int v = 3;
        int length = 10;
        vector <double> probabilities = {0.2, 0.3, 0.1, 0.2, 0.1, 0.1};
        double res = calcul_daudin(6, 10, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value new: "<<res);

        double res2 = calcul_daudinOLD(6, 10, probabilities, -v, u);
        BOOST_TEST_MESSAGE("daudin p-value old: "<<res2);
    }
BOOST_AUTO_TEST_SUITE_END()
*/

BOOST_AUTO_TEST_SUITE( cm_exact )
    BOOST_AUTO_TEST_CASE( cm_exact_1 ){
        int u = 1;
        int v = 1;
        MatrixXd mat(3, 3);
        mat << 0.2, 0.3, 0.5,
                0.3, 0.4, 0.3,
                0.2, 0.4, 0.4;
        int length = 100;
        int localscore = 25;
        int memorysize = (localscore * (u+v+1))*(localscore * (u+v+1))*8;
        BOOST_TEST_MESSAGE("memory size: "<<memorysize/1000 << " KB");
        double res = mh_markov(localscore, mat, length, -v, u);
        BOOST_TEST_MESSAGE("cm p-value: "<<res);

    }
    BOOST_AUTO_TEST_CASE( cm_exact_2 ){
        int u = 1;
        int v = 1;
        MatrixXd mat(3, 3);
        mat << 0.2, 0.3, 0.5,
                0.3, 0.4, 0.3,
                0.2, 0.4, 0.4;
        int length = 200;
        int localscore = 100;
        int memorysize = (localscore * (u+v+1))*(localscore * (u+v+1))*8;
        BOOST_TEST_MESSAGE("memory size: "<<memorysize/1000 << " KB");
        double res = mh_markov(localscore, mat, length, -v, u);
        BOOST_TEST_MESSAGE("cm p-value: "<<res);

    }
BOOST_AUTO_TEST_SUITE_END()
