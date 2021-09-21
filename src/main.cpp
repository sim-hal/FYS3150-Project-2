#include <iostream>
#include <armadillo>
#include <math.h>
#include "BeamProblem.hpp"
#include "JacobiAlgorithm.hpp"

using namespace std;

bool eigenvecs_agree(const arma::mat &computed, const arma::mat &expected, int N){
    const double TOL = 1e-6;
    for (int i = 0; i < N; i++){
        arma::vec v = computed.col(i);
        arma::vec u = expected.col(i);
        if (arma::any(arma::abs(v - u) > TOL) && arma::any(arma::abs(v + u) > TOL))
            return false;
    }
    return true;
}

bool eigenvals_agree(const arma::vec &eigenvals1, const arma::vec &eigenvals2){
    const double TOL = 1e-6;
    return arma::all(arma::abs(eigenvals1 - eigenvals2) < TOL);
}

void small_example() {
    const int n = 7;
    BeamProblem small_problem(n);
    small_problem.fill_analytical();
    small_problem.compute_with_armadillo();
    small_problem.compute_with_jacobi();
    cout << eigenvecs_agree(small_problem.arma_eigenvecs, small_problem.analytical_eigenvecs, small_problem.N) << endl;
    cout << eigenvecs_agree(small_problem.jacobi_eigenvecs, small_problem.analytical_eigenvecs, small_problem.N) << endl;
    cout << eigenvals_agree(small_problem.arma_eigenvals, small_problem.analytical_eigenvals) << endl;
    cout << eigenvals_agree(small_problem.jacobi_eigenvals, small_problem.analytical_eigenvals) << endl;
    small_problem.jacobi_eigenvecs.print();
    cout << "------------" << endl;
    small_problem.analytical_eigenvecs.print();
}

bool test_max_offdiag_symmetric(){
    const double TOL = 1e-6;
    arma::mat A = {
        {1., 0., 0., .5},
        {0., 1., -.7, 0.},
        {0., -.7, 1., 0.},
        {.5, 0., 0., 1},
    };
    double expected_val = -.7;
    int expected_k = 1;
    int expected_l = 2;
    int computed_k;
    int computed_l;
    double computed_val = J::max_offdiag_symmetric(A, computed_k, computed_l, 4);
    return abs(computed_val - expected_val) < TOL && expected_k == computed_k && expected_l == computed_l;
}

int main() {
    small_example();
    cout << test_max_offdiag_symmetric() << endl;
    return 0;
}

