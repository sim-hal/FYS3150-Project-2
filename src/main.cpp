#include <iostream>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <cassert>
#include "project2/BeamProblem.hpp"
#include "project2/JacobiAlgorithm.hpp"

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
    const int n = 10;
    BeamProblem small_problem(n);
    small_problem.fill_analytical();
    small_problem.compute_with_armadillo();
    small_problem.compute_with_jacobi();

    assert (eigenvecs_agree(small_problem.arma_eigenvecs, small_problem.analytical_eigenvecs, small_problem.N) == 1);
    assert (eigenvecs_agree(small_problem.jacobi_eigenvecs, small_problem.analytical_eigenvecs, small_problem.N) == 1);
    assert (eigenvals_agree(small_problem.arma_eigenvals, small_problem.analytical_eigenvals) == 1);
    assert (eigenvals_agree(small_problem.jacobi_eigenvals, small_problem.analytical_eigenvals) == 1);
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
    assert (test_max_offdiag_symmetric() == 1);
    ofstream outfile("output/iters.csv");
    //J::estimate_complexity(100, outfile);
    outfile.close();
    // Problem size 10
    ofstream outfile10("output/n10.csv");
    BeamProblem problem10(10);
    problem10.compute_with_jacobi();
    problem10.fill_analytical();
    problem10.write_solutions_to_file(outfile10);
    outfile10.close();
    // Problem size 100
    ofstream outfile100("output/n100.csv");
    BeamProblem problem100(100);
    problem100.compute_with_jacobi();
    problem100.fill_analytical();
    problem100.write_solutions_to_file(outfile100);
    outfile100.close();

    return 0;
}

