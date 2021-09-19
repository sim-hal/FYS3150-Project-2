#include <iostream>
#include <armadillo>
#include <math.h>
#include "BeamProblem.hpp"

using namespace std;


bool eigenvals_agree(arma::vec eigenvals1, arma::vec eigenvals2){
    const double tol = 1e-6;
    return arma::all(arma::abs(eigenvals1 - eigenvals2) < tol);
}

void small_example() {
    const int n = 7;
    BeamProblem small_problem(n);
    small_problem.fill_analytical();
    small_problem.compute_with_armadillo();
    cout << small_problem.eigenvecs_agree() << endl;

}

int main() {
    small_example();
    return 0;
}

