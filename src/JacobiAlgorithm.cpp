#include <iostream>
#include <armadillo>
#include <math.h>
#include "JacobiAlgorithm.hpp"

double J::max_offdiag_symmetric(const arma::mat& A, int& k, int &l, int N) {
    double gr = A(0, 1);
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++){
            double v = A(i, j);
            if (abs(v) > abs(gr)){
                gr = v;
                k = i;
                l = j;
            }
        }
    return gr;
}