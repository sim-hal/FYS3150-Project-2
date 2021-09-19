#include <iostream>
#include <armadillo>
#include <math.h>
#include "BeamProblem.hpp"

BeamProblem::BeamProblem(int n) {
                /*
            * n is the number of steps in the discretization
            */
            N = n - 1;
            h = 1. / n;
            a = -1. / (h * h);
            d = 2. / (h * h);

           // fill the tridiagonal 
            A.set_size(N, N);
            A.fill(0);
            A.diag(0).fill(d);
            A.diag(-1).fill(a);
            A.diag(1).fill(a);
}

void BeamProblem::fill_analytical() {
    analytical_eigenvecs.set_size(N, N);
    analytical_eigenvals.set_size(N);
    for (int i = 0; i < N; i++) {
        analytical_eigenvals(i) = d + 2 * a * cos(((i + 1) * M_PI) / (N + 1));
        for (int j = 0; j < N; j++)
            analytical_eigenvecs(i, j) = sin(((i + 1) * (j + 1) * M_PI) / (N + 1));
        }
    analytical_eigenvecs = arma::normalise(analytical_eigenvecs); 
}

void BeamProblem::compute_with_armadillo() {
    arma::eig_sym(arma_eigenvals, arma_eigenvecs, A);
}

