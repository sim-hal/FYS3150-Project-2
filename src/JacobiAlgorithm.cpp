#include <iostream>
#include <armadillo>
#include <math.h>
#include <utility>
#include <fstream>
#include "project2/JacobiAlgorithm.hpp"
#include "project2/BeamProblem.hpp"

double J::max_offdiag_symmetric(const arma::mat& A, int& k, int &l, int N) {
    k = 0;
    l = 1;
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

void J::jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int N){
    double s, c;
    double a_kl = A(k, l);
    double a_ll = A(l, l);
    double a_kk = A(k, k);
    if (a_kl != 0.){
        double tau = (a_ll - a_kk)/(2*a_kl);
        double t = (tau > 0) ? 1.0/(tau + sqrt(1.0 + tau*tau)) : -1.0/( -tau + sqrt(1.0 + tau*tau));
        c = 1/sqrt(1+t*t);
        s = c*t;
    }
    else {
        c = 1.;
        s = 0.;
    }
    A(k, k) = c*c*a_kk - 2.0*c*s*A(k, l) + s*s*a_ll;
    A(l, l) = s*s*a_kk + 2.0*c*s*A(k, l) + c*c*a_ll;
    A(k, l) = 0.0;
    A(l, k) = 0.0;

    for (int i = 0; i < N; i++){
        if (i != l && i != k) {
            double a_ik = A(i, k);
            double a_il = A(i, l);
            A(i, k) = c*a_ik - s*a_il;
            A(k, i) = A(i, k);
            A(i, l) = c*a_il + s*a_ik;
            A(l, i) = A(i, l);
        }
        double r_ik = R(i, k);
        double r_il = R(i, l);
        R(i, k) = c*r_ik - s*r_il;
        R(i, l) = c*r_il + s*r_ik;
    }
}

void sort_by_eigenvalues(arma::mat& eigenvectors, arma::vec& eigenvalues, int N){
    arma::uvec indices = arma::sort_index(eigenvalues);

    for (int i = N - 1; 0 < i; i--) {
        int new_idx = indices[i];

        while (new_idx > i) {
            new_idx = indices[new_idx];
        }

        std::swap(eigenvalues[i], eigenvalues[new_idx]);

        for (auto j = 0; j < N; j++) {
            std::swap(eigenvectors(j, i), eigenvectors(j, new_idx));
        }
    }
}

void J::jacobi_eigensolver(arma::mat A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
const int maxiter, int& iters, int N){
    double max_element;
    int k, l;
    eigenvectors.eye(N , N);
    do {
        max_element = J::max_offdiag_symmetric(A, k, l, N);
        J::jacobi_rotate(A, eigenvectors, k, l, N);
        iters++;
    } while (max_element * max_element > eps && iters < maxiter);
    eigenvalues = A.diag(0);
    sort_by_eigenvalues(eigenvectors, eigenvalues, N);
    if (maxiter < iters){
        std::cout << "Warning - did not converge" << std::endl;
    }
}

void J::estimate_complexity(int upto, std::ofstream &outfile){
    outfile << "N," << "Iterations" << std::endl;
    for (int n = 7; n < upto; n++){
        BeamProblem problem(n);
        int iters = 0;
        problem.compute_with_jacobi(iters);
        outfile << n - 1 << "," << iters << std::endl;
    }
}
