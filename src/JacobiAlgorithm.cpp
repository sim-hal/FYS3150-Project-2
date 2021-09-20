#include <iostream>
#include <armadillo>
#include <math.h>
#include "JacobiAlgorithm.hpp"

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
        double t, tau;
        tau = (a_ll - a_kk)/(2*a_kl);
        t = (tau > 0) ? 1.0/(tau + sqrt(1.0 + tau*tau)) : -1.0/( -tau + sqrt(1.0 + tau*tau));
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

void J::jacobi_eigensolver(arma::mat A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
const int maxiter, bool& converged, int N){
    double max_element;
    int k, l;
    int iters = 0;
    eigenvectors.eye(N , N);
    do {
        max_element = J::max_offdiag_symmetric(A, k, l, N);
        J::jacobi_rotate(A, eigenvectors, k, l, N);
        iters++;
    } while (max_element * max_element > eps && iters < maxiter);
    converged = iters < maxiter;
    eigenvalues = A.diag(0);
}