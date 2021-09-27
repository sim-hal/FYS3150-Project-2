#ifndef JACOBIALGORITHM_H
#define JACOBIALGORITHM_H


namespace J{
    double max_offdiag_symmetric(const arma::mat& A, int& k, int &l, int N);
    void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int N);
    void jacobi_eigensolver(const arma::mat A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iters, int N);
    void estimate_complexity(int upto, std::ofstream &outfile);
}

#endif