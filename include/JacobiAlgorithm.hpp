#ifndef JACOBIALGORITHM_H
#define JACOBIALGORITHM_H

namespace J{
    double max_offdiag_symmetric(const arma::mat& A, int& k, int &l, int N);
}

#endif