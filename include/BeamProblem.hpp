#ifndef BEAMPROBLEM_H
#define BEAMPROBLEM_H

class BeamProblem{
    private:
        double h;
        double a;
        double d;
        const double TOL = 1e-6;
        arma::mat A;
    public:
        int N;
        arma::mat arma_eigenvecs;
        arma::vec arma_eigenvals;
        arma::mat analytical_eigenvecs;
        arma::vec analytical_eigenvals;
        BeamProblem(int n);
        void fill_analytical();
        void compute_with_armadillo();
};

#endif