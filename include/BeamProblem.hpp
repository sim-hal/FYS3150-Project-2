#ifndef BEAMPROBLEM_H
#define BEAMPROBLEM_H

class BeamProblem{
    private:
        int N;
        double h;
        double a;
        double d;
        const double TOL = 1e-6;
        arma::mat A;
    public:
        arma::mat computed_eigenvecs;
        arma::vec computed_eigenvals;
        arma::mat analytical_eigenvecs;
        arma::vec analytical_eigenvals;
        BeamProblem(int n);
        void fill_analytical();
        void compute_with_armadillo();
        bool eigenvecs_agree();
};

#endif