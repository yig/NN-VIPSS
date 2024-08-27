#include "spectral_clustering.h"
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

void BuildLaplacianMat(const arma::sp_mat& adj_mat, arma::sp_mat& lap_mat)
{
    arma::sp_mat degree_mat = arma::sum(adj_mat);
    lap_mat = degree_mat - adj_mat;
}

void DeregularLaplacians(const arma::sp_mat& lap_mat)
{
    arma::eigs_opts opts;
    opts.maxiter = 10000;
    arma::vec eig_val;
    arma::mat eig_vec;
    arma::eigs_sym(eig_val, eig_vec, lap_mat, 2, "sm", opts);

}