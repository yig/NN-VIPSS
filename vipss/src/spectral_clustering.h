#include <armadillo>

void BuildLaplacianMat(const arma::sp_mat& adj_mat, const arma::sp_mat& degree_mat, arma::sp_mat& lap_mat);

void DeregularLaplacians(const arma::sp_mat& lap_mat);