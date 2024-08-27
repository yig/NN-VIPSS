#include <armadillo>

void BuildLaplacianMat(const arma::sp_imat& adj_mat, arma::sp_mat& lap_imat);

void DeregularLaplacians(const arma::sp_imat& lap_mat);

void TestSpectralClustering(const arma::sp_imat& adj_mat);
