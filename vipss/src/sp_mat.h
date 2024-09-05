#pragma once
#include <Eigen/Sparse>
#include <armadillo>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<int> SpiMat;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::Triplet<int> TripletInt;
typedef Eigen::SparseVector<int> SpiVec;

SpiMat ConvertToEigenSp(const arma::sp_imat& inmat);
SpiMat ConvertToEigenSp(const arma::sp_umat& inmat);