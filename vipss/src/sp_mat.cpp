#include "sp_mat.h"

SpiMat ConvertToEigenSp(const arma::sp_imat& inmat)
{
    int cols = int(inmat.n_cols);
    int rows = int(inmat.n_rows);
    SpiMat new_mat(rows, cols);
    std::vector<TripletInt> eles;
    arma::sp_imat::const_iterator it     = inmat.begin();
    arma::sp_imat::const_iterator it_end = inmat.end();
    for(; it != it_end; ++it)
    {
        eles.push_back(TripletInt(it.row(), it.col(), *it));
    }
    new_mat.setFromTriplets(eles.begin(), eles.end());
    return new_mat;
} 

SpiMat ConvertToEigenSp(const arma::sp_umat& inmat)
{
    int cols = inmat.n_cols;
    int rows = inmat.n_rows;
    SpiMat new_mat(rows, cols);
    std::vector<TripletInt> eles;
    arma::sp_umat::const_iterator it     = inmat.begin();
    arma::sp_umat::const_iterator it_end = inmat.end();
    for(; it != it_end; ++it)
    {
        eles.push_back(TripletInt(it.row(), it.col(), *it));
    }
    new_mat.setFromTriplets(eles.begin(), eles.end());
    return new_mat;
} 