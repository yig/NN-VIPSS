#include "spectral_clustering.h"
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

void BuildLaplacianMat(const arma::sp_imat& adj_mat, arma::sp_imat& lap_mat)
{
    size_t rows = adj_mat.n_rows;
    arma::sp_irowvec degree_vec = arma::sum(adj_mat, 0);
    arma::sp_imat degree_mat = arma::diagmat(degree_vec);
    // arma::sp_imat degree_mat(rows, rows);
    // degree_mat.eye();
    // degree_mat = degree_mat * degree_vec.t();
    lap_mat = degree_mat - adj_mat;
}

void DeregularLaplacians(const arma::sp_imat& lap_mat)
{
    arma::eigs_opts opts;
    opts.maxiter = 50;
    arma::vec eig_val;
    arma::mat eig_vec;
    
    arma::sp_mat laplacian_mat = arma::conv_to<arma::sp_mat>::from(lap_mat); 
    
    arma::eigs_sym(eig_val, eig_vec, laplacian_mat, 2, "sm", opts);

    // std::cout << "eigen vals : " << eig_val << std::endl;
}

void TestSpectralClustering(const arma::sp_imat& adj_mat)
{
    arma::sp_imat lap_mat;
    auto t0 = Clock::now();
    BuildLaplacianMat( adj_mat,  lap_mat);
    auto t1 = Clock::now();
    DeregularLaplacians(lap_mat);
    auto t2 = Clock::now();

    double build_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    double eigen_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;

    
    std::cout << "------- test spectral Clustering build time : " << build_time << " -- eigen time : " << eigen_time << std::endl;
} 