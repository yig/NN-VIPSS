#pragma once

#include <iostream>
#include <vector>
#include <armadillo>
#include <set>

class PTSampler
{
    public:
    std::set<size_t> sample_ids;
    arma::vec min_dist_vec;
    std::vector<double> inpts;
    double cur_min_dist; 

    public:

    void init(const std::vector<double>& pts);

    double squaredDistance(double x, double y, double z, double px, double py, double pz);

    double squaredDistance(const std::vector<double>& pts, size_t i, size_t j);

    void FurthestSamplePointCloud(size_t sample_num);

    void SplitSamplePts(std::vector<double>& out_key_pts, std::vector<double>& out_auxi_pts);

    void OutputNormals(const std::vector<double>& in_normals, std::vector<double>& out_normals);

    

};


std::vector<double> CreateSpherePoints(double cx, double cy, double cz, double radius, int pt_num);



