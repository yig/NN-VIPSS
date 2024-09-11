#include <iostream>
#include <vector>
#include <set>
#include <limits>

#include <string>
#include "sample.h"
#include "readers.h"

// namespace PTSample {


void PTSampler::init(const std::vector<double>& pts)
{
    inpts = pts;
    size_t pt_num = inpts.size() / 3;
    double MAXDV = std::numeric_limits<double>::max();
    min_dist_vec = arma::vec(pt_num, arma::fill::ones);
    min_dist_vec *= MAXDV;
    sample_ids.clear();

}

double PTSampler::squaredDistance(double x, double y, double z, double px, double py, double pz)
{
    double dx = x - px;
    double dy = y - py;
    double dz = z - pz;
    return dx * dx + dy * dy + dz * dz;
}

double PTSampler::squaredDistance(const std::vector<double>& pts, size_t i, size_t j)
{
    double dx = pts[3 * i] - pts[3 * j];
    double dy = pts[3 * i + 1] - pts[3 * j + 1];
    double dz = pts[3 * i + 2] - pts[3 * j + 2];

    return dx* dx + dy * dy + dz * dz;
}

void PTSampler::SplitSamplePts( std::vector<double>& out_key_pts, std::vector<double>& out_auxi_pts)
{
    out_key_pts.clear();
    out_auxi_pts.clear();
    size_t pt_num = inpts.size() / 3;

    for(size_t i = 0; i < pt_num; ++i)
    {
        if(sample_ids.find(i) != sample_ids.end())
        {
            out_key_pts.push_back(inpts[3* i]);
            out_key_pts.push_back(inpts[3* i + 1]);
            out_key_pts.push_back(inpts[3* i + 2]);
        } else {
            out_auxi_pts.push_back(inpts[3* i]);
            out_auxi_pts.push_back(inpts[3* i + 1]);
            out_auxi_pts.push_back(inpts[3* i + 2]);
        }
    }
}

void PTSampler::OutputNormals(const std::vector<double>& in_normals, std::vector<double>& out_normals)
{
    out_normals.clear();
    size_t pt_num = inpts.size() / 3;
    for(size_t i = 0; i < pt_num; ++i)
    {
        if(sample_ids.find(i) != sample_ids.end())
        {
            out_normals.push_back(in_normals[3*i]);
            out_normals.push_back(in_normals[3*i + 1]);
            out_normals.push_back(in_normals[3*i + 2]);
        }
    }
}

void PTSampler::FurthestSamplePointCloud(size_t sample_num)
{
    size_t pt_num = inpts.size() / 3;
    sample_num = sample_num > pt_num ?  pt_num : sample_num;
    
    size_t can_id = 0;
    if(!sample_ids.empty())
    {
        can_id = min_dist_vec.index_max();
    }
    // cout << "can id " << can_id << endl;
    while(sample_ids.size() < sample_num)
    {
        sample_ids.insert(can_id);
        cur_min_dist = min_dist_vec[can_id];
        // size_t r_id = sample_ids.size() - 1; 
        for(size_t i = 0; i < pt_num; ++i)
        {
            double dist = squaredDistance(inpts, i, can_id);
            // cout << "idst  id " << i << " dist " <<dist << endl;
            // cout << "min_dist_vec " << min_dist_vec[i] << endl;
            min_dist_vec[i] = min_dist_vec[i] < dist ? min_dist_vec[i] : dist;
        }
        can_id = min_dist_vec.index_max();
        
    }
}

// }

// https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
std::vector<double> CreateSpherePoints(double cx, double cy, double cz, double radius, int pt_num)
{
    double M_PI_ = 3.1415926535897932384626433832795028842;
    double a = 4 * M_PI_ * 1.0 * 1.0 / double(pt_num);
    double d = sqrt(a);
    int Mv = int(M_PI_ / d);
    double dv = M_PI_ / double(Mv);
    double df = a / dv;
    std::vector<double> sphere_pts;
    for(int i = 0; i < Mv; ++i)
    {
        double V = M_PI_ * (i + 0.5) / double(Mv);
        double Mf = int(2 * M_PI_ * sin(V)/df);
        for(int j =0; j < Mf; ++j)
        {
            double fi = 2 * M_PI_ * j / double(Mf);
            double x = cx + radius * sin(V) * cos(fi);
            double y = cy + radius * sin(V) * sin(fi);
            double z = cz + radius * cos(V);
            sphere_pts.push_back(x);
            sphere_pts.push_back(y);
            sphere_pts.push_back(z);
        }
    }
    return sphere_pts;   
}
