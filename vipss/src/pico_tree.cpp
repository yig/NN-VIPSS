#include "pico_tree.h"
#include <iostream>

void PicoTree::Init(const std::vector<double>& in_pts)
{
    points_.clear();
    for(size_t i = 0; i < in_pts.size()/3; ++i)
    {
        std::array<data_type, 3> pt{in_pts[3*i], in_pts[3*i + 1], in_pts[3*i + 2]};
        points_.push_back(pt);
    }
    kdtree_ = std::make_shared<pico_tree::KdTree<Space>>(points_, max_leaf_size);
    return;
}

int PicoTree::SearchNearestPt(double x, double y, double z)
{
    data_type query[3] = {x, y, z};
    Nn result;
    kdtree_->SearchNn(query, result);
    // double min_dist = 1.0;
    // int min_id = 0;
    // for(int i = 0; i < points_.size(); ++i)
    // {
    //     double dx = points_[i][0] - x;
    //     double dy = points_[i][1] - y;
    //     double dz = points_[i][2] - z;

    //     double dist = (dx * dx + dy * dy + dz * dz);
    //     min_dist = min_dist < dist ? min_dist : dist;
    //     min_id = min_dist < dist ? min_id : i;
    // }
    // std::cout << " caled min dist : " << min_dist << std::endl;
    // std::cout << " searched min dist : " << result.distance << std::endl;
    // std::cout << " min id " << min_id << std::endl;
    // std::cout << " search min id " << result.index << std::endl;
    int pt_id = result.index;
    return pt_id;
}



double PicoTree::NearestPtDist(double x, double y, double z)
{
    data_type query[3] = {x, y, z};
    Nn result;
    kdtree_->SearchNn(query, result);
    return result.distance;
}