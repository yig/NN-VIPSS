#pragma once
#include <pico_tree/internal/point.hpp>
#include <pico_tree/array_traits.hpp>
#include <pico_tree/kd_tree.hpp>
// Provides support for std::vector.
#include <pico_tree/vector_traits.hpp>
#include <memory>


class PicoTree {
    using data_type = double;
    using Space = std::vector<std::array<data_type, 3>>;
    using Nn = pico_tree::Neighbor<int, data_type>;
    int max_leaf_size = 8;

public:
    PicoTree(){};
    ~PicoTree(){};

void Init(const std::vector<double>& in_pts);
int SearchNearestPt(double x, double y, double z);
double NearestPtDist(double x, double y, double z);

public:
    // pico_tree::KdTree<Space> kdtree_;
    std::shared_ptr<pico_tree::KdTree<Space>>  kdtree_;
    std::vector<std::array<data_type, 3>> points_;

private:

};