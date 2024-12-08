#pragma once
#include "rbf_api.h"
#include "readers.h"
#include "local_vipss.hpp"

namespace test_vipss_timing {

typedef std::array<double, 3> Point;
typedef std::chrono::high_resolution_clock Clock;

void cal_bbox(const std::vector<double>&pts, Point& min_corner, Point& max_corner);
void generate_test_volume_pts(const std::vector<double>& inpts, std::vector<Point>& grid_pts, int volume_res = 50);
void generate_test_volume_pts_all(const std::string& data_path, std::vector<double>& inpts, 
                                std::vector<Point>& grid_pts, int volume_res = 50);

void test_vipss_timing(const std::string& data_path );
void test_local_vipss(const std::string& in_path);
void visual_distval_pt(const std::string& in_path, int volume_res );

}



