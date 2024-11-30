#pragma once
#include <string>
#include <fstream>

struct VP_STATS {
    // time statistics for normal initializaiton
    double init_normal_total_time_ = 0; 
    double init_cluster_normal_time_ = 0;
    double cal_cluster_neigbor_scores_time_ = 0;
    double build_normal_MST_time_ = 0;
    double normal_flip_time_ = 0;

    // time statistics for optimzaiton 
    double cal_cluster_J_total_time_ = 0;
    double add_J_ele_to_triplet_vector_time_ = 0;
    double build_eigen_final_h_from_tris_time_ = 0;
    double take_h_sub_block_time_ = 0;
    double build_H_total_time_ = 0;
    double opt_solver_time_ = 0;
    int opt_func_call_num_ = 0;

    // time statistics for surface with OMP 
    // double build_per_cluster_hrbf_total_time_ = 0;
    double neighbor_search_time_ = 0;
    double cal_nn_coordinate_and_hbrf_time_ = 0;
    double surface_total_time_ = 0;
    // double generate_voroi_data_time_ = 0;
    int voxel_cal_num = 0;
    double generate_voro_data_time_ = 0;
    double build_nn_rbf_time_ = 0;
    int nn_evaluate_count_ = 0;
    double tetgen_triangulation_time_ = 0;
    int pt_num_ = 0;
    int average_neighbor_num_ = 0;
    int average_cluster_size_ = 0;
    int octree_dummy_pt_num_ = 0;
};

extern VP_STATS G_VP_stats; 
