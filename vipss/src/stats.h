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
    double build_per_cluster_hrbf_total_time_ = 0;
    double neighbor_search_time_ = 0;
    double cal_nn_coordinate_and_hbrf_time_ = 0;
    double surface_total_time_ = 0;
    double generate_voroi_data_time_ = 0;
    int voxel_cal_num = 0;
};

extern VP_STATS G_VP_stats; 

// void WriteStatsLog(const std::string& path, const VP_STATS& vp_stats)
// {
//     std::ofstream log_file;
//     log_file.open(path);
//     if(log_file)
//     {
//         log_file << "time statistics for normal initializaiton : " << std::endl;

//         log_file << "cal cluster init normal with partial vipps time : " << vp_stats.init_cluster_normal_time_ << std::endl;
//         log_file << "cal cluster scores time : " << vp_stats.cal_cluster_neigbor_scores_time_ << std::endl;
//         log_file << "build normal MST tree time : " << vp_stats.build_normal_MST_time_ << std::endl;
//         log_file << "normal flip with MST tree time : " << vp_stats.normal_flip_time_ << std::endl;

//         log_file << "     " << std::endl;
//         log_file << "time statistics for optimization : " << std::endl;
//         log_file << "cal cluster J total time : " << vp_stats.cal_cluster_J_total_time_ << std::endl;
//         log_file << "add J ele to tris vector time : " << vp_stats.add_J_ele_to_triplet_vector_time_ << std::endl;
//         log_file << "get final H sub block time : " << vp_stats.take_h_sub_block_time_ << std::endl;
//         log_file << "build H total time : " << vp_stats.build_H_total_time_ << std::endl;
//         log_file << "optimization total time : " << vp_stats.build_H_total_time_ << std::endl;
//         log_file << "opt func call total num : " << vp_stats.opt_func_call_num_ << std::endl;

//         log_file << "     " << std::endl;
//         log_file << "time statistics for surface : " << std::endl;
//         log_file << "build cluster HRBF time : " << vp_stats.build_per_cluster_hrbf_total_time_ << std::endl;
//         log_file << "natural neighbor n search time : " << vp_stats.neighbor_search_time_ << std::endl;
//         log_file << "cal nn coords and HRBF time : " << vp_stats.cal_nn_coordinate_and_hbrf_time_ << std::endl;
//         log_file << "natural neighbor surfacing total time : " << vp_stats.surface_total_time_ << std::endl;
//     } 
// }