#pragma once
#include "kdtree.hpp"
#include "sample.h"
#include "rbf_api.h"

class IncreVipss {

    public: 
        IncreVipss(){};
        ~ IncreVipss(){};

        void Run();
        void Init();
        void InitRBFPara();
        void SampleIncrePts();

    public:
        PTSampler pt_sampler_;
        // MKdtree::MyKDTree sample_tree_;
        RBF_API vipss_api_;

        std::string input_path_ = "";
        std::string out_dir_ = "./";
        std::vector<double> in_pts_;

        double lambda_ = 0.0;
        double cur_sample_dist_ = 0.1;
        double residual_beta_ = 0.5;
        size_t init_sample_num_ = 10;
        size_t incremental_num_ = 10;
        size_t iter_num_ = 7;
        size_t n_voxel_line_ = 100;
        double add_pt_dist_scale_ = 0.1;
        double pt_dist_scale_decay_ = 0.95;

        std::vector<double> key_pts_;
        std::vector<double> auxi_pts_;
        std::vector<double> key_init_normals_;

        bool save_init_sample_pts_ = false;
        bool save_normals_ = true;
        bool is_surfacing_ = true;

};