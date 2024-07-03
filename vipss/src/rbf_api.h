#pragma once
#include "rbfcore.h"

class RBF_API{
    public:
        RBF_API(){};
        ~RBF_API(){};
        void Set_RBF_PARA();
        void run_vipss(std::vector<double> &Vs);
        void run_vipss(std::vector<double> &Vs, size_t key_ptn);
        void run_vipss(std::vector<double> &Vs, std::vector<double> &Vn);
        void build_unit_vipss(std::vector<double> &Vs); 


    public:
        RBF_Paras para_;
        bool is_surfacing_ = true;
        bool is_outputtime_ = false;
        double user_lambda_ = 0.0; 
        int n_voxel_line_ = 80;
        std::vector<double> normals_;
        // size_t key_ptn_ = 0;
        RBF_Core rbf_core_;
        std::string outpath_ = "./vipss";
        double u_v_time = 0;

};