#pragma once
#include "rbfcore.h"
#include <memory>

class RBF_API{
    public:
        RBF_API(){};
        ~RBF_API(){};
        void Set_RBF_PARA();
        void run_vipss(std::vector<double> &Vs);
        void run_vipss(std::vector<double> &Vs, size_t key_ptn);
        void run_vipss_for_incremental(std::vector<double> &Vs, size_t key_ptn);
        void run_vipss(std::vector<double> &Vs, std::vector<double> &Vn);
        void run_vipss(std::vector<double> &Vs, std::vector<double> &Vn, const std::vector<double>& s_vals);
        void build_unit_vipss(std::vector<double> &Vs); 
        void build_unit_vipss(std::vector<double> &Vs, size_t key_npt); 
        void BuildLocalVipssVec();

        void build_cluster_hrbf(std::vector<double> &Vs, std::vector<double> &Vn, 
                                 const std::vector<double>& s_vals, std::shared_ptr<RBF_Core> rbf_ptr);

        void build_cluster_hrbf_surface(std::shared_ptr<RBF_Core> rbf_ptr, const std::string& mesh_path);

    public:

        RBF_Paras para_;
        bool is_surfacing_ = false;
        bool is_outputtime_ = false;
        double user_lambda_ = 0.0; 
        double unit_lambda_ = 0.0; 
        int n_voxel_line_ = 80;
        std::vector<double> normals_;
        // size_t key_ptn_ = 0;
        RBF_Core rbf_core_;
        std::string outpath_ = "./vipss";
        double u_v_time = 0;
        std::vector<size_t> p_ids_;
        // std::vector<double> dist_vals_; 
        arma::vec auxi_dist_vec_; 
        std::vector<double> key_opt_normals_;
        std::vector<double> pre_out_normals_;
        bool use_input_normal_ = false;
        double incre_vipss_beta_ = 1.0;

        std::vector<std::shared_ptr<RBF_Core>> local_rbf_vec;

};

