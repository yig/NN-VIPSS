#pragma once
#include "rbfcore.h"
#include "voronoi_gen.h"
#include "local_vipss.hpp"
#include "Solver.h"

enum AXI_PlANE {
        XYZ,
        XZY,
        YZX
    };

enum HRBF_SURFACE_TYPE{
    GLOBAL_HRBF,
    LOCAL_HRBF_NN
};


// void WriteStatsLog(const std::string& path, const VP_STATS& vp_stats);


class VIPSSUnit {

    public:
        VIPSSUnit(){};
        ~VIPSSUnit(){};

        // void GetLocalVipssClusters();
        void InitPtNormalWithLocalVipss();
        // void BuildVipssUnitMatrixP();
        void OptUnitVipssNormalSimple();
        void OptUnitVipssNormalDirectSimple();
        void OptUnitVipssNormalDirect();
        void OptUnitVipssNormal();
        void ReconSurface();
        void Run();
        void test_partial_vipss();
        void BuildTetCentersMap();
        void BuildLocalHRBFPerNode();

    public:
        std::string file_name_;
        std::string data_dir_; 
        //VoronoiGen voro_gen_;
        LocalVipss local_vipss_;
        RBF_API rbf_api_;
        arma::sp_imat adjacent_mat_;
        arma::vec sg_;
        
        // std::vector<arma::mat> local_M_vec_;
        // std::vector<arma::mat> local_J_vec_;
        // std::vector<arma::sp_mat> unit_matrix_vec_;
        // double user_lambda_ = 0;
        double user_lambda_ = 0;
        arma::sp_mat Final_H_;
        size_t npt_; 
        Solution_Struct solver_;
        std::vector<double> initnormals_;
        std::vector<double> newnormals_;
        std::vector<double> s_func_vals_;
        const double M_PI_ = 3.14159265358979323846;
        bool use_hrbf_surface_ = false;
        int volume_dim_ = 100;
        size_t countopt_ = 0;
        std::vector<arma::sp_mat> temp_Hs_;
        size_t constraint_count_ = 0; 
        AXI_PlANE axi_plane_ = AXI_PlANE::XYZ;

        bool init_with_cluster_merge_ = false;
        double merge_angle_ = 40;
        static int opt_func_count;
        

        std::vector<std::shared_ptr<RBF_Core*>> node_rbf_vec;
        std::vector<double>finalMesh_v_;
        std::vector<uint>finalMesh_fv_;

        HRBF_SURFACE_TYPE hrbf_type_;

};