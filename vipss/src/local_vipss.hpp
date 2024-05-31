#pragma once
#include "rbfcore.h"
#include "voronoi_gen.h"



class LocalVipss {

    typedef tetgenmesh::point P3tr;

    public:
        LocalVipss() {};
        ~ LocalVipss() {};

        void Init(const std::string & path);
        std::vector<size_t> GetClusterPtIds(size_t cluster_id);
        std::vector<size_t> GetClusterCoreIds(size_t cluster_id);
        std::vector<double> GetClusterVerticesFromIds(const std::vector<size_t>& pt_ids); 
        std::vector<double> GetClusterVertices(size_t cluster_id);

        size_t GetClusterIdFromCorePtId(const size_t pid);

        void CalculateClusterNormals(size_t cluster_id);
        void InitNormalWithVipss();
        void UpdateClusterNormals();

        double CalculateClusterPairScore(size_t c_a, size_t c_b);

        void BuildClusterAdjacentMat();
        void CalculateSingleClusterNeiScores(size_t i);

        void CalculateClusterNeiScores();
        void CalculateClusterScores();
        void MergeClusters();
        void UpdateClusterScoreMat();
        void OuputPtN(const std::string& out_path, bool orient_normal = false);
        void flipClusterNormalsByScores();

        void flipClusterNormalsByMinST();
        void GetClusterCenters();
        void BuildClusterMST();
        void FlipClusterNormalsByMST();
        void WriteVipssTimeLog();
        void SaveCluster();
        void SaveClusterPts(const std::string& path,
                            const std::vector<P3tr>& key_pts, 
                            const std::vector<P3tr>& nei_pts);
        void Run();
    
    public:
        void AppendRow(arma::sp_imat& in_mat,  arma::sp_irowvec& append_row);
        void AppendCol(arma::sp_imat& in_mat,  arma::sp_icolvec& append_row);

        double CalculateScores(std::vector<arma::vec3>& a_normals, std::vector<arma::vec3>& b_normals);

    public:
        VoronoiGen voro_gen_;
        arma::sp_imat adjacent_mat_;
        arma::sp_imat cluster_cores_mat_;
        arma::sp_imat cluster_adjacent_mat_;
        arma::sp_imat cluster_MST_mat_;

        arma::sp_mat cluster_scores_mat_;
        arma::sp_imat cluster_adjacent_pt_mat_;
        
        arma::sp_mat cluster_normal_x_;
        arma::sp_mat cluster_normal_y_;
        arma::sp_mat cluster_normal_z_;

        std::vector<tetgenmesh::point> points_; 
        std::vector<double> normals_;

        std::vector<tetgenmesh::point> cluster_centers_;

        size_t pt_num_;
        std::vector<std::pair<int, double>> cluster_id_scores_;

        RBF_API vipss_api_;

        double angle_threshold_= 30;
        size_t merged_cluster_size_ = 1;

        std::string out_dir_ = "./";
        std::string filename_ = "local_vipss"; 

        bool flip_normal_ = false;
        std::vector<double> out_pts_;
        std::vector<double> out_normals_;

        double user_lambda_ = 0.0;

        std::vector<std::pair<size_t, double>> vipss_time_stats_;
        std::vector<std::vector<std::pair<size_t, double>>> cluster_ptn_vipss_time_stats_;

};