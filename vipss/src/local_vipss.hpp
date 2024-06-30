#pragma once
#include "rbfcore.h"
#include "voronoi_gen.h"



class LocalVipss {

    typedef tetgenmesh::point P3tr;

    public:
        LocalVipss() {};
        ~ LocalVipss() {};

        void Init(const std::string & path);
        inline std::vector<size_t> GetClusterPtIds(size_t cluster_id) const;
        void GetInitClusterPtIds(size_t cluster_id, 
        std::vector<double>& pts, std::vector<size_t>& pt_ids);

        inline std::vector<size_t> GetClusterCoreIds(size_t cluster_id) const;
        inline std::vector<double> GetClusterVerticesFromIds(const std::vector<size_t>& pt_ids) const; 
        std::vector<double> GetClusterVertices(size_t cluster_id) const;

        size_t GetClusterIdFromCorePtId(const size_t pid);
        void CalculateClusterNormals(size_t cluster_id);
        void InitNormalWithVipss();
        void UpdateClusterNormals();

        inline double CalculateClusterPairScore(size_t c_a, size_t c_b, bool& flip) const;
        // double CalculateClusterPairScore(size_t c_a, size_t c_b);
        void BuildClusterAdjacentMat();
        void CalculateSingleClusterNeiScores(size_t i);

        void CalculateClusterNeiScores(bool is_init = false);
        void CalculateClusterScores();
        void MergeClusters();
        void UpdateClusterScoreMat();
        void OuputPtN(const std::string& out_path, bool orient_normal = false);
        void FlipClusterNormalsByScores();

        void FlipClusterNormalsByMinST();
        void GetClusterCenters();
        void BuildClusterMST();
        void FlipClusterNormalsByMST();
        void WriteVipssTimeLog();
        void SaveClusterCorePts(const std::string& path);
        void SaveCluster();
        void SaveClusterPts(const std::string& path,
                            const std::vector<P3tr>& key_pts, 
                            const std::vector<P3tr>& nei_pts);
        
        void SaveClusterCorePts(const std::string& path,
                            const std::vector<std::vector<P3tr>>& key_pts);
        void Run();
        void InitNormals();

    public:
        inline void AppendRow(arma::sp_imat& in_mat,  arma::sp_irowvec& append_row);
        inline void AppendCol(arma::sp_imat& in_mat,  arma::sp_icolvec& append_row);
        void BuidClusterCoresPtIds();
        void UpdateClusterCoresPtIds();

        void InitSingleClusterNeiScores(size_t i);
        inline double CalculateScores(const std::vector<arma::vec3>& a_normals, const std::vector<arma::vec3>& b_normals) const;
        inline double CalculateScores(const arma::mat& a_normals, const arma::mat& b_normals) const;

        inline bool IsFlipNormal(const arma::mat& a_normals, const arma::mat& b_normals) const;
        inline bool FlipClusterNormal(size_t c_a, size_t c_b) const;

    private:

        inline void ShedCols(arma::sp_imat& in_mat, const std::vector<arma::uword>& delete_ids);
        inline void ShedCols(arma::sp_mat& in_mat, const std::vector<arma::uword>& delete_ids);
    
    public:
        VoronoiGen voro_gen_;
        arma::sp_imat adjacent_mat_;
        arma::sp_imat cluster_cores_mat_;
        arma::sp_imat cluster_adjacent_mat_;
        arma::sp_imat cluster_MST_mat_;

        arma::sp_mat cluster_scores_mat_;
        arma::sp_imat cluster_adjacent_pt_mat_;
        arma::sp_imat cluster_adjacent_flip_mat_;
        
        arma::sp_mat cluster_normal_x_;
        arma::sp_mat cluster_normal_y_;
        arma::sp_mat cluster_normal_z_;


        std::vector<tetgenmesh::point> points_; 
        std::vector<double> normals_;
        std::vector<std::vector<size_t>> cluster_core_pt_ids_;
        std::vector<tetgenmesh::point> cluster_centers_;
        size_t pt_num_;
        std::vector<std::pair<int, double>> cluster_id_scores_;
        // arma::vec cluster_scores_vec_;
        // std::vector<double> cluster_scores_vec_;
        arma::sp_colvec cluster_scores_vec_;
        std::set<size_t> update_score_cluster_ids_;
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
        const double M_PI2  = 2*acos(0.0);
        const double Anlge_PI_Rate = 180 / M_PI2;

        size_t max_iter_ = 30;
        bool use_hrbf_surface_ = false;
        int volume_dim_ = 100;
        int top_k_ = 5;
};