#pragma once
#include <unordered_set>
#include "rbfcore.h"
#include "voronoi_gen.h"
#include "sp_mat.h"


class C_Edege{
    public:
        C_Edege(size_t c_a, size_t c_b)
        {
            c_a_ = c_a;
            c_b_ = c_b; 
            // if(c_a > c_b)
            // {
            //     c_a_ = c_b;
            //     c_b_ = c_a;
            // }
            e_id_ = std::to_string(c_a_) + "_" + std::to_string(c_b_);
        };
        ~C_Edege () {};

    public:
        size_t c_a_;
        size_t c_b_;
        double score_;
        std::string e_id_;
};

struct SP_BBOX{
    double min_corner[3];
    double max_corner[3];
};

class LocalVipss {

    typedef tetgenmesh::point P3tr;
    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::SparseMatrix<int> SpiMat;
    typedef Eigen::Triplet<double> Triplet;
    

    public:
        LocalVipss() {};
        ~ LocalVipss() {};

        void Init(const std::string & path);
        inline std::vector<size_t> GetClusterPtIds(size_t cluster_id) const;
        void GetInitClusterPtIds(size_t cluster_id, 
        std::vector<double>& pts, std::vector<size_t>& pt_ids);
        void BuildHRBFPerNode();
        inline std::vector<double> GetClusterNormalsFromIds
                            (const std::vector<size_t>& pt_ids, const std::vector<double>& all_normals) const;

        inline std::vector<double> GetClusterSvalsFromIds(const std::vector<size_t>& pt_ids, 
                        const std::vector<double>& all_svals) const; 


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
        void MergeGoodNeighborClusters();
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
        void InitNormalsWithMerge();

        void GroupPtsWithVolume();

    public:
        inline void AppendRow(arma::sp_umat& in_mat,  arma::sp_urowvec& append_row);
        inline void AppendCol(arma::sp_umat& in_mat,  arma::sp_ucolvec& append_row);
        void BuidClusterCoresPtIds();
        void UpdateClusterCoresPtIds();

        void InitSingleClusterNeiScores(size_t i);
        inline double CalculateScores(const std::vector<arma::vec3>& a_normals, const std::vector<arma::vec3>& b_normals) const;
        inline double CalculateScores(const arma::mat& a_normals, const arma::mat& b_normals) const;

        inline bool IsFlipNormal(const arma::mat& a_normals, const arma::mat& b_normals) const;
        inline bool FlipClusterNormal(size_t c_a, size_t c_b) const;
        void BuildCurrentH(const arma::sp_mat& unit_cluster_mat, const arma::sp_mat&J_m, double lambda);
        inline void AddClusterHMatrix(const std::vector<size_t>& p_ids, const arma::mat& J_m,size_t npt);
        void CalHVecSum();
        void BuildMatrixH();

        void SampleClusterPts();
        void TestInsertPt();
        void TestVoronoiPts();
        // double DistanceFunction(const std::vector<tetgenmesh::point>& nn_pts);
        double NodeDistanceFunction(const tetgenmesh::point nn_pt, const tetgenmesh::point cur_pt);
        double NatureNeighborDistanceFunction(const tetgenmesh::point cur_pt);
        double NatureNeighborDistanceFunctionOMP(const tetgenmesh::point cur_pt);
        static double NNDistFunction(const R3Pt &in_pt);  
        void SetThis();      
        void VisualFuncValues(double (*function)(const R3Pt &in_pt), const VoroPlane& plane,
                              const std::string& dist_fuc_color_path);
        
        void testNNPtDist();

        static int DistCallNum;
        static double DistCallTime;

        // void CalClusterDegrees();
        // void MergeHRBFClustersWithDegree();
        // void MergeClusterPairs(std::vector<arma::uword>& merged_cluster_ids);
        void GroupClustersWithDegree();

        void SaveGroupPtsWithColor(const std::string& path);
        void PtPCA(std::vector<double>& pts);
        void OptimizeAdjacentMat();

        void InitAdjacentData();
        void MergeHRBFClustersWithMap();
        void BuildHRBFPerCluster();

    private:

        inline void ShedCols(arma::sp_umat& in_mat, const std::vector<arma::uword>& delete_ids);
        inline void ShedCols(arma::sp_mat& in_mat, const std::vector<arma::uword>& delete_ids);
    
    public:
        bool is_group_cluster_ = false;
        VoronoiGen voro_gen_;
        RBF_API vipss_api_;

        SpiMat cluster_adjacent_emat_;
        SpiMat cluster_adjacent_pt_emat_;
        SpiMat cluster_cores_emat_;
        SpiMat valid_pt_diag_emat_;

        std::vector<std::unordered_set<size_t>> cluster_adjacent_ids_;
        std::vector<std::unordered_set<size_t>> cluster_pt_ids_;
        std::vector<std::vector<size_t> > cluster_core_pt_ids_vec_; 
        arma::uvec cluster_core_pt_nums_;
        arma::vec cluster_volume_vals_;
        std::vector<size_t> cluster_id_map_;
        std::vector<size_t> valid_cluster_dist_map_;

        arma::sp_umat adjacent_mat_;
        arma::sp_umat cluster_cores_mat_;
        arma::sp_umat cluster_valid_cores_mat_;
        arma::sp_umat cluster_adjacent_mat_;
        
        arma::sp_umat cluster_adjacent_mat_opt_;
        arma::sp_umat cluster_MST_mat_;
        arma::sp_mat cluster_scores_mat_;
        arma::sp_umat cluster_adjacent_pt_mat_;
        arma::sp_umat cluster_valid_adjacent_pt_mat_;
        arma::sp_umat cluster_adjacent_share_pt_mat_;
        arma::sp_umat cluster_adjacent_flip_mat_;
        arma::sp_mat cluster_normal_x_;
        arma::sp_mat cluster_normal_y_;
        arma::sp_mat cluster_normal_z_;
        arma::sp_mat final_H_; 
        arma::sp_mat final_H_temp_;
        arma::vec cluster_degrees_; 
        arma::sp_umat valid_pt_diag_mat_;
        arma::uvec cluster_valid_sign_vec_;

        std::vector<tetgenmesh::point> points_; 
        std::vector<double> normals_;
        std::vector<double> s_vals_;
        std::vector<std::vector<size_t>> cluster_core_pt_ids_;
        std::vector<tetgenmesh::point> cluster_centers_;
        size_t pt_num_;
        std::vector<std::pair<int, double>> cluster_id_scores_;
        arma::sp_colvec cluster_scores_vec_;
        std::set<size_t> update_score_cluster_ids_;
        std::vector<double> out_pts_;
        std::vector<double> out_normals_;
        std::vector<std::pair<size_t, double>> vipss_time_stats_;
        std::vector<std::vector<std::pair<size_t, double>>> cluster_ptn_vipss_time_stats_;

        std::vector<arma::sp_mat> cluster_pt_mat_vec_; 
        std::vector<arma::mat> cluster_J_mat_vec_;
        std::vector<arma::sp_mat> temp_H_vec_;
        const size_t temp_H_max_num_ = 2048;
        // std::vector<std::vector<size_t>> cluster_pt_ids_vec_;
        std::vector< std::unordered_map<int,double>> temp_Hmat_values_;
        SpMat final_h_eigen_;
        // std::unordered_map<long ,double> h_pos_value_map_;
        std::vector<Triplet> h_ele_triplets_;
        std::set<P3tr> sample_cluster_pts_;
        std::vector<std::shared_ptr<RBF_Core>> node_rbf_vec_;
        std::vector<double>finalMesh_v_;
        std::vector<uint>finalMesh_fv_;
        arma::vec nn_dist_vec_;
        arma::vec nn_volume_vec_;
        
    public:
        int max_group_iter_ = 12;
        int max_group_pt_num_ = 1024;
        bool flip_normal_ = false;
        double user_lambda_ = 0.0;
        double unit_lambda_ = 0.0;
        const double M_PI2  = 2*acos(0.0);
        const double Anlge_PI_Rate = 180 / M_PI2;
        size_t max_iter_ = 30;
        bool use_hrbf_surface_ = false;
        int volume_dim_ = 100;
        int top_k_ = 5;
        double angle_threshold_= 45;
        size_t merged_cluster_size_ = 1;
        std::string out_dir_ = "./";
        std::string filename_ = "local_vipss"; 
        double build_j_time_total_ = 0;
        double build_h_time_total_ = 0;
        double build_h_time_total2_ = 0;
        std::vector<tetgenmesh::point> insert_pts_;

        double pass_time_sum_ = 0;
        double search_nn_time_sum_ = 0;
        double dist_time_sum_ = 0;
        size_t dist_call_num_ = 0;
        size_t in_cluster_surface_pt_count = 0;
};