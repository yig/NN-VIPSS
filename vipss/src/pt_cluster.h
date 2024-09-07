#pragma once
#include <set>
#include <unordered_map>
#include "tetgen.h"
#include <armadillo>
#include "rbf_api.h"


namespace cluster
{
    
    typedef arma::vec3 vec3d;
    typedef tetgenmesh::point P3d_Ptr;
    typedef std::set<P3d_Ptr> P_Set; 
    typedef std::unordered_map<P3d_Ptr, vec3d> PtN_MAP;
    
    class PtCluster{
        public:
            PtCluster(){};
            ~PtCluster(){};
            void UpdateAdjacentClusters(); 
            void MergeCluster(PtCluster* cluster);
            std::vector<P3d_Ptr> GetClusterPts();
            void CalcualteNormals();
            std::vector<vec3d> EstimateNormals(RBF_API& rbf_api);

            void CalculateScoresWithAdjacentClusters();
            void ReplaceAdjacentClustersWithMergedClusters(PtCluster* merged_cluster);

            void OutputPtWithColor(const std::string& out_path);
            
        public: 
            typedef std::unordered_map<PtCluster*, double> ClusterScoreMap;
            P_Set key_pts_;
            P_Set group_pts_;
            std::set<PtCluster*> adjacent_clusters_;
            ClusterScoreMap adj_cluster_scores_;
            PtN_MAP pt_normal_map_;
            bool merged_ = false;
            PtCluster* merged_cluster_;
            double score_;
            size_t scored_times_ = 0;
            bool need_update_ = false;
            RBF_API* vipss_api_;

    };
    double CalculateClusterPairScore(PtCluster& a, PtCluster &b);
    double CalculateScores(std::vector<arma::vec3>& a_normals, std::vector<arma::vec3>& b_normals);


    

}  


