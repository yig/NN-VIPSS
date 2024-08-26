#pragma once
#include <vector>
#include <string>
#include "pt_cluster.h"
#include "rbf_api.h"
#include "voronoi_gen.h"
#include <unordered_map>

class PtVipss{
        public:
            PtVipss(){};
            ~PtVipss(){};
            void Init(const std::string& path);
            void InitClusters();
            void UpdateClusters();
            // void InitClusterNormals();
            void MergeClusters();
            void OutPutNormals(const std::string& out_path);
            void UpdatePtsAndNormals();
            void Run();

        public:
            typedef cluster::PtCluster PtCluster;
            typedef cluster::P3d_Ptr P3d_Ptr;

            std::unordered_map<P3d_Ptr, PtCluster> pt_cluster_map_;
            std::vector<PtCluster*> pt_clusters_;
            std::vector<P3d_Ptr> pts_; 
            std::set<PtCluster*> cluster_set_;
            std::vector<cluster::vec3d> normals_;
            cluster::PtN_MAP all_pt_normal_map_;
            VoronoiGen voro_gen_;
            size_t pt_num_;
            RBF_API vipss_api_;
            double angle_threshold_ = 30.0;
            int merge_pairs_num_ = 1; // 
            std::string out_dir_ = "./";
            std::string out_name_ = "cluster_output";
            
    };