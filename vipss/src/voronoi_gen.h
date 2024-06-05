#pragma once
#include "tetgen.h"
#include <string>
#include <vector>
#include "rbf_api.h"
#include <unordered_map>
#include <set>
#include "pt_cluster.h"

typedef std::set<tetgenmesh::point> P_Set; 
typedef std::unordered_map<tetgenmesh::point, arma::vec3> PtNCluster;


// struct PtCluster{
//     P_Set key_pts;
//     P_Set group_pts;
//     std::set<PtCluster*> adjacent_clusters;
//     PtNCluster pt_normal_map;
//     bool merged = false;
//     PtCluster* merged_cluster;

//     void UpdateAdjacentClusters(); 
// };

inline double PtDistance(double* p1, double* p2)
{
  return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
              (p2[1] - p1[1]) * (p2[1] - p1[1]) +
              (p2[2] - p1[2]) * (p2[2] - p1[2]));
}

class VoronoiGen{

    inline double PtDistance(double* p1, double* p2)
    {
    return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
                (p2[1] - p1[1]) * (p2[1] - p1[1]) +
                (p2[2] - p1[2]) * (p2[2] - p1[2]));
    }

    public:
        typedef cluster::PtCluster PtCluster;
        VoronoiGen() {};
        ~VoronoiGen() {};
        void loadData(const std::string& path );
        void loadData(const std::vector<double>& in_pts);
        void InitMesh();
        void Tetrahedralize();


        void BuildAdjecentMat();
        void GetVertexStar(tetgenmesh::point& p_st, std::set<tetgenmesh::point>& candid_pts, int level);
        void EstimateNormals();
        void CalcualteNormalWithVIPSS(std::vector<double>& vts, std::vector<double>& normal);
        void OrientPtNormals();
        void SavePtVn(const std::string& path, bool orient_normal = false);
        bool IsGoodNormal();
        void BuildPtIdMap();
        

        std::vector<double> ConvertPtAndNeighborsToVect(std::set<tetgenmesh::point>& candid_pts);
        void BuildPtNCluster(P_Set& pset, const std::vector<double>& normals, PtNCluster& pt_cluster);
        void CalculateScores();
        void SavePtVnColor(const std::string& path, bool orient_normal);
        void CalculatePtColors();
        void UpdateVtAndVn();
        void ScaleNormalByScore();
        void MergeCluster();

        void Run();

        std::unordered_map<tetgenmesh::point, PtCluster> BuildAllClusters(); 


    public:
        
        size_t pt_num_;
        tetgenio tetIO_;
        tetgenmesh tetMesh_;
        tetgenbehavior tet_beha_;
        tetgenmesh::arraypool *tetlist_, *ptlist_;
        tetgenmesh::arraypool *tetlist2_, *ptlist2_;
        tetgenmesh::arraypool *tetlist3_, *ptlist3_;

        std::vector<tetgenmesh::point> points_;

        RBF_API vipss_api_;
        tetgenmesh::point dummypoint_;
        
        PtNCluster pt_normal_map_;
        std::string filename_;
        std::string out_dir_;
        std::set<tetgenmesh::point> candidate_pts_;
             
        arma::sp_mat pt_score_mat_;
        arma::sp_mat pt_dist_mat_;
        arma::sp_imat pt_adjecent_mat_;

        std::unordered_map<PtCluster*, std::set<PtCluster*>> cluster_adjacent_map_;

        
        std::unordered_map<tetgenmesh::point, size_t> point_id_map_;
        std::unordered_map<tetgenmesh::point, PtCluster> point_cluster_map_;
        std::unordered_map<tetgenmesh::point, P_Set> point_cluster_pts_map_;
        std::unordered_map<tetgenmesh::point, PtNCluster> point_cluster_normal_map_;
        std::unordered_map<tetgenmesh::point, double> pt_score_map_;

        std::vector<double> vertices_;
        std::vector<double> normals_;
        std::vector<uint8_t> colors_;
        double min_angle_ = 30.0;
        

    private:
        clock_t ts_[6];
        clock_t tv_[2];


};