#pragma once
#include "tetgen.h"
#include <string>
#include <vector>
#include "rbf_api.h"
#include <unordered_map>
#include <set>
#include "pt_cluster.h"
#include "pico_tree.h"
// #include <eigen3/Eigen/Dense>

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

struct VoroPlane{
    double px, py, pz;
    double nx, ny, nz;
    double intPx, intPy, intPz;
    VoroPlane(){};

    VoroPlane(double in_px, double in_py, double in_pz, double in_nx, double in_ny, double in_nz)
    {
        px = in_px; py = in_py; pz = in_pz;
        nx = in_nx; ny = in_ny; nz = in_nz; 
    };
};




inline double PtDistance(double* p1, double* p2)
{
  return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
              (p2[1] - p1[1]) * (p2[1] - p1[1]) +
              (p2[2] - p1[2]) * (p2[2] - p1[2]));
}

class VoronoiGen{

    enum EdgeSign {
        IN_UNION,
        INTERSECT,
        OUT_UNION,
        IN_VALID
    };

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
        void GenerateVoroData();

        void Run();
        void InsertPt(tetgenmesh::point pt);
        std::unordered_map<tetgenmesh::point, PtCluster> BuildAllClusters(); 

        void GetVoroCellEdgeList(tetgenmesh::point nei_pt, std::vector<tetgenio::voroedge*>& edge_list);
        void GetVoroCellPtIdList(tetgenmesh::point nei_pt, std::set<int>& vcell_pt_list);
        void GetVoroCellPtAndEdgeIdList(tetgenmesh::point nei_pt, std::set<int>& vcell_pt_list, std::set<int>& vcell_edge_list);

        void CalVoroCellPtSign(const VoroPlane &plane, const std::set<int>& vcell_pt_list, std::vector<int>& positive_pids);
        // void CalVoroEdgeIntersection(tetgenmesh::point nei_pt, const VoroPlane& plane, std::vector<double>& intersect_pts);
        void CalVoroEdgeIntersection(const std::set<int>& cell_eids, std::vector<double>& intersect_pts);
        void GetVoroCellPtList(std::vector<tetgenio::voroedge*>& edge_list, std::vector<double>& pt_list);
        std::vector<double> GetPolygonSplitPts(tetgenmesh::point cur_pt, tetgenmesh::point nei_pt);
        bool CalPlaneEdgeIntersection(tetgenio::voroedge* edge, VoroPlane &plane);
        double CalSplitPolyhedronVolume(std::vector<double>& split_pts);
        double CalUnionCellVolume(tetgenmesh::point in_pt, tetgenmesh::point nei_pt);
        // double CalVoroFaceVolume(tetgenio::facet* facet,  )
        void CalVoroFaceNormal(tetgenmesh::point in_pt, tetgenmesh::point nei_pt);
        void GetTruncatedCellPts(tetgenmesh::point in_pt, tetgenmesh::point nei_pt);
        void InitVoronoiDataForCellVolume();
        double CalTruncatedCellVolume(tetgenmesh::point in_pt, tetgenmesh::point nei_pt);

        void OutputVoronisMesh();
        void InsertBoundryPts();
        void GetVoronoiNeiPts(tetgenmesh::point pt, std::vector<tetgenmesh::point>& candid_pts);
        void BuildTetMeshTetCenterMap();
        void BuildPicoTree();
        tetgenmesh::tetrahedron* GetClosetTet(double x, double y, double z);
        double CalTruncatedCellVolumePass(tetgenmesh::point in_pt, tetgenmesh::point nei_pt);
                
    public:
        size_t pt_num_;
        tetgenio tetIO_;
        tetgenio tetInsertIO_;
        tetgenmesh tetMesh_;
        tetgenbehavior tet_beha_;
        tetgenmesh::arraypool *tetlist_, *ptlist_;
        tetgenmesh::arraypool *tetlist2_, *ptlist2_;
        tetgenmesh::arraypool *tetlist3_, *ptlist3_;
        std::vector<double> tet_center_pts_;
        std::unordered_map<int, tetgenmesh::tetrahedron*> tc_pt_tet_map_;
        std::vector<tetgenmesh::point> points_;
        // std::unordered_map<int, double> vcell_pt_sign_values_;
        std::vector<double> vpt_sign_vals_;
        std::vector<int> edge_insect_symbols_;
        std::vector<double> edge_insect_pts_;
        std::vector<arma::vec3> vcell_face_centers_; 
        std::vector<arma::vec3> vcell_face_normals_; 
        double truncated_cell_center_[3];
        double truncated_face_center_[3];

        static const int MAX_ELEMENT_ = 10000;
        double pt_sign_vals_[MAX_ELEMENT_];
        int edge_sign_vals_[MAX_ELEMENT_];
        int face_sign_vals_[MAX_ELEMENT_];
        int edge_intersect_pt_ids_[MAX_ELEMENT_];
        int edge_positive_pt_ids_[MAX_ELEMENT_];
        int face_insersect_pt_ids_[MAX_ELEMENT_];
        double face_centers_[MAX_ELEMENT_];
        double intersect_pts_[MAX_ELEMENT_];
        double tet_list_pts_[MAX_ELEMENT_ *3];

        RBF_API vipss_api_;
        tetgenmesh::point dummypoint_;
        tetgenio voronoi_data_;
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
        PicoTree pTree_;
        
    private:
        clock_t ts_[6];
        clock_t tv_[2];


};