#pragma once
#include "tetgen.h"
#include <string>
#include <vector>
#include "rbf_api.h"
#include <unordered_map>
#include <set>

class VoronoiGen{

    public:
        VoronoiGen() {};
        ~VoronoiGen() {};
        void loadData(const std::string& path );
        void InitMesh();
        void Tetrahedralize();
        void GetVertexStar(tetgenmesh::point p_st, std::set<tetgenmesh::point>& candid_pts, int level);
        void EstimateNormals();
        void CalcualteNormalWithVIPSS(std::vector<double>& vts, std::vector<double>& normal);
        void OrientPtNormals();
        void SavePtVn(const std::string& path);
        bool IsGoodNormal();


    public:
        tetgenio tetIO_;
        tetgenmesh tetMesh_;
        tetgenbehavior tet_beha_;
        tetgenmesh::arraypool *tetlist_, *ptlist_;
        tetgenmesh::arraypool *tetlist2_, *ptlist2_;
        tetgenmesh::arraypool *tetlist3_, *ptlist3_;
        tetgenmesh::memorypool *points_;

        RBF_API vipss_api_;
        tetgenmesh::point dummypoint_;
        std::unordered_map<tetgenmesh::point, arma::vec3> pt_normal_map_;
        std::string filename_;
        std::string out_dir_;
        std::set<tetgenmesh::point> candidate_pts_;
        

    private:
        clock_t ts_[6];
        clock_t tv_[2];


};