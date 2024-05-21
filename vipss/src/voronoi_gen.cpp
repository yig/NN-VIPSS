#include "voronoi_gen.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "readers.h"
#include "orient_normal.h"
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

REAL cps = (REAL) CLOCKS_PER_SEC;

void VoronoiGen::loadData(const std::string& path )
{

    tetIO_.load_node((char*)path.c_str());
    // tetIO_.load_ply((char*)path.c_str());
    printf("load data contain %d points! \n", tetIO_.numberofpoints);
    // tetIO_.pointlist;
}

void VoronoiGen::Tetrahedralize()
{
    tetgenmesh& m = tetMesh_;
    tetgenbehavior *b = &tet_beha_;
    tetgenio *in = &tetIO_;


    clock_t tv[13], ts[6]; // Timing informations (defined in time.h)
    

    tv[0] = clock();
    b->voroout = 1;
    m.b = b;
    m.in = in;
    // m.addin = addin;

    // if (b->metric && bgmin && (bgmin->numberofpoints > 0)) {
    //     m.bgm = new tetgenmesh(); // Create an empty background mesh.
    //     m.bgm->b = b;
    //     m.bgm->in = bgmin;
    // }

    m.initializepools();
    m.transfernodes();
    tv[1] = clock();
    if (b->refine) { // -r
        m.reconstructmesh();
    } else { // -p
        m.incrementaldelaunay(ts[0]);
    }

    tv[2] = clock();

    if (!b->quiet) {
        if (b->refine) {
        printf("Mesh reconstruction seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);
        } else {
        printf("Delaunay seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);
        if (b->verbose) {
            printf("  Point sorting seconds:  %g\n", ((REAL)(ts[0]-tv[1])) / cps);

        }
        }
    }
}

void VoronoiGen::GetVertexStar(tetgenmesh::point p_st, 
            std::set<tetgenmesh::point>& candid_pts, int level = 1)
{
    tetlist_->restart();
    ptlist_->restart();
    tetMesh_.getvertexstar(1, p_st, tetlist_,  ptlist_, NULL);
    // candid_pts.insert(p_st);

    for (size_t i = 0; i < ptlist_->objects; i++) {
        auto pt = * (tetgenmesh::point *) ptlist_-> lookup(i);
        candid_pts.insert(pt);
        if(level >= 2)
        {
            tetlist2_->restart();
            ptlist2_->restart();
            tetMesh_.getvertexstar(1, pt, tetlist2_,  ptlist2_, NULL);
            for (size_t j = 0; j < ptlist2_->objects; j+=2)
            {
                auto pt2 = * (tetgenmesh::point *) ptlist2_-> lookup(j);
                candid_pts.insert(pt2);
                if(level>=3)
                {
                    tetlist3_->restart();
                    ptlist3_->restart();
                    tetMesh_.getvertexstar(1, pt2, tetlist3_,  ptlist3_, NULL);
                    for (size_t k = 0; k < ptlist3_->objects; k++)
                    {
                        auto pt3 = * (tetgenmesh::point *) ptlist3_-> lookup(k);
                        candid_pts.insert(pt3);
                    }
                }
            }
        }
    }


    
}

void VoronoiGen::CalcualteNormalWithVIPSS(std::vector<double>& vts, std::vector<double>& normal)
{
    // RBF_API vipss_api;
    vipss_api_.Set_RBF_PARA();
    vipss_api_.run_vipss(vts);

}

std::vector<uint8_t> cal_colors(size_t ptn)
{
    vector<uint8_t> colors;
    colors.push_back(255);
    colors.push_back(0);
    colors.push_back(0);
    for (size_t i = 0; i < ptn-1; i++) {
        colors.push_back(0);
        colors.push_back(255);
        colors.push_back(0);
    }
    return colors;
}

void VoronoiGen::EstimateNormals()
{
    printf("start to estimate normals \n");
    vipss_api_.Set_RBF_PARA();
    vipss_api_.is_surfacing_ = false;
    std::vector<double> pts;
    std::vector<double> ptns;
    clock_t tv[7];
    int id = 0;
    tetMesh_.points->traversalinit();
    // printf("in ploop" );
    tetgenmesh::point ploop = tetMesh_.pointtraverse();
    // printf("in ploop : %f %f %f \n", ploop[0], ploop[1], ploop[2]);
    // ploop = tetMesh_.pointtraverse();
    tv[1] = clock();
    auto t1 = Clock::now();
    while(ploop != (tetgenmesh::point)NULL)
    {
        pts.push_back(ploop[0]);
        pts.push_back(ploop[1]);
        pts.push_back(ploop[2]);

        candidate_pts_.clear();
        GetVertexStar(ploop, candidate_pts_, 1);
        std::vector<double> vts;
        vts.push_back(ploop[0]);
        vts.push_back(ploop[1]);
        vts.push_back(ploop[2]);
        for(auto cur_pt : candidate_pts_)
        {
            if(cur_pt == ploop) continue;
            vts.push_back(cur_pt[0]);
            vts.push_back(cur_pt[1]);
            vts.push_back(cur_pt[2]);
        }
        vipss_api_.run_vipss(vts);
        arma::vec3 new_normal;
        new_normal[0] = vipss_api_.normals_[0];
        new_normal[1] = vipss_api_.normals_[1];
        new_normal[2] = vipss_api_.normals_[2];
        pt_normal_map_[ploop] = new_normal;

        ptns.push_back(vipss_api_.normals_[0]);
        ptns.push_back(vipss_api_.normals_[1]);
        ptns.push_back(vipss_api_.normals_[2]);

        ploop = tetMesh_.pointtraverse();

    }
    tv[2] = clock();
    
    auto t2 = Clock::now();
    double estimate_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
    printf("normal estimation time : %g \n", estimate_time);

    std::string out_normal_path = out_dir_ + filename_ + "_origin_ns";
    SavePtVn(out_normal_path);

    std::string out_normal_path1 = out_dir_ + filename_ + "_origin_ns1";
    writePLYFile_VN(out_normal_path1, pts, ptns);

    auto t3 = Clock::now();
    ORIENT::OrientPointNormals(pts, ptns);
    auto t4 = Clock::now();
    double orient_time = std::chrono::nanoseconds(t4- t3).count()/1e9;
    printf("Orient normal time : %g \n", orient_time);

    std::string out_normal_path2 = out_dir_ + filename_ + "_oriented_ns";
    writePLYFile_VN(out_normal_path2, pts, ptns);
    
    // tv[2] = clock();
    // OrientPtNormals();

    // tv[3] = clock();
    // printf("Orient Normal seconds:  %g\n", ((REAL)(tv[3]-tv[2])) / cps);

    // out_normal_path = out_dir_ + filename_ +  "_oriented_ns";
    // SavePtVn(out_normal_path);
    
}

void VoronoiGen::OrientPtNormals()
{
    
    tetMesh_.points->traversalinit();
    tetgenmesh::point ploop = tetMesh_.pointtraverse();
    while(ploop != (tetgenmesh::point)NULL)
    {
        if(pt_normal_map_.find(ploop) == pt_normal_map_.end())
        {
            ploop = tetMesh_.pointtraverse();
            continue;
        }
        auto p_n = pt_normal_map_[ploop];
        tetlist_->restart();
        ptlist_->restart();
        tetMesh_.getvertexstar(1, ploop, tetlist_,  ptlist_, NULL);

        int p_count = 0;
        int n_count = 0;

        for (size_t i = 0; i < ptlist_->objects; i++) {
            auto pt = * (tetgenmesh::point *) ptlist_-> lookup(i);
            if(pt_normal_map_.find(ploop) != pt_normal_map_.end())
            {
                auto cur_n = pt_normal_map_[ploop];
                if(arma::dot(cur_n, p_n) > 0) 
                {
                    p_count ++;
                } else {
                    n_count ++;
                }
            }
        }
        if(p_count < n_count)
        {
            p_n *= -1.0;
            pt_normal_map_[ploop] = p_n;
        }
        ploop = tetMesh_.pointtraverse();
    }
}

void VoronoiGen::InitMesh()
{
    Tetrahedralize();
    // dummypoint = (tetgenmesh::point) new char[3];
    // Initialize all fields of this point.
    // dummypoint_[0] = 0.0;
    // dummypoint_[1] = 0.0;
    // dummypoint_[2] = 0.0;
    
    tetlist_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::triface), 10);
    ptlist_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::point), 10);

    tetlist2_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::triface), 10);
    ptlist2_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::point), 10);

    tetlist3_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::triface), 10);
    ptlist3_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::point), 10);
}




void VoronoiGen::SavePtVn(const std::string& path)
{
    std::vector<double> pts;
    std::vector<double> ptns;
    for(auto &ele : pt_normal_map_)
    {
        // printf("point : %f %f %f \n", ele.first[0], ele.first[1], ele.first[2]);
        // printf("normal : %f %f %f \n", ele.second[0], ele.second[1], ele.second[2]);
        if(ele.first == (tetgenmesh::point)NULL) continue;
        pts.push_back(ele.first[0]);
        pts.push_back(ele.first[1]);
        pts.push_back(ele.first[2]);

        ptns.push_back(ele.second[0]);
        ptns.push_back(ele.second[1]);
        ptns.push_back(ele.second[2]);
        
    }
    // std::string out_normal_path = "data/vipss_estimated_ns";
    writePLYFile_VN(path, pts, ptns);
    printf("successfully save points and normals to file : %s", path.c_str());
}