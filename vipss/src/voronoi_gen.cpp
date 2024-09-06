#include "voronoi_gen.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "readers.h"
#include "orient_normal.h"
#include <chrono>
#include <math.h>
#include <cmath>
// #include <libqhull/libqhull.h>
// #include <libqhullcpp/Qhull.h>
// #include <libqhullcpp/RboxPoints.h>
// #include <libqhull/qhull_a.h>
#include "sample.h"
#include "omp.h"

typedef std::chrono::high_resolution_clock Clock;

REAL cps = (REAL) CLOCKS_PER_SEC;
// double M_PI  = 2*acos(0.0);

double CalVolumeQhull(std::vector<double>& pts)
{
    // orgQhull::Coordinates newPts(pts);
    // orgQhull::RboxPoints inpts;
    // inpts.setDimension(3);
    // inpts.append(pts); 

    // const char* cmd = "Fx";
    // orgQhull::Qhull q_convex(inpts, cmd);
    double volume = 0;
    // double volume = q_convex.volume();

    // printf("qhull volume : %Lf \n", volume);
    // std::cout << " qhull volume " << volume << std::endl;
    return volume;
}


// void PtCluster::UpdateAdjacentClusters()
// {
//     std::set<PtCluster*> new_clusters;
//     for(auto &cluster : adjacent_clusters)
//     {
//         if(cluster->merged)
//         {
//             if(cluster->merged_cluster != this)
//             {
//                 new_clusters.insert(cluster->merged_cluster);
//             }
//         } else {
//             new_clusters.insert(cluster);
//         }
//     }
// }

void VoronoiGen::loadData(const std::string& path )
{

    tetIO_.load_node((char*)path.c_str());
    // tetIO_.load_ply((char*)path.c_str());
    printf("load data contain %d points! \n", tetIO_.numberofpoints);
    // tetIO_.pointlist;
}

void VoronoiGen::loadData(const std::vector<double>& in_pts)
{
    tetIO_.load_node((double*)in_pts.data(), in_pts.size()/3);
}

void VoronoiGen::InsertSphereBoundryPts()
{
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();

    double max_x = std::numeric_limits<double>::min();
    double max_y = std::numeric_limits<double>::min();
    double max_z = std::numeric_limits<double>::min();

    auto& in_pts = tetIO_.pointlist;
    for(size_t i =0; i < tetIO_.numberofpoints; ++i)
    {
        min_x = min_x < in_pts[3*i]     ? min_x : in_pts[3*i];
        min_y = min_y < in_pts[3*i + 1] ? min_y : in_pts[3*i + 1];
        min_z = min_z < in_pts[3*i + 2] ? min_z : in_pts[3*i + 2];

        max_x = max_x > in_pts[3*i]     ? max_x : in_pts[3*i];
        max_y = max_y > in_pts[3*i + 1] ? max_y : in_pts[3*i + 1];
        max_z = max_z > in_pts[3*i + 2] ? max_z : in_pts[3*i + 2];
    }
    double cx = (min_x + max_x) / 2.0;
    double cy = (min_y + max_y) / 2.0;
    double cz = (min_z + max_z) / 2.0;
    
    double scale = 2.0;
    double dx = (max_x - min_x) / 2.0 * scale;
    double dy = (max_y - min_y) / 2.0 * scale;
    double dz = (max_z - min_z) / 2.0 * scale;
    double radius = std::max(dx, std::max(dy, dz)); 
    int pt_num = 32;
    auto sphere_pts = CreateSpherePoints(cx, cy, cz, radius, pt_num);
    pt_num = sphere_pts.size()/3;
    printf("sphere pt size : %ld \n",sphere_pts.size()/3 );
    std::string out_sphere_path = out_dir_ + "sphere_boundry.xyz";
    writeXYZ(out_sphere_path, sphere_pts);
    printf("successfully write sphere boundary pts to file : %s\n", out_sphere_path.c_str());
    // printf(out_sphere_path);
    std::vector<tetgenmesh::point> boundary_pts;

    for(size_t i = 0; i < pt_num; ++i)
    {
        tetgenmesh::point newpt;
        auto& temp_mesh = tetMesh_;
        temp_mesh.makepoint(&newpt, tetgenmesh::UNUSEDVERTEX);
        newpt[0] = sphere_pts[i * 3];
        newpt[1] = sphere_pts[i * 3 + 1];
        newpt[2] = sphere_pts[i * 3 + 2];
        tetgenmesh::insertvertexflags ivf;
        ivf.bowywat = 1; // Use Bowyer-Watson algorithm
        ivf.lawson = 0;
        tetgenmesh::triface searchtet;
        searchtet.tet = NULL;
        ivf.iloc = (int) tetgenmesh::UNKNOWN;
        auto t0 = Clock::now();
        tetgenmesh::face *splitsh = NULL;
        tetgenmesh::face *splitseg = NULL;
        // temp_mesh.setpointtype(newpt, tetgenmesh::UNUSEDVERTEX);
        if(temp_mesh.insertpoint(newpt, &searchtet, splitsh, splitseg,  &ivf))
        {
            // printf(" %ld insertion succeeded ! \n", i);
        }
        boundary_pts.push_back(newpt);
        // temp_mesh.setpointtype(newpt, tetgenmesh::UNUSEDVERTEX);
        // printf(" insertion pt type : %d! \n", temp_mesh.pointtype(newpt));
        // printf("insert pt id : %ld \n", i);
    }

    for(auto pt : boundary_pts)
    {
        tetMesh_.setpointtype(pt, tetgenmesh::UNUSEDVERTEX);
    }

}

void VoronoiGen::InsertBoundryPts()
{
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();

    double max_x = std::numeric_limits<double>::min();
    double max_y = std::numeric_limits<double>::min();
    double max_z = std::numeric_limits<double>::min();

    auto& in_pts = tetIO_.pointlist;
    for(size_t i =0; i < tetIO_.numberofpoints; ++i)
    {
        min_x = min_x < in_pts[3*i]     ? min_x : in_pts[3*i];
        min_y = min_y < in_pts[3*i + 1] ? min_y : in_pts[3*i + 1];
        min_z = min_z < in_pts[3*i + 2] ? min_z : in_pts[3*i + 2];

        max_x = max_x > in_pts[3*i]     ? max_x : in_pts[3*i];
        max_y = max_y > in_pts[3*i + 1] ? max_y : in_pts[3*i + 1];
        max_z = max_z > in_pts[3*i + 2] ? max_z : in_pts[3*i + 2];
    }
    // bbox final_scale = scale + 1 
    double scale = 2.0;
    double dx = (max_x - min_x) / 2.0 * scale;
    double dy = (max_y - min_y) / 2.0 * scale;
    double dz = (max_z - min_z) / 2.0 * scale;

    min_x -= dx; min_y -= dy; min_z -= dz;
    max_x += dx; max_y += dy; max_z += dz;

    //     5-------6
    //    /|      /|
    //   4-|-----7 |
    //   | 1-----|-2
    //   |/      |/
    //   0-------3 
    double box_pts[14][3];
    box_pts[0][0] = min_x; box_pts[0][1] = min_y; box_pts[0][2] = min_z;
    box_pts[1][0] = min_x; box_pts[1][1] = max_y; box_pts[1][2] = min_z;
    box_pts[2][0] = max_x; box_pts[2][1] = max_y; box_pts[2][2] = min_z;
    box_pts[3][0] = max_x; box_pts[3][1] = min_y; box_pts[3][2] = min_z;

    box_pts[4][0] = min_x; box_pts[4][1] = min_y; box_pts[4][2] = max_z;
    box_pts[5][0] = min_x; box_pts[5][1] = max_y; box_pts[5][2] = max_z;
    box_pts[6][0] = max_x; box_pts[6][1] = max_y; box_pts[6][2] = max_z;
    box_pts[7][0] = max_x; box_pts[7][1] = min_y; box_pts[7][2] = max_z;

    // face centers 
    // f0 : 0 1 2 3
    // f1 : 0 3 7 4
    // f2 : 3 2 6 7
    // f3 : 2 1 5 6
    // f4 : 0 1 5 4
    // f5 : 4 5 6 7
if(1)
{
    box_pts[8][0] = (box_pts[0][0] + box_pts[1][0] + box_pts[2][0] + box_pts[3][0])/4.0;
    box_pts[8][1] = (box_pts[0][1] + box_pts[1][1] + box_pts[2][1] + box_pts[3][1])/4.0;
    box_pts[8][2] = (box_pts[0][2] + box_pts[1][2] + box_pts[2][2] + box_pts[3][2])/4.0;

    box_pts[9][0] = (box_pts[0][0] + box_pts[3][0] + box_pts[7][0] + box_pts[4][0])/4.0;
    box_pts[9][1] = (box_pts[0][1] + box_pts[3][1] + box_pts[7][1] + box_pts[4][1])/4.0;
    box_pts[9][2] = (box_pts[0][2] + box_pts[3][2] + box_pts[7][2] + box_pts[4][2])/4.0;

    box_pts[10][0] = (box_pts[3][0] + box_pts[2][0] + box_pts[6][0] + box_pts[7][0])/4.0;
    box_pts[10][1] = (box_pts[3][1] + box_pts[2][1] + box_pts[6][1] + box_pts[7][1])/4.0;
    box_pts[10][2] = (box_pts[3][2] + box_pts[2][2] + box_pts[6][2] + box_pts[7][2])/4.0;

    box_pts[11][0] = (box_pts[2][0] + box_pts[1][0] + box_pts[5][0] + box_pts[6][0])/4.0;
    box_pts[11][1] = (box_pts[2][1] + box_pts[1][1] + box_pts[5][1] + box_pts[6][1])/4.0;
    box_pts[11][2] = (box_pts[2][2] + box_pts[1][2] + box_pts[5][2] + box_pts[6][2])/4.0;
    
    box_pts[12][0] = (box_pts[0][0] + box_pts[1][0] + box_pts[5][0] + box_pts[4][0])/4.0;
    box_pts[12][1] = (box_pts[0][1] + box_pts[1][1] + box_pts[5][1] + box_pts[4][1])/4.0;
    box_pts[12][2] = (box_pts[0][2] + box_pts[1][2] + box_pts[5][2] + box_pts[4][2])/4.0;

    box_pts[13][0] = (box_pts[4][0] + box_pts[5][0] + box_pts[7][0] + box_pts[6][0])/4.0;
    box_pts[13][1] = (box_pts[4][1] + box_pts[5][1] + box_pts[7][1] + box_pts[6][1])/4.0;
    box_pts[13][2] = (box_pts[4][2] + box_pts[5][2] + box_pts[7][2] + box_pts[6][2])/4.0;
}

    for(size_t i = 0; i < 14; ++i)
    {
        tetgenmesh::point newpt;
        auto& temp_mesh = tetMesh_;
        temp_mesh.makepoint(&newpt, tetgenmesh::UNUSEDVERTEX);
        newpt[0] = box_pts[i][0];
        newpt[1] = box_pts[i][1];
        newpt[2] = box_pts[i][2];
        tetgenmesh::insertvertexflags ivf;
        ivf.bowywat = 1; // Use Bowyer-Watson algorithm
        ivf.lawson = 0;
        tetgenmesh::triface searchtet;
        searchtet.tet = NULL;
        ivf.iloc = (int) tetgenmesh::UNKNOWN;
        auto t0 = Clock::now();
        tetgenmesh::face *splitsh = NULL;
        tetgenmesh::face *splitseg = NULL;
        // temp_mesh.setpointtype(newpt, tetgenmesh::UNUSEDVERTEX);
        if(temp_mesh.insertpoint(newpt, &searchtet, splitsh, splitseg,  &ivf))
        {
            // printf(" %ld insertion succeeded ! \n", i);
        }
        temp_mesh.setpointtype(newpt, tetgenmesh::UNUSEDVERTEX);
        // printf(" insertion pt type : %d! \n", temp_mesh.pointtype(newpt));
    }
}

void VoronoiGen::Tetrahedralize()
{
    tetgenmesh& m = tetMesh_;
    tetgenbehavior *b = &tet_beha_;
    tetgenio *in = &tetIO_;

    clock_t tv[13], ts[6]; // Timing informations (defined in time.h)
    printf("start to Tetrahedralize \n");
    tv[0] = clock();
    b->voroout = 0;
    // b->refine = 1;
    // b->plc = 1;
    m.b = b;
    m.in = in;
    // m.addin = addin;
    // addin != NULL
    m.initializepools();
    m.transfernodes();
    tv[1] = clock();
    m.incrementaldelaunay(ts[0]);

    BuildPtIdMap();
    printf("finsh BuildPtIdMap \n");
    // InsertBoundryPts();
    // GenerateVoroData();
    // m.facetverticeslist
    // m.meshsurface();
    // m.recoverboundary(ts[0]);
    // m.reconstructmesh();
    // m.outvoronoi(&voronoi_data_);
    // voronoi_data_.vpointlist
    // m.generate_voronoi_cell(&voronoi_data_);
    // printf("voronoi pt num : %d \n", voronoi_data_.numberofvpoints);
    // printf("voronoi edge num : %d \n", voronoi_data_.numberofvedges);
    // printf("voronoi facet num : %d \n", voronoi_data_.numberofvfacets);
    // printf("voronoi cell num : %d \n", voronoi_data_.numberofvcells);

    // printf("finsh generate_voronoi_cell \n");

    // InitVoronoiDataForCellVolume();
    
    // OutputVoronisMesh();
    
    // printf("finsh to Tetrahedralize \n");
}

void VoronoiGen::GenerateVoroData()
{
    // InsertBoundryPts();
// 
    auto t0 = Clock::now();
    InsertSphereBoundryPts();
    

    // printf("finsh InsertBoundryPts \n");
    tetMesh_.generate_voronoi_cell(&voronoi_data_);
    

    printf("voronoi pt num : %d \n", voronoi_data_.numberofvpoints);
    // printf("voronoi edge num : %d \n", voronoi_data_.numberofvedges);
    // printf("voronoi facet num : %d \n", voronoi_data_.numberofvfacets);
    printf("voronoi cell num : %d \n", voronoi_data_.numberofvcells);
    // printf("finsh generate_voronoi_cell \n");
    // InitVoronoiDataForCellVolume();
    // return;
    printf("finsh InitVoronoiDataForCellVolume \n");
    auto t1 = Clock::now();
    // PrecomputeVoroData();
    auto t2 = Clock::now();
    double build_vo_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("-----------compute voronoi data time : %f \n", build_vo_time);

    // double prevo_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
    // printf("-----------pre iter voronoi data time : %f \n", prevo_time);

    if(0)
    {
        OutputVoronisMesh();
        std::string out_mesh_path = out_dir_ + "tet_mesh.obj";
        OutputTetMesh(out_mesh_path);
    }
    // printf("finsh to Tetrahedralize \n");
}

void VoronoiGen::GetVoroCellEdgeList(tetgenmesh::point nei_pt, std::vector<tetgenio::voroedge*>& edge_list)
{
    // printf("pt size 0000 : %ld \n", points_.size());
    auto vc_id   = point_id_map_[nei_pt];
    auto v_cell  = voronoi_data_.vcelllist[vc_id];
    // printf("pt size 000000 : %ld \n", points_.size());
    auto vf_list = voronoi_data_.vfacetlist;
    auto ve_list = voronoi_data_.vedgelist;
    // printf("pt size 0001 : %ld \n", points_.size());
    if(v_cell == NULL) return;
    int f_num = v_cell[0];
    // printf("f_num : %ld \n", f_num);
    std::set<tetgenio::voroedge*> cell_edges;
    
    for(size_t i = 0; i < f_num; ++i)
    {
        size_t f_id = v_cell[i + 1];
        // printf("v_cell f id 0002 : %ld \n", f_id);
        auto facet = vf_list[f_id];
        // printf("pt size 0002 : %ld \n", points_.size());
        size_t e_num = facet.elist[0];
        for(size_t j = 0; j < e_num; ++j)
        {
            auto& ve = ve_list[facet.elist[1 + j]];
            cell_edges.insert(&ve);
        }
    }
    edge_list = std::vector<tetgenio::voroedge*>(cell_edges.begin(), cell_edges.end());
}

void VoronoiGen::GetVoroCellPtList(std::vector<tetgenio::voroedge*>& edge_list, std::vector<double>& pt_list)
{
    auto& vplist = voronoi_data_.vpointlist;
    std::set<int> pidset;
    for(auto edge : edge_list)
    {
        // auto pt1 = vplist[edge->v1];
        pidset.insert(edge->v1);
        pidset.insert(edge->v2);
    }
    for(auto pid : pidset)
    {
        pt_list.push_back(vplist[3 * pid]);
        pt_list.push_back(vplist[3 * pid + 1]);
        pt_list.push_back(vplist[3 * pid + 2]);
    }
}

void VoronoiGen::GetVoroCellPtAndEdgeIdList(tetgenmesh::point nei_pt, std::set<int>& vcell_pt_list, std::set<int>& vcell_edge_list)
{
    auto vc_id   = point_id_map_[nei_pt];
    auto v_cell  = voronoi_data_.vcelllist[vc_id];
    auto vf_list = voronoi_data_.vfacetlist;
    auto ve_list = voronoi_data_.vedgelist;
    // printf("pt size 0001 : %ld \n", points_.size());
    if(v_cell == NULL) return;
    int f_num = v_cell[0];
    std::set<int> p_ids;
    for(size_t i = 0; i < f_num; ++i)
    {
        size_t f_id = v_cell[i + 1];
        auto facet = vf_list[f_id];
        size_t e_num = facet.elist[0];
        arma::vec3 f_center;
        for(size_t j = 0; j < e_num; ++j)
        {
            int e_id = facet.elist[1 + j];
            vcell_edge_list.insert(e_id);
            const auto& ve = ve_list[e_id];
            vcell_pt_list.insert(ve.v1);
            vcell_pt_list.insert(ve.v2);
        }
    }
}

void VoronoiGen::CalVoroEdgeIntersection(const std::set<int>& cell_eids, std::vector<double>& intersect_pts)
{

    auto ve_list = voronoi_data_.vedgelist;
    auto vp_list = voronoi_data_.vpointlist;
    // printf("pt size 0001 : %ld \n", points_.size());
    // printf("start to cal facet center 000\n");
    // printf("vpt_sign_vals_ size %ld \n", vpt_sign_vals_.size());

    for(auto e_id : cell_eids)
    {
        auto& ve = ve_list[e_id];
        double t1 = vpt_sign_vals_[ve.v1];
        double t2 = vpt_sign_vals_[ve.v2];
        if(t1 * t2 <= 0)
        {
            double sum = abs(t1) + abs(t2);
            if(sum > 0) 
            {
                t1 = abs(t1) / sum;
                t2 = abs(t2) / sum;
            }
            double p1_x = vp_list[3 * ve.v1]; 
            double p1_y = vp_list[3 * ve.v1 + 1]; 
            double p1_z = vp_list[3 * ve.v1 + 2];
            double p2_x = vp_list[3 * ve.v2]; 
            double p2_y = vp_list[3 * ve.v2 + 1]; 
            double p2_z = vp_list[3 * ve.v2 + 2];
            double ix = t2 * p1_x + t1 * p2_x;
            double iy = t2 * p1_y + t1 * p2_y;
            double iz = t2 * p1_z + t1 * p2_z;
            intersect_pts.push_back(ix);
            intersect_pts.push_back(iy);
            intersect_pts.push_back(iz);
        }
    }
}

void VoronoiGen::CalVoroFaceNormal(tetgenmesh::point in_pt, tetgenmesh::point nei_pt)
{
    auto vc_id   = point_id_map_[nei_pt];
    auto v_cell  = voronoi_data_.vcelllist[vc_id];
    auto vf_list = voronoi_data_.vfacetlist;
    auto ve_list = voronoi_data_.vedgelist;
    auto vp_list = voronoi_data_.vpointlist;

    if(v_cell == NULL) return;

    int f_num = v_cell[0];
    std::set<int> p_ids;
    for(size_t i = 0; i < f_num; ++i)
    {
        size_t f_id = v_cell[i + 1];
        auto facet = vf_list[f_id];
        size_t e_num = facet.elist[0];
        if(e_num >= 2)
        {
            arma::vec3 e_n1;
            arma::vec3 e_n2;
            int v1 = ve_list[facet.elist[1 + 0]].v1;
            int v2 = ve_list[facet.elist[1 + 0]].v2;
            e_n1[0] = vp_list[3 * v1]     - vp_list[3 * v2];
            e_n1[1] = vp_list[3 * v1 + 1] - vp_list[3 * v2 + 1];
            e_n1[2] = vp_list[3 * v1 + 2] - vp_list[3 * v2 + 2];
            v1 = ve_list[facet.elist[1 + 1]].v1;
            v2 = ve_list[facet.elist[1 + 1]].v2;
            e_n2[0] = vp_list[3 * v1]     - vp_list[3 * v2];
            e_n2[1] = vp_list[3 * v1 + 1] - vp_list[3 * v2 + 1];
            e_n2[2] = vp_list[3 * v1 + 2] - vp_list[3 * v2 + 2];

            arma::vec3 normal = arma::cross(e_n1, e_n2);
            normal = arma::normalise(normal);
        }
        
    }
}

void VoronoiGen::InitVoronoiDataForCellVolume()
{
    // size_t vp_num = voronoi_data_.numberofvpoints;
    // size_t ve_num = voronoi_data_.numberofvedges;
    // size_t vf_num = voronoi_data_.numberofvfacets;

    // vpt_sign_vals_.resize(vp_num);
    // edge_insect_symbols_.resize(ve_num, -1);
    // edge_insect_pts_.resize(ve_num * 3);
    // vcell_face_normals_.resize(vf_num);
    // vcell_face_centers_.resize(vf_num);
    // auto vp_list = voronoi_data_.vpointlist;
    // auto ve_list = voronoi_data_.vedgelist;
    // auto vf_list = voronoi_data_.vfacetlist;

    // int f_num = voronoi_data_.numberofvfacets;
    // for(size_t i = 0; i < f_num; ++i)
    // {
    //     auto facet = vf_list[i];
    //     size_t e_num = facet.elist[0];
    //     arma::vec3 f_center;
    //     std::set<int> facet_p_ids;
    //     for(size_t j = 0; j < e_num; ++j)
    //     {
    //         auto& ve = ve_list[facet.elist[1 + j]];
    //         facet_p_ids.insert(ve.v1);
    //         facet_p_ids.insert(ve.v2);
    //     }
    //     for(auto id : facet_p_ids)
    //     {
    //         f_center[0] += vp_list[3 *id];
    //         f_center[1] += vp_list[3 *id + 1];
    //         f_center[2] += vp_list[3 *id + 2];
    //     }
    //     if(!facet_p_ids.empty())
    //     {
    //         f_center /= double(facet_p_ids.size());
    //     }
    //     vcell_face_centers_[i] = f_center;
    // }
}

void VoronoiGen::CalVoroCellPtSign(const VoroPlane &plane, 
                                   const std::set<int>& vcell_pt_list, 
                                   std::vector<int>& positive_pids)
{
    // printf("vpt_sign_vals_ size : %ld \n", vpt_sign_vals_.size());
    // vcell_pt_sign_values_.clear();
    auto vp_list = voronoi_data_.vpointlist;
    for(auto p_id : vcell_pt_list)
    {
        double dx = vp_list[3 * p_id]     - plane.px;
        double dy = vp_list[3 * p_id + 1] - plane.py;
        double dz = vp_list[3 * p_id + 2] - plane.pz;
        double proj = dx * plane.nx + dy * plane.ny + dz * plane.nz;
        vpt_sign_vals_[p_id] = proj; 
        if(proj > 0)
        {
            positive_pids.push_back(p_id);
        }   
    } 
}


bool VoronoiGen::CalPlaneEdgeIntersection(tetgenio::voroedge* edge, VoroPlane &plane)
{
    auto vp_list = voronoi_data_.vpointlist;
    auto ax = vp_list[3 * edge->v1];
    auto ay = vp_list[3 * edge->v1 + 1];
    auto az = vp_list[3 * edge->v1 + 2];

    auto bx = vp_list[3 * edge->v2];
    auto by = vp_list[3 * edge->v2 + 1];
    auto bz = vp_list[3 * edge->v2 + 2];

    double ab_x = bx - ax;
    double ab_y = by - ay;
    double ab_z = bz - az;

    double ap_x =  ax - plane.px;
    double ap_y =  ay - plane.py;
    double ap_z =  az - plane.pz;

    double ap_dot = ap_x * plane.nx + ap_y * plane.ny + ap_z * plane.nz;
    double abp = ab_x * plane.nx  + ab_y * plane.ny + ab_z * plane.nz;
    if(abs(abp) < 1e-16) return false;
    double t = - ap_dot / abp;
    if(t <0 || t > 1.0) return false;

    plane.intPx = ax + ab_x * t;
    plane.intPy = ay + ab_y * t;
    plane.intPz = az + ab_z * t;
    return true;
}

std::vector<double> VoronoiGen::GetPolygonSplitPts(tetgenmesh::point cur_pt, tetgenmesh::point nei_pt)
{
    std::vector<tetgenio::voroedge*> vc_edges;
    GetVoroCellEdgeList(nei_pt, vc_edges);

    std::vector<double> vc_pts;
    GetVoroCellPtList(vc_edges, vc_pts);

    double mid_x = (cur_pt[0] + nei_pt[0])/2;
    double mid_y = (cur_pt[1] + nei_pt[1])/2;
    double mid_z = (cur_pt[2] + nei_pt[2])/2;

    double p_nx = cur_pt[0] - nei_pt[0];
    double p_ny = cur_pt[1] - nei_pt[1];
    double p_nz = cur_pt[2] - nei_pt[2];
    double n_len = sqrt(p_nx * p_nx + p_ny * p_ny + p_nz * p_nz);
    // n_len = 1e-14 > n_len ? 1e-14 : n_len; 
    p_nx /= n_len;
    p_ny /= n_len;
    p_nz /= n_len;

    std::vector<double> split_pts;
    for(size_t i = 0; i < vc_pts.size()/3; ++i)
    {
        double vp_nx = vc_pts[3 * i]     - mid_x;
        double vp_ny = vc_pts[3 * i + 1] - mid_y;
        double vp_nz = vc_pts[3 * i + 2] - mid_z;
        double sign = vp_nx * p_nx + vp_ny * p_ny + vp_nz * p_nz;
        if(sign > 0)
        {
            split_pts.push_back(vc_pts[3 * i]);
            split_pts.push_back(vc_pts[3 * i + 1]);
            split_pts.push_back(vc_pts[3 * i + 2]);
        }
    }
    VoroPlane v_plane(mid_x, mid_y, mid_z, p_nx, p_ny, p_nz);
    
    for(auto edge : vc_edges)
    {
        if(CalPlaneEdgeIntersection(edge, v_plane))
        {
            // double dx = v_plane.intPx - mid_x;
            // double dy = v_plane.intPy - mid_y;
            // double dz = v_plane.intPz - mid_z;
            // double res = v_plane.nx * dx + v_plane.ny * dy + v_plane.nz * dz;
            // printf(" res -------------- %f \n", res);
            split_pts.push_back(v_plane.intPx);
            split_pts.push_back(v_plane.intPy);
            split_pts.push_back(v_plane.intPz);
            // myfile_split << "v " << v_plane.intPx << " " << v_plane.intPy << " " << v_plane.intPz << std::endl;
            // v_id++;
        }
    }
    return split_pts;
}

inline double CalPtLen(double ax, double ay, double az, double bx, double by, double bz)
{
    double dx = ax - bx;
    double dy = ay - by;
    double dz = az - bz;
    double len = sqrt(dx * dx + dy * dy + dz * dz);
    return len;
}

double CalTetrahedronVolume(double *pa, double* pb, double* pc, double* pd)
{
    double len_ab = CalPtLen(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2]);
    double len_bc = CalPtLen(pb[0], pb[1], pb[2], pc[0], pc[1], pc[2]);
    double len_ac = CalPtLen(pa[0], pa[1], pa[2], pc[0], pc[1], pc[2]);

    double s = (len_ab + len_ac + len_bc) / 2;
    double area = sqrt(s * (s - len_ab) *(s - len_ac) *(s- len_bc));

    double ab_x = pa[0] - pb[0];
    double ab_y = pa[1] - pb[1];
    double ab_z = pa[2] - pb[2];

    double bc_x = pb[0] - pc[0];
    double bc_y = pb[1] - pc[1];
    double bc_z = pb[2] - pc[2];

    double nx = ab_y * bc_z - ab_z * bc_y;
    double ny = ab_z * bc_x - ab_x * bc_z;
    double nz = ab_x * bc_y - ab_y * bc_x;

    double n_len = sqrt(nx * nx + ny * ny + nz * nz);

    if(n_len < 1e-16) return 0;
    nx /= n_len;
    ny /= n_len;
    nz /= n_len;

    double ad_x = pd[0] - pa[0];
    double ad_y = pd[1] - pa[1];
    double ad_z = pd[2] - pa[2];
    double proj = nx * ad_x + ny * ad_y + nz * ad_z;
    double volume = abs(proj) * area / 3.0;
    return volume;
}


inline double CalTetrahedronVolumeDet(const double* pa, const double* pb, const double* pc, const double* pd)
{
    double a00 = pa[0] - pd[0];
    double a10 = pb[0] - pd[0];
    double a20 = pc[0] - pd[0];

    double a01 = pa[1] - pd[1];
    double a11 = pb[1] - pd[1];
    double a21 = pc[1] - pd[1];

    double a02 = pa[2] - pd[2];
    double a12 = pb[2] - pd[2];
    double a22 = pc[2] - pd[2];

    double det1 = a00 * a11 * a22 + a01 * a12 * a20 + a02 * a10 * a21;
    double det2 = a02 * a11 * a20 + a01 * a10 * a22 + a00 * a12 * a21;
    double det = abs(det2 - det1) / 6.0;
    return det;
}

// inline double CalTetrahedronVolumeDet(double* pa, double* pb, double* pc, double* pd)
// {
//     double a00 = pa[0] - pd[0];
//     double a10 = pb[0] - pd[0];
//     double a20 = pc[0] - pd[0];

//     double a01 = pa[1] - pd[1];
//     double a11 = pb[1] - pd[1];
//     double a21 = pc[1] - pd[1];

//     double a02 = pa[2] - pd[2];
//     double a12 = pb[2] - pd[2];
//     double a22 = pc[2] - pd[2];

//     double det1 = a00 * a11 * a22 + a01 * a12 * a20 + a02 * a10 * a21;
//     double det2 = a02 * a11 * a20 + a01 * a10 * a22 + a00 * a12 * a21;
//     double det = abs(det2 - det1) / 6.0;
//     return det;
// }

double CalTetrahedronVolume(tetgenmesh::tetrahedron* searchtet)
{
    // int j = (searchtet->ver & 3); // The current vertex index.
    auto pa = (tetgenmesh::point) searchtet[4 + 0];
    auto pb = (tetgenmesh::point) searchtet[4 + 1];
    auto pc = (tetgenmesh::point) searchtet[4 + 2];
    auto pd = (tetgenmesh::point) searchtet[4 + 3];

    double volume = CalTetrahedronVolumeDet(pa, pb, pc, pd);
    // double volume = CalTetrahedronVolume(pa, pb, pc, pd);
    return volume;
}



double VoronoiGen::CalSplitPolyhedronVolume(std::vector<double>& split_pts)
{
    // printf("............. CalSplitPolyhedronVolume \n");
    tetgenio temp_io;
    temp_io.load_node((double*)split_pts.data(), split_pts.size()/3);
    tetgenmesh temp_mesh;
    temp_mesh.in = &temp_io;
    tet_beha_.quiet = 1;
    temp_mesh.b = &tet_beha_;
    temp_mesh.b->epsilon = 1e-18;
    // printf("............. CalSplitPolyhedronVolume 000\n");
    temp_mesh.initializepools();
    // printf("............. CalSplitPolyhedronVolume 0000000\n");
    try{
        temp_mesh.transfernodes();
    } catch(int e)
    {
        std::cout << "error code " << e << std::endl;
        return 0;
    }
    
    clock_t t_val;
    // printf("............. CalSplitPolyhedronVolume 001\n");
    double volume_sum = 0;
    try{
        temp_mesh.incrementaldelaunay(t_val);
    }
    catch(int e)
    {
        std::cout << "error code " << e << std::endl;
        return 0;
    }
    temp_mesh.tetrahedrons->traversalinit();

    
    // printf("......tet num :  %d \n", temp_mesh.tetrahedrons->items);
    // printf(" out path : %s \n", out_path);
    tetgenmesh::triface tetloop;
    tetloop.tet = (tetgenmesh::tetrahedron *) temp_mesh.alltetrahedrontraverse();
    int tid = 0;
    while(tetloop.tet != (tetgenmesh::tetrahedron *) NULL)
    {
        // auto pa = temp_mesh.org(tetloop);
        // auto pb = temp_mesh.dest(tetloop);
        // auto pc = temp_mesh.apex(tetloop);
        // auto pd = temp_mesh.oppo(tetloop);

        if (temp_mesh.ishulltet(tetloop)) 
        {
            tetloop.tet = (tetgenmesh::tetrahedron *) temp_mesh.alltetrahedrontraverse();
            continue; 
        }
        // std::string out_path = "../out/tet_" + std::to_string(tid) + ".obj";
        // ofstream myfile;
        // myfile.open(out_path);
        // for(size_t i = 0; i < 4; ++i)
        // {
        //     auto cur_pt = (tetgenmesh::point) tetloop.tet[4 + i];
        //     // printf("pt type : %d \n", temp_mesh.pointtype(cur_pt));
        //     myfile << "v " << pa[0] << " " << pa[1] << " " << pa[2] << std::endl;
        //     myfile << "v " << pb[0] << " " << pb[1] << " " << pb[2] << std::endl;
        //     myfile << "v " << pc[0] << " " << pc[1] << " " << pc[2] << std::endl;
        //     myfile << "v " << pd[0] << " " << pd[1] << " " << pd[2] << std::endl;
        // }
        // myfile << "f 1 2 3" << std::endl;
        // myfile << "f 1 2 4" << std::endl;
        // myfile << "f 1 3 4" << std::endl;
        // myfile << "f 2 3 4" << std::endl;
        // myfile.close();         
        // myfile << ploop[4]
        double volume = CalTetrahedronVolume(tetloop.tet);
        volume_sum += volume;
        tid++;
        tetloop.tet = (tetgenmesh::tetrahedron *) temp_mesh.alltetrahedrontraverse();
        // if(temp_mesh. )
    }
    // printf("-------------- tid : %d \n", tid);
    
    // std::string out_path1 = "../out/tet_pts.obj";
    // ofstream myfile;
    // myfile.open(out_path1);
    // temp_mesh.points->traversalinit();
    // auto pt = (tetgenmesh::point)temp_mesh.points->traverse();
    // while(pt != NULL)
    // {
    //     myfile << "v " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
    //     pt = (tetgenmesh::point)temp_mesh.points->traverse();
    // }
    return volume_sum;
}

double VoronoiGen::CalUnionCellVolume(tetgenmesh::point in_pt, tetgenmesh::point nei_pt)
{
    // std::vector<tetgenio::voroedge*> split_polygon_edges;
    // GetVoroCellEdgeList(nei_pt, split_polygon_edges);
    auto split_pts = GetPolygonSplitPts(in_pt, nei_pt);
    // std::ofstream out_file;
    // std:;string out_path = "../out/split_pts.obj";
    // out_file.open(out_path);
    // for(size_t i = 0; i < split_pts.size()/3; ++i)
    // {
    //     out_file <<"v "<< split_pts[3*i] <<" " << split_pts[3 *i + 1] << " " << split_pts[3 * i + 2] << std::endl;
    // }
    // out_file.close();
    

    double volume = CalSplitPolyhedronVolume(split_pts);
    // double volume = CalVolumeQhull(split_pts);
    return volume;
}

double VoronoiGen::CalTruncatedCellVolume(tetgenmesh::point in_pt, tetgenmesh::point nei_pt)
{
    double mid_x = (in_pt[0] + nei_pt[0])/2.0;
    double mid_y = (in_pt[1] + nei_pt[1])/2.0;
    double mid_z = (in_pt[2] + nei_pt[2])/2.0;
    double p_nx  = (in_pt[0] - nei_pt[0]);
    double p_ny  = (in_pt[1] - nei_pt[1]);
    double p_nz  = (in_pt[2] - nei_pt[2]);
    double pn_len = sqrt(p_nx * p_nx + p_ny * p_ny + p_nz * p_nz);
    if(pn_len > 0) 
    {
        p_nx /= pn_len;
        p_ny /= pn_len;
        p_nz /= pn_len;
    }

    // printf("start to retrive voronois cell pt and edges\n");

    VoroPlane v_plane(mid_x, mid_y, mid_z, p_nx, p_ny, p_nz);
    std::set<int> vcell_pids;
    std::set<int> vcell_eids;
    GetVoroCellPtAndEdgeIdList(nei_pt, vcell_pids, vcell_eids);

    // printf("start to cal CalVoroCellPtSign \n");
    std::vector<int> truncated_cell_pt_list;
    CalVoroCellPtSign(v_plane, vcell_pids, truncated_cell_pt_list);
    // printf("truncated cell pt size : %ld \n", truncated_cell_pt_list.size());

    std::vector<double> truncated_cell_pts;
    auto vp_list = voronoi_data_.vpointlist;
    for(auto p_id : truncated_cell_pt_list)
    {
        truncated_cell_pts.push_back(vp_list[3 * p_id]);
        truncated_cell_pts.push_back(vp_list[3 * p_id + 1]);
        truncated_cell_pts.push_back(vp_list[3 * p_id + 2]);
    }
    CalVoroEdgeIntersection(vcell_eids, truncated_cell_pts);

    // printf("truncated_cell_pts size : %ld \n", truncated_cell_pts.size()/3);
    // double volume = 0;
    // try{
    // double volume = CalVolumeQhull(truncated_cell_pts);
    // auto t000 = Clock::now();
    if(truncated_cell_pts.size() < 4) return 0;
    
    double volume = CalSplitPolyhedronVolume(truncated_cell_pts);
    // double volume = CalVolumeQhull(truncated_cell_pts);
    // auto t001 = Clock::now();
    // double cal_volume_time = std::chrono::nanoseconds(t001 - t000).count()/1e9;
    // printf(" cal_volume_time with tet: %g \n", cal_volume_time);
    // }
    // catch(const std::runtime_error& e)
    // {
    //     std::cout << e.what() << std::endl;
    //     return volume;
    // }
    // std::cout << " tet volume " << volume << std::endl;
    // printf("tet volume : %Lf \n", volume);
    // catch(const std::exception& e)  
    // double volume = CalSplitPolyhedronVolume(truncated_cell_pts);
    return volume;
}

// int tet_vol_count = 0;

void VoroPlane::SavePlane(const std::string& outpath)
{
    arma::vec3 ax_1 = {1, 1, 0}; 
    ax_1[2] = -(nx * ax_1[0] + ny * ax_1[1]) / nz;
    double len = sqrt(ax_1[0] * ax_1[0] + ax_1[1] * ax_1[1] + ax_1[2] * ax_1[2]);
    ax_1[0] /= len;
    ax_1[1] /= len;
    ax_1[2] /= len;
    arma::vec3 pv{nx, ny, nz};
    arma::vec3 ax_2 = arma::cross(ax_1, pv);
    std::ofstream out_file;
    out_file.open(outpath);
    arma::vec3 ori{px, py, pz};
    arma::vec3 pa = ori + ax_1 * 2;
    arma::vec3 pb = ori + ax_2 * 2;
    arma::vec3 pc = ori - ax_1 * 2;
    arma::vec3 pd = ori - ax_2 * 2;

    out_file << "v " << pa[0] << " " << pa[1] << " " << pa[2] << std::endl;
    out_file << "v " << pb[0] << " " << pb[1] << " " << pb[2] << std::endl;
    out_file << "v " << pc[0] << " " << pc[1] << " " << pc[2] << std::endl;
    out_file << "v " << pd[0] << " " << pd[1] << " " << pd[2] << std::endl;

    out_file << "f 1 2 3" << std::endl;
    out_file << "f 1 4 3" << std::endl;
    out_file.close();
}

double VoronoiGen::CalTruncatedCellVolumePass(tetgenmesh::point in_pt, tetgenmesh::point nei_pt)
{
    double plane_mid_x = (in_pt[0] + nei_pt[0])/2.0;
    double plane_mid_y = (in_pt[1] + nei_pt[1])/2.0;
    double plane_mid_z = (in_pt[2] + nei_pt[2])/2.0;
    double p_nx  = (in_pt[0] - nei_pt[0]);
    double p_ny  = (in_pt[1] - nei_pt[1]);
    double p_nz  = (in_pt[2] - nei_pt[2]);
    double pn_len = sqrt(p_nx * p_nx + p_ny * p_ny + p_nz * p_nz);
    if(pn_len > 0) 
    {
        p_nx /= pn_len;
        p_ny /= pn_len;
        p_nz /= pn_len;
    }
    // VoroPlane v_plane(plane_mid_x, plane_mid_y, plane_mid_z, p_nx, p_ny, p_nz);
    // SavePlane(v_plane);
    if(point_id_map_.find(nei_pt) == point_id_map_.end())
    return 0;
    auto vc_id   = point_id_map_[nei_pt];
    auto v_cell  = voronoi_data_.vcelllist[vc_id];
    auto vf_list = voronoi_data_.vfacetlist;
    auto ve_list = voronoi_data_.vedgelist;
    auto vp_list = voronoi_data_.vpointlist;
    // printf("pt size 0001 : %ld \n", points_.size());
    if(v_cell == NULL) return 0;
    int f_num = v_cell[0];
    std::set<int> p_ids;
    int e_count = 0;
    double cell_center[3] = {0, 0, 0};
    int valid_count = 0;
    int intersect_count = 0; 
    double trunc_facet_center_[3] = {0, 0, 0};
    // std::set<int> ptid_set;
    for(size_t i = 0; i < f_num; ++i)
    {
        // printf("face id : %ld \n", i);
        size_t f_id = v_cell[i + 1];
        const auto& facet = vf_list[f_id];
        size_t e_num = facet.elist[0];
        arma::vec3 f_center;

        double face_center[3] = {0, 0, 0};
        int positive_count = 0;
        int p_count = 0;
        int face_edge_inter_count = 0;
        int valid_face_count = 0;
        for(size_t j = 0; j < e_num; ++j)
        {
            // printf("edge id : %ld \n", j);
            // edge_insect_symbols_[e_count] = -1;
            size_t e_id = facet.elist[1 + j];
            // printf("e_id %ld \n", e_id);
            const auto& ve = ve_list[e_id];
            if(ve.v1 != -1 && ve.v2 != -1)
            {
                double dx = vp_list[3 * ve.v1]     - plane_mid_x;
                double dy = vp_list[3 * ve.v1 + 1] - plane_mid_y;
                double dz = vp_list[3 * ve.v1 + 2] - plane_mid_z;
                double proj1 = dx * p_nx + dy * p_ny + dz * p_nz;

                dx = vp_list[3 * ve.v2]     - plane_mid_x;
                dy = vp_list[3 * ve.v2 + 1] - plane_mid_y;
                dz = vp_list[3 * ve.v2 + 2] - plane_mid_z;
                double proj2 = dx * p_nx + dy * p_ny + dz * p_nz;
                // printf("proj1 %f, proj2 %f \n", proj1, proj2);
                if(proj1 * proj2 < 0 )
                {
                    double r1 = abs(proj1) /(abs(proj1) + abs(proj2));
                    double r2 = abs(proj2) /(abs(proj1) + abs(proj2));

                    double inter_x = r2 * vp_list[3 * ve.v1]     + r1 * vp_list[3 * ve.v2];
                    double inter_y = r2 * vp_list[3 * ve.v1 + 1] + r1 * vp_list[3 * ve.v2 + 1];
                    double inter_z = r2 * vp_list[3 * ve.v1 + 2] + r1 * vp_list[3 * ve.v2 + 2];
                    
                    intersect_pts_[3 * intersect_count]     = inter_x;
                    intersect_pts_[3 * intersect_count + 1] = inter_y;
                    intersect_pts_[3 * intersect_count + 2] = inter_z;
                    // edge_insect_pts_[3*e_count]     = r2 * vp_list[3 * ve.v1]     + r1 * vp_list[3 * ve.v2];
                    // edge_insect_pts_[3*e_count + 1] = r2 * vp_list[3 * ve.v1 + 1] + r1 * vp_list[3 * ve.v2 + 1];
                    // edge_insect_pts_[3*e_count + 2] = r2 * vp_list[3 * ve.v1 + 2] + r1 * vp_list[3 * ve.v2 + 2];
                    int inter_id = 2 * i + face_edge_inter_count;
                    face_insersect_pt_ids_[inter_id] = intersect_count;
                    edge_intersect_pt_ids_[e_count]  = intersect_count;
                    edge_sign_vals_[e_count] = INTERSECT;
                    // edge_intersect_pt_ids_[e_count] 
                    // face_insersect_pt_ids_[2 * i + intersect_count] = e_count;
                    int pos_id = proj1 >= 0 ? ve.v1 : ve.v2;
                    edge_positive_pt_ids_[e_count] = pos_id;
                    // ptid_set.insert(pos_id);

                    face_center[0] += (vp_list[3 * pos_id]     + inter_x); 
                    face_center[1] += (vp_list[3 * pos_id + 1] + inter_y);
                    face_center[2] += (vp_list[3 * pos_id + 2] + inter_z);

                    trunc_facet_center_[0] += inter_x;
                    trunc_facet_center_[1] += inter_y;
                    trunc_facet_center_[2] += inter_z;

                    face_edge_inter_count ++;
                    // printf(" face_edge_inter_count -------- : %d \n", face_edge_inter_count);
                    intersect_count++;
                    p_count += 2;
                } 
                else if(proj1 >= 0 && proj2 >= 0)
                {
                    // printf("----------- edge count %d proj 1 %f, proj 2 %f \n", e_count, proj1, proj2);
                    edge_sign_vals_[e_count] = IN_UNION;
                    face_center[0] += (vp_list[3 * ve.v1]     + vp_list[3 * ve.v2]); 
                    face_center[1] += (vp_list[3 * ve.v1 + 1] + vp_list[3 * ve.v2 + 1]);
                    face_center[2] += (vp_list[3 * ve.v1 + 2] + vp_list[3 * ve.v2 + 2]);
                    p_count += 2;
                    positive_count ++;
                } else {
                    edge_sign_vals_[e_count] = OUT_UNION;
                }
            } else {
                edge_sign_vals_[e_count] = IN_VALID;
                return 0;
            }
            e_count ++;
        }
        if(p_count> 0)
        {
            face_centers_[3*i]     = face_center[0] / double(p_count);
            face_centers_[3*i + 1] = face_center[1] / double(p_count);
            face_centers_[3*i + 2] = face_center[2] / double(p_count);
        }
        // printf(" --------------- intersect_count : %d \n", intersect_count);
        if(face_edge_inter_count > 0)
        {
            face_sign_vals_[i] = INTERSECT;
            cell_center[0] += face_centers_[3*i];
            cell_center[1] += face_centers_[3*i + 1];
            cell_center[2] += face_centers_[3*i + 2];
            valid_count ++;
        } else if(positive_count > 0) {
            face_sign_vals_[i] = IN_UNION;
            cell_center[0] += face_centers_[3*i];
            cell_center[1] += face_centers_[3*i + 1];
            cell_center[2] += face_centers_[3*i + 2];
            valid_count ++;
        } else {
            face_sign_vals_[i] = OUT_UNION;
        }
        // face_sign_vals_[i] = intersect_count == 2 ? INTERSECT : (positive_count > 0 ? IN_UNION : OUT_UNION);
    }
    if(intersect_count > 0)
    {
        trunc_facet_center_[0] /= (double(intersect_count) );
        trunc_facet_center_[1] /= (double(intersect_count) );
        trunc_facet_center_[2] /= (double(intersect_count) );
        cell_center[0] +=  trunc_facet_center_[0];
        cell_center[1] +=  trunc_facet_center_[1];
        cell_center[2] +=  trunc_facet_center_[2];

        valid_count ++;
    }
    
    
    if(valid_count == 0) return 0;
    cell_center[0] /= (double(valid_count) );
    cell_center[1] /= (double(valid_count) );
    cell_center[2] /= (double(valid_count) );

    e_count = 0;
    double cell_volume = 0; 
    int new_tet_count = 0; 
    for(size_t i = 0; i < f_num; ++i)
    {
        size_t f_id = v_cell[i + 1];
        const auto& facet = vf_list[f_id];
        size_t e_num = facet.elist[0];
        if(face_sign_vals_[i] == OUT_UNION || face_sign_vals_[i] == IN_VALID)  
        {
            e_count += e_num;
            continue;
        }
        
        for(size_t j = 0; j < e_num; ++j)
        {
            // const auto& cur_edge = ve_list[facet.elist[j + 1]];

            size_t e_id = facet.elist[1 + j];
            const auto& ve = ve_list[e_id];

            if(edge_sign_vals_[e_count] == IN_UNION)
            {
                double volume = CalTetrahedronVolumeDet(&vp_list[3 * ve.v1], 
                                                        &vp_list[3 * ve.v2], 
                                                        &face_centers_[3 * i],
                                                        &cell_center[0]);
                cell_volume += volume;
                new_tet_count++;
            } else if(edge_sign_vals_[e_count] == INTERSECT)
            {
                int positivie_id = edge_positive_pt_ids_[e_count]; 
                int intersect_id = edge_intersect_pt_ids_[e_count];
                double volume = CalTetrahedronVolumeDet(&vp_list[3 * positivie_id], 
                                                        &intersect_pts_[3 * intersect_id],
                                                        &face_centers_[3 * i],
                                                        &cell_center[0]);
                cell_volume += volume;
                new_tet_count++;
            }
            e_count++;
        }
        if(face_sign_vals_[i] == INTERSECT)
        {
            int intersect_id1 = face_insersect_pt_ids_[2 * i];
            int intersect_id2 = face_insersect_pt_ids_[2 * i + 1];
            double volume = CalTetrahedronVolumeDet(&intersect_pts_[3 * intersect_id1], 
                                                    &intersect_pts_[3 * intersect_id2], 
                                                    &face_centers_[3 * i],
                                                    &cell_center[0]);
            cell_volume += volume;

            volume = CalTetrahedronVolumeDet(&intersect_pts_[3 * intersect_id1], 
                                                &intersect_pts_[3 * intersect_id2], 
                                                &trunc_facet_center_[0],
                                                &cell_center[0]);
            cell_volume += volume;
            new_tet_count += 2;
        }
    }
    // printf("------------ new_tet_count : %d \n", new_tet_count);
if(0)
{
    std::ofstream out_file;
    std::string out_path = "../out/centers.obj";
    out_file.open(out_path);

    for(size_t i = 0; i < f_num; ++i)
    {
        if(face_sign_vals_[i] != OUT_UNION) 
        {
            out_file << "v " << face_centers_[3*i] << " " << face_centers_[3*i + 1] 
                        << " " << face_centers_[3*i + 2] << " 0 0 0" << std::endl;
        }
    }
    out_file << "v " << trunc_facet_center_[0] << " " << trunc_facet_center_[1] 
                        << " " << trunc_facet_center_[2] << " 0.2 0.1 0.4" << std::endl;
    out_file << "v " << cell_center[0] << " " << cell_center[1] 
                        << " " << cell_center[2] << " 1 0 0" << std::endl;

    for(size_t i = 0; i < intersect_count; ++i)
    {
        out_file << "v " << intersect_pts_[3*i] << " " << intersect_pts_[3*i + 1] 
                        << " " << intersect_pts_[3*i + 2] << " 0 0 1" << std::endl;
    }
    out_file.close();

    // out_path = "../out/poly_pts.obj";
    // out_file.open(out_path);
    // for(auto p_id : ptid_set)
    // {
    //     out_file << "v " << vp_list[3*p_id] << " " 
    //                      << vp_list[3*p_id + 1] << " " 
    //                      << vp_list[3*p_id + 2] << std::endl;
    // }
    // out_file.close();

    out_path = "../out/polyhedron.obj";
    out_file.open(out_path);

    auto out_e_count = 0;
    for(size_t i = 0; i < f_num; ++i)
    {
        size_t f_id = v_cell[i + 1];
        auto facet = vf_list[f_id];
        size_t e_num = facet.elist[0];
        for(size_t j = 0; j < e_num; ++j)
        {
            int e_id = facet.elist[j + 1];
            auto ve = ve_list[e_id];
            out_file << "v " << vp_list[3*ve.v1] << " " 
                         << vp_list[3*ve.v1 + 1] << " " 
                         << vp_list[3*ve.v1 + 2] << std::endl;

            out_file << "v " << vp_list[3*ve.v2] << " " 
                         << vp_list[3*ve.v2 + 1] << " " 
                         << vp_list[3*ve.v2 + 2] << std::endl;
            out_file << "l " << std::to_string(1 + out_e_count) << " " << std::to_string(2 + out_e_count) << std::endl;            
            out_e_count += 2;

        }
    }
    out_file.close();
}


    return cell_volume;
}



double VoronoiGen::CalTruncatedCellVolumePassOMP(tetgenmesh::point in_pt, tetgenmesh::point nei_pt, int thread_id)
{
    double plane_mid_x = (in_pt[0] + nei_pt[0])/2.0;
    double plane_mid_y = (in_pt[1] + nei_pt[1])/2.0;
    double plane_mid_z = (in_pt[2] + nei_pt[2])/2.0;
    double p_nx  = (in_pt[0] - nei_pt[0]);
    double p_ny  = (in_pt[1] - nei_pt[1]);
    double p_nz  = (in_pt[2] - nei_pt[2]);
    double pn_len = sqrt(p_nx * p_nx + p_ny * p_ny + p_nz * p_nz);
    if(pn_len > 0) 
    {
        p_nx /= pn_len;
        p_ny /= pn_len;
        p_nz /= pn_len;
    }
    // VoroPlane v_plane(plane_mid_x, plane_mid_y, plane_mid_z, p_nx, p_ny, p_nz);
    // SavePlane(v_plane);
    // if(point_id_map_.find(nei_pt) == point_id_map_.end())
    // return 0;
    const auto vc_id   = point_id_map_[nei_pt];
    const auto v_cell  = voronoi_data_.vcelllist[vc_id];
    const auto vf_list = voronoi_data_.vfacetlist;
    const auto ve_list = voronoi_data_.vedgelist;
    const auto vp_list = voronoi_data_.vpointlist;
    // printf("pt size 0001 : %ld \n", points_.size());
    if(v_cell == NULL) return 0;
    const int f_num = v_cell[0];

    int cur_opm_id = thread_id;
    double* intersect_pts = intersect_pts_omp_[cur_opm_id];
    double* trunc_center = truncated_centers_omp_[cur_opm_id];
    double** trunc_tets = truncated_tets_omp_[cur_opm_id];
    int* face_tet_count = face_tet_count_omp_[cur_opm_id];

    int intersect_total_count = 0;
    int tet_total_count = 0;
    for(size_t i = 0; i < f_num; ++i)
    {
        size_t f_id = v_cell[i + 1];
        const auto& facet = vf_list[f_id];
        size_t e_num = facet.elist[0];
        int tet_count = 0;
        int intersect_count = 0; 
        for(size_t j = 0; j < e_num; ++j)
        {
            // printf("edge id : %ld \n", j);
            size_t e_id = facet.elist[1 + j];
            const auto& ve = ve_list[e_id];
            if(ve.v1 != -1 && ve.v2 != -1)
            {
                double dx = vp_list[3 * ve.v1]     - plane_mid_x;
                double dy = vp_list[3 * ve.v1 + 1] - plane_mid_y;
                double dz = vp_list[3 * ve.v1 + 2] - plane_mid_z;
                double proj1 = dx * p_nx + dy * p_ny + dz * p_nz;

                dx = vp_list[3 * ve.v2]     - plane_mid_x;
                dy = vp_list[3 * ve.v2 + 1] - plane_mid_y;
                dz = vp_list[3 * ve.v2 + 2] - plane_mid_z;
                double proj2 = dx * p_nx + dy * p_ny + dz * p_nz;
                if(proj1 * proj2 < 0 )
                {
                    double r1 = abs(proj1) /(abs(proj1) + abs(proj2));
                    double r2 = abs(proj2) /(abs(proj1) + abs(proj2));
                    double inter_x = r2 * vp_list[3 * ve.v1]     + r1 * vp_list[3 * ve.v2];
                    double inter_y = r2 * vp_list[3 * ve.v1 + 1] + r1 * vp_list[3 * ve.v2 + 1];
                    double inter_z = r2 * vp_list[3 * ve.v1 + 2] + r1 * vp_list[3 * ve.v2 + 2];
                    intersect_pts[3 * (intersect_total_count + intersect_count)]     = inter_x;
                    intersect_pts[3 * (intersect_total_count + intersect_count) + 1] = inter_y;
                    intersect_pts[3 * (intersect_total_count + intersect_count) + 2] = inter_z;
                    int posi_id = proj1 > 0 ? ve.v1 : ve.v2;
                    trunc_tets[(tet_total_count + tet_count)*2] = &vp_list[3 * posi_id];
                    trunc_tets[(tet_count + tet_total_count)*2 + 1] = &intersect_pts[3 * (intersect_total_count + intersect_count)];

                    intersect_count++;
                    tet_count ++;

                } 
                else if(proj1 >= 0 && proj2 >= 0)
                {
                    trunc_tets[(tet_total_count + tet_count)*2] = &vp_list[3 * ve.v2];
                    trunc_tets[(tet_total_count + tet_count)*2 + 1] = &vp_list[3 * ve.v1];
                    tet_count ++;
                } 
            } 
        }
        // printf("intersect_count num : %d \n", intersect_count);
        if(intersect_count > 0)
        {
            trunc_tets[(tet_total_count + tet_count)*2] = &intersect_pts[3 * (intersect_total_count)];
            trunc_tets[(tet_total_count + tet_count)*2 + 1] = &intersect_pts[3 * intersect_total_count + 3];
            
            intersect_total_count += 2;
            tet_count ++;
        }
        // printf("tet_count num : %d \n", tet_count);
        face_tet_count[i] = tet_count;
        tet_total_count += tet_count;
    }
    double cell_volume = 0;

    if(intersect_total_count > 0)
    {
        const double* trunc_pt = intersect_pts; 
        int cur_tet_id = 0;
        for(size_t i = 0; i < f_num; ++i)
        {
            const int f_tet_num = face_tet_count[i];
            const auto f_pt = trunc_tets[cur_tet_id*2];
            // printf("cur face f_tet_num : %d \n", f_tet_num);
            for(int ti = 0; ti < f_tet_num; ++ti )
            {
                double volume = CalTetrahedronVolumeDet(trunc_tets[cur_tet_id*2], 
                                                        trunc_tets[cur_tet_id*2 + 1], 
                                                        f_pt,
                                                        trunc_pt);
                cell_volume += volume;
                cur_tet_id ++;
            }
        }
        return cell_volume;
    }
    return 0;
}


double VoronoiGen::CalTruncatedCellVolumePass2(tetgenmesh::point in_pt, tetgenmesh::point nei_pt)
{
    double plane_mid_x = (in_pt[0] + nei_pt[0])/2.0;
    double plane_mid_y = (in_pt[1] + nei_pt[1])/2.0;
    double plane_mid_z = (in_pt[2] + nei_pt[2])/2.0;
    double p_nx  = (in_pt[0] - nei_pt[0]);
    double p_ny  = (in_pt[1] - nei_pt[1]);
    double p_nz  = (in_pt[2] - nei_pt[2]);
    double pn_len = sqrt(p_nx * p_nx + p_ny * p_ny + p_nz * p_nz);
    if(pn_len > 0) 
    {
        p_nx /= pn_len;
        p_ny /= pn_len;
        p_nz /= pn_len;
    }
    // VoroPlane v_plane(plane_mid_x, plane_mid_y, plane_mid_z, p_nx, p_ny, p_nz);
    // SavePlane(v_plane);
    // if(point_id_map_.find(nei_pt) == point_id_map_.end())
    // return 0;
    auto vc_id = point_id_map_[nei_pt];
    auto v_cell  = voronoi_data_.vcelllist[vc_id];
    auto vf_list = voronoi_data_.vfacetlist;
    auto ve_list = voronoi_data_.vedgelist;
    auto vp_list = voronoi_data_.vpointlist;
    // printf("pt size 0001 : %ld \n", points_.size());
    if(v_cell == NULL) return 0;
    int f_num = v_cell[0];
    int e_count = 0;
    // double cell_center[3] = {0, 0, 0};
    int valid_count = 0;
    int intersect_count = 0; 
    double trunc_facet_center_[3] = {0, 0, 0};
    const auto& vc_pids = vorocell_pids_[vc_id];
    const auto& eids = vorocell_eids_[vc_id];

if(1)
{
    size_t t_i;
    // #pragma omp parallel for private(t_i) shared(vc_pids, v_sign_vals_, vp_list)
    for( t_i = 0; t_i < vc_pids.size(); ++t_i)
    {
        int cur_pid = vc_pids[t_i];
        double dx = vp_list[3 * cur_pid]     - plane_mid_x;
        double dy = vp_list[3 * cur_pid + 1] - plane_mid_y;
        double dz = vp_list[3 * cur_pid + 2] - plane_mid_z;
        v_sign_vals_[cur_pid] = dx * p_nx + dy * p_ny + dz * p_nz;
    }
}

    
    // std::cout <<std::endl;
    // return 0;
    // printf("v_sign_vals_ size  %ld \n", v_sign_vals_.size());
    // std::cout << " eids : " ;

    // std::string intersect_pt_path = out_dir_ + std::to_string(vc_id) + "intersecet_pts.obj";
    // std::ofstream inter_pt_file;
    // inter_pt_file.open(intersect_pt_path);

    for(size_t i = 0; i < eids.size(); ++i)
    {
        int cur_eid = eids[i];
        const auto& ve = ve_list[cur_eid];
        if(ve.v1 != -1 && ve.v2 != -1)
        {
            double proj1 = v_sign_vals_[ve.v1];
            double proj2 = v_sign_vals_[ve.v2];
            if(proj1 * proj2 <= 0)
            {
                double r1 = abs(proj1) /(abs(proj1) + abs(proj2));
                double r2 = abs(proj2) /(abs(proj1) + abs(proj2));
                double inter_x = r2 * vp_list[3 * ve.v1]     + r1 * vp_list[3 * ve.v2];
                double inter_y = r2 * vp_list[3 * ve.v1 + 1] + r1 * vp_list[3 * ve.v2 + 1];
                double inter_z = r2 * vp_list[3 * ve.v1 + 2] + r1 * vp_list[3 * ve.v2 + 2];
                ve_intersect_pts_[3 * cur_eid]     = inter_x;
                ve_intersect_pts_[3 * cur_eid + 1] = inter_y;
                ve_intersect_pts_[3 * cur_eid + 2] = inter_z;

                trunc_facet_center_[0] +=  inter_x;
                trunc_facet_center_[1] +=  inter_y;
                trunc_facet_center_[2] +=  inter_z;
                intersect_count ++;
                e_sign_vals_[cur_eid] = INTERSECT;
            } 
            else if(proj1 >0 && proj2 > 0)
            {
                e_sign_vals_[cur_eid] = IN_UNION;
            } else {
                e_sign_vals_[cur_eid] = OUT_UNION;
            }

            // std::cout << "e_sign_vals_[cur_eid] " << e_sign_vals_[cur_eid] << std::endl;
        }    
    }
    // return 0;
    // inter_pt_file.close();
    // std::cout << std::endl;
    // printf("e_sign_vals_ size  %ld \n", e_sign_vals_.size());
    // printf("intersect_count   %d \n", intersect_count)

    if(intersect_count > 0)
    {
        trunc_facet_center_[0] /= (double(intersect_count) );
        trunc_facet_center_[1] /= (double(intersect_count) );
        trunc_facet_center_[2] /= (double(intersect_count) );
        valid_count ++;
    }

    // printf("valid_count   %d \n", valid_count);
    
    // int max_triangles = 30000;
    // triagnles.resize(max_triangles);
    int tri_count = 0;
    double volume_sum = 0;
    
    for(size_t i = 0; i < f_num; ++i)
    {
        size_t f_id = v_cell[i + 1];
        const auto& facet = vf_list[f_id];
        size_t e_num = facet.elist[0];
        bool find_first_vid= false;
        int start_id = 0;
        std::vector<double*> intersect_triagnles;
        // int intersect_cout = 0;
        for(size_t j = 0; j < e_num; ++j)
        {
            size_t e_id = facet.elist[1 + j];
            const auto& ve = ve_list[e_id];
            if(e_sign_vals_[e_id] == IN_UNION)
            {
                if(find_first_vid)
                {
                    // tet_pts_[tri_count * 4*3 ] = vp_list[3 * start_id];
                    // tet_pts_[tri_count * 4*3 + 1] = vp_list[3 * start_id + 1];
                    // tet_pts_[tri_count * 4*3 + 2] = vp_list[3 * start_id + 2];

                    // tet_pts_[tri_count * 4*3 + 3] = vp_list[3 * ve.v1];
                    // tet_pts_[tri_count * 4*3 + 4] = vp_list[3 * ve.v1 + 1];
                    // tet_pts_[tri_count * 4*3 + 5] = vp_list[3 * ve.v1 + 2];

                    // tet_pts_[tri_count * 4*3 + 6] = vp_list[3 * ve.v2];
                    // tet_pts_[tri_count * 4*3 + 7] = vp_list[3 * ve.v2 + 1];
                    // tet_pts_[tri_count * 4*3 + 8] = vp_list[3 * ve.v2 + 2];

                    // tet_pts_[tri_count * 4*3 + 9] = trunc_facet_center_[0];
                    // tet_pts_[tri_count * 4*3 + 10] = trunc_facet_center_[1];
                    // tet_pts_[tri_count * 4*3 + 11] = trunc_facet_center_[2];

                    // triagnles_[tri_count *3]     = &vp_list[3 * start_id];
                    // triagnles_[tri_count *3 + 1] = &vp_list[3 * ve.v1];
                    // triagnles_[tri_count *3 + 2] = &vp_list[3 * ve.v2];
                    volume_sum += CalTetrahedronVolumeDet(&vp_list[3 * start_id],
                                                &vp_list[3 * ve.v1],
                                                &vp_list[3 * ve.v2],
                                                trunc_facet_center_);
                    tri_count ++;
                } else {
                    find_first_vid = true;
                    start_id = ve.v2;
                    
                }
            } else if(e_sign_vals_[e_id] == INTERSECT)
            {
                if(find_first_vid)
                {
                //    triagnles_[tri_count *3] = &vp_list[3 * start_id];
                    // tet_pts_[tri_count * 4*3 ] = vp_list[3 * start_id];
                    // tet_pts_[tri_count * 4*3 + 1] = vp_list[3 * start_id + 1];
                    // tet_pts_[tri_count * 4*3 + 2] = vp_list[3 * start_id + 2];

                    if(v_sign_vals_[ve.v1] > 0)
                    {
                        // triagnles_[tri_count *3 + 1] = &vp_list[3 * ve.v1];
                        // tet_pts_[tri_count * 4*3 + 3] = vp_list[3 * ve.v1];
                        // tet_pts_[tri_count * 4*3 + 4] = vp_list[3 * ve.v1 + 1];
                        // tet_pts_[tri_count * 4*3 + 5] = vp_list[3 * ve.v1 + 2];
                        volume_sum += CalTetrahedronVolumeDet(&vp_list[3 * start_id],
                                                &vp_list[3 * ve.v1],
                                                &ve_intersect_pts_[3 * e_id],
                                                trunc_facet_center_);
                    } else {
                        // tet_pts_[tri_count * 4*3 + 3] = vp_list[3 * ve.v2];
                        // tet_pts_[tri_count * 4*3 + 4] = vp_list[3 * ve.v2 + 1];
                        // tet_pts_[tri_count * 4*3 + 5] = vp_list[3 * ve.v2 + 2];
                        // triagnles_[tri_count *3 + 1] = &vp_list[3 * ve.v2];
                        volume_sum += CalTetrahedronVolumeDet(&vp_list[3 * start_id],
                                                &vp_list[3 * ve.v2],
                                                &ve_intersect_pts_[3 * e_id],
                                                trunc_facet_center_);
                    }
                    // tet_pts_[tri_count * 4*3 + 6] = ve_intersect_pts_[3 * e_id];
                    // tet_pts_[tri_count * 4*3 + 7] = ve_intersect_pts_[3 * e_id + 1];
                    // tet_pts_[tri_count * 4*3 + 8] = ve_intersect_pts_[3 * e_id + 2];

                    // triagnles_[tri_count *3 + 2] = &ve_intersect_pts_[3 * e_id];
                    tri_count ++;
                } else {
                    if(v_sign_vals_[ve.v1] > 0) 
                    {
                        start_id = ve.v1;
                    } else {
                        start_id = ve.v2;
                    }
                    find_first_vid = true;
                }
                intersect_triagnles.push_back(&ve_intersect_pts_[3 * e_id]);
            } 
        }

        
        if(find_first_vid)
        {
            // std::cout << " intersect_triagnles.size()  "<< intersect_triagnles.size()  << std::endl;
            if(intersect_triagnles.size() == 2)
            {
                // triagnles_[tri_count *3] = intersect_triagnles[0];
                // triagnles_[tri_count *3 + 1] = intersect_triagnles[1];
                // triagnles_[tri_count *3 + 2] = &vp_list[3 * start_id];

                double volume = CalTetrahedronVolumeDet(intersect_triagnles[0],
                                                intersect_triagnles[1],
                                                &vp_list[3 * start_id],
                                                trunc_facet_center_);
                volume_sum += volume;
                // tri_count ++;
            }
        } 
    }
    
    

    // printf("triangle size : %ld \n", triagnles.size());
    // printf("tri_count size : %d \n", tri_count);
// #pragma omp parallel for
//     for(size_t i = 0; i < tri_count; ++i)
//     {
//         // double volume = CalTetrahedronVolumeDet(triagnles_[3*i],
//         //                                         triagnles_[3*i + 1],
//         //                                         triagnles_[3*i + 2],
//         //                                         trunc_facet_center_);

//         double volume = CalTetrahedronVolumeDet(&tet_pts_[12*i], &tet_pts_[12*i+ 3],  &tet_pts_[12*i+ 6],  &tet_pts_[12*i+ 9]);                                        

        
//         volume_vec_[i] = volume;
//     }

// for(size_t i = 0; i < tri_count; ++i)
// {
//     volume_sum += volume_vec_[i];
// }

    return volume_sum;
}

void VoronoiGen::OutputTetMesh(const std::string& tetmesh_path)
{
    std::ofstream mesh_file;
    mesh_file.open(tetmesh_path);
    // for(auto pt : points_)
    // {
    //     mesh_file << "v " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl; 
    // }

    tetMesh_.points->traversalinit();
    tetgenmesh::point ptloop;
    ptloop = tetMesh_.pointtraverse();
    while (ptloop != (tetgenmesh::point) NULL) {
        mesh_file << "v " << ptloop[0] << " " << ptloop[1] << " " << ptloop[2] << std::endl; 
        ptloop = tetMesh_.pointtraverse();
    }

    tetgenmesh::triface tetloop;
    tetMesh_.tetrahedrons->traversalinit();
    tetloop.tet = tetMesh_.alltetrahedrontraverse();
    size_t t_count = 0;
    while(tetloop.tet != (tetgenmesh::tetrahedron*) NULL)
    {
        auto pts = (tetgenmesh::point *) tetloop.tet;
        if(pts[4] == (tetgenmesh::point) NULL || pts[5] == (tetgenmesh::point) NULL ||
          pts[6] == (tetgenmesh::point) NULL || pts[7] == (tetgenmesh::point) NULL) 
        {
            tetloop.tet = tetMesh_.alltetrahedrontraverse();
            continue;
        }
      
        int idx0 = tetMesh_.pointmark(pts[4]); 
        int idx1 = tetMesh_.pointmark(pts[5]); 
        int idx2 = tetMesh_.pointmark(pts[6]); 
        int idx3 = tetMesh_.pointmark(pts[7]); 

        if(idx0 == -1 || idx1 == -1 || idx2 == -1 || idx3 == -1)
        {
            tetloop.tet = tetMesh_.alltetrahedrontraverse();
            continue;
        }

        mesh_file << "l " << idx0 + 1 << " " << idx1 + 1 << std::endl;
        mesh_file << "l " << idx0 + 1 << " " << idx2 + 1 << std::endl;
        mesh_file << "l " << idx0 + 1 << " " << idx3 + 1 << std::endl;
        mesh_file << "l " << idx1 + 1 << " " << idx2 + 1 << std::endl;
        mesh_file << "l " << idx1 + 1 << " " << idx3 + 1 << std::endl;
        mesh_file << "l " << idx2 + 1 << " " << idx3 + 1 << std::endl;
        tetloop.tet = tetMesh_.alltetrahedrontraverse();
    }

    // tetMesh_.subsegs->traversalinit();
    // tetgenmesh::face segloop;
    // segloop.sh = tetMesh_.shellfacetraverse(tetMesh_.subsegs);
    // while (segloop.sh != (tetgenmesh::shellface *) NULL) {
    //     auto torg = tetMesh_.sorg(segloop);
    //     auto tdest = tetMesh_.sdest(segloop);
    //     int idx0 = tetMesh_.pointmark(torg); 
    //     int idx1 = tetMesh_.pointmark(tdest);
    //     mesh_file << "l " << idx0 + 1 << " " << idx1 + 1 << std::endl; 
    //     segloop.sh = tetMesh_.shellfacetraverse(tetMesh_.subsegs);
    // }
    mesh_file.close();
}


void VoronoiGen::OutputVoronisMesh()
{
    std::string out_path = out_dir_ + "voronoi_cell.obj";
    std::string out_path2 = out_dir_ + "voronoi_split_pts.obj";
    std::cout << " voronoi out put size : " << out_path << std::endl;
    // std::cout << " voronoi tet size : " << tetMesh_.tetrahedrons->items << std::endl;
    ofstream myfile;
    myfile.open(out_path);
    // ofstream myfile_split;
    // myfile_split.open(out_path2);
    size_t vp_num = voronoi_data_.numberofvpoints;
    // std::cout << " voronoi point size : " << vp_num << std::endl;
    auto& vplist = voronoi_data_.vpointlist;
    for(size_t i = 0; i < vp_num; ++i)
    {
        myfile << "v " << vplist[3*i] << " " << vplist[3*i + 1] << " " << vplist[3*i + 2] << std::endl;
    }
    size_t e_num = voronoi_data_.numberofvedges;
    auto& v_edges = voronoi_data_.vedgelist;
    // int incre_pt_count = 0;
    // for(size_t i =0; i < e_num; ++i)
    // {
    //     if( v_edges[i].v2 == -1)
    //     {
    //         double nx =  vplist[3*v_edges[i].v1]     + v_edges[i].vnormal[0] * 10;
    //         double ny =  vplist[3*v_edges[i].v1 + 1] + v_edges[i].vnormal[1] * 10;
    //         double nz =  vplist[3*v_edges[i].v1 + 2] + v_edges[i].vnormal[2] * 10;
    //         myfile << "v " << nx << " " << ny << " " << nz << std::endl;
    //         v_edges[i].v2 = vp_num + incre_pt_count;
    //         incre_pt_count ++;
    //     }
    // }

    for(size_t i =0; i < e_num; ++i)
    {
        if(v_edges[i].v1 == -1 || v_edges[i].v2 == -1) continue;
        myfile << "l " << v_edges[i].v1 + 1 << " " << v_edges[i].v2 + 1 << std::endl;
    }
    myfile.close();

    // size_t cell_num = voronoi_data_.numberofvcells;
    // size_t face_num = voronoi_data_.numberofvfacets;
    // auto v_faces = voronoi_data_.vfacetlist;
    // auto v_cells = voronoi_data_.vcelllist;
    // printf("pt size : %ld \n", points_.size());
    // size_t f_vcount = 0;
    // for(size_t i = 0; i < face_num; ++i)
    // {
    //     f_vcount += v_faces[i].elist[0];
    // }
    // f_vcount += face_num * 3;
    // printf("voro face memory size : %ld \n", f_vcount * 4);

    // size_t c_fcount = 0;
    // for(size_t i = 0; i < cell_num; ++i)
    // {
    //     if(v_cells[i])
    //     {
    //         c_fcount += v_cells[i][0];
    //     }
    // }
    // c_fcount += cell_num ;
    // printf("voro cell memory size : %ld \n", c_fcount * 4);

    // tetMesh_.points->traversalinit();
    // tetgenmesh::point ploop = tetMesh_.pointtraverse();
    // while(ploop != (tetgenmesh::point)NULL)
    // {
    //     // printf("pt size 00 : %ld \n", points_.size());
    //     const auto p_type = tetMesh_.pointtype(ploop);
    //     if ((p_type == tetgenmesh::UNUSEDVERTEX)) 
    //     {
    //         ploop = tetMesh_.pointtraverse();
    //         continue;
    //     } 
    //     std::vector<tetgenio::voroedge*> edge_list;
    //     GetVoroCellEdgeList(ploop, edge_list);
    //     std::vector<double> pt_list;
    //     GetVoroCellPtList(edge_list, pt_list);

    //     VoroPlane v_plane;
    //     v_plane.nx = 0.5;
    //     v_plane.ny = 0.5;
    //     v_plane.nz = 0.5;
    //     v_plane.px = ploop[0] + 0.05;
    //     v_plane.py = ploop[1] + 0.05;
    //     v_plane.pz = ploop[2] + 0.05;
    //     double cur_pt[3];

    //     cur_pt[0] = v_plane.px;
    //     cur_pt[1] = v_plane.py;
    //     cur_pt[2] = v_plane.pz;
    //     auto split_pts = GetPolygonSplitPts(cur_pt, ploop);

    //     // CalSplitPolyhedronVolume(split_pts);

    //     // int v_id = split_pts.size()/3;
    //     // for(auto edge : edge_list)
    //     // {
    //     //     if(CalPlaneEdgeIntersection(edge, v_plane))
    //     //     {
    //     //         split_pts.push_back(v_plane.intPx);
    //     //         split_pts.push_back(v_plane.intPy);
    //     //         split_pts.push_back(v_plane.intPz);
    //     //         // myfile_split << "v " << v_plane.intPx << " " << v_plane.intPy << " " << v_plane.intPz << std::endl;
    //     //         // v_id++;
    //     //     }
    //     // }

    //     // for(size_t p_id = 0; p_id < split_pts.size()/3; ++p_id)
    //     // {
    //     //     myfile_split << "v " << split_pts[3*p_id]<< " " 
    //     //     << split_pts[3*p_id + 1] << " " <<  split_pts[3*p_id + 2] << std::endl;
    //     // }

    //     // myfile_split << "f ";
    //     // for(size_t v_i = 0; v_i < v_id; ++v_i)
    //     // {
    //     //     myfile_split << v_i + 1 << " ";
    //     // }
    //     // myfile_split <<std::endl;
    //     // myfile_split.close();
      
    //     // printf("pt size 02 : %ld \n", points_.size());
    //     // for(auto edge : edge_list)
    //     // {
    //     //     myfile << "l " << edge->v1 + 1 << " " << edge->v2 + 1 << std::endl;
    //     // }
    //     // printf("pt size : %ld \n", points_.size());
    //     ploop = tetMesh_.pointtraverse();
    //     // break;
    // }

    // for(auto& n_pt : points_)
    // {
    //     std::vector<tetgenio::voroedge*> edge_list;
    //     GetVoroCellEdgeList(n_pt, edge_list);
    //     for(auto edge : edge_list)
    //     {
    //         myfile << "l " << edge->v1 + 1 << " " << edge->v2 + 1 << std::endl;
    //     }
    //     // break;
    // }
    
    // for(voronoi_data_.vpointlist)

}



void VoronoiGen::PrecomputeVoroData()
{
    auto vf_list = voronoi_data_.vfacetlist;
    auto ve_list = voronoi_data_.vedgelist;
    auto vp_list = voronoi_data_.vpointlist;
    auto vc_list = voronoi_data_.vcelllist;
    // vids_status_.resize(voronoi_data_.numberofvpoints, false);
    v_sign_vals_.resize(voronoi_data_.numberofvpoints, 0);

    // eids_status_.resize(voronoi_data_.numberofvedges, false);
    e_sign_vals_.resize(voronoi_data_.numberofvedges);
    ve_intersect_pts_.resize(voronoi_data_.numberofvedges * 3);
    int vc_num = voronoi_data_.numberofvcells;

    vorocell_pids_.resize(vc_num);
    vorocell_eids_.resize(vc_num);

    for(int vci =0; vci < vc_num; ++vci)
    {
        auto v_cell  = voronoi_data_.vcelllist[vci];
        // printf("pt size 0001 : %ld \n", points_.size());
        if(v_cell == NULL) continue;
        int f_num = v_cell[0];
        std::set<int> vcell_pt_list;
        std::set<int> vcell_edge_list ;
        for(size_t i = 0; i < f_num; ++i)
        {
            size_t f_id = v_cell[i + 1];
            auto facet = vf_list[f_id];
            size_t e_num = facet.elist[0];
            arma::vec3 f_center;
            for(size_t j = 0; j < e_num; ++j)
            {
                int e_id = facet.elist[1 + j];
                vcell_edge_list.insert(e_id);
                const auto& ve = ve_list[e_id];

                
                vcell_pt_list.insert(ve.v1);
                vcell_pt_list.insert(ve.v2);
            }
        }
        vorocell_pids_[vci] = std::vector<int>(vcell_pt_list.begin(), vcell_pt_list.end());
        vorocell_eids_[vci] = std::vector<int>(vcell_edge_list.begin(), vcell_edge_list.end());;
        
    }   
    // for(auto vc_ids : vorocell_pids_)
    // {
    //     for(auto id : vc_ids)
    //     {
    //         if(id > voronoi_data_.numberofvpoints)
    //         {
    //             std::cout << " too large ids : " << id << std::endl;
    //         }
    //     }
    // }
    

    // for(int i =0; i < vc_num; ++i)
    // {
    //     eids_status_.resize(voronoi_data_.numberofvedges, false);
    //     vids_status_.resize(voronoi_data_.numberofvpoints, false);
    //     // std::vector<int> vids;
    //     // std::vector<int> eids;
    //     auto vcell = vc_list[i];
    //     if(vcell != NULL)
    //     {
    //         int f_num = vcell[0];
    //         for(size_t fi = 0; fi < f_num; ++fi)
    //         {
    //             int f_id = vcell[fi + 1];
    //             const auto& facet = vf_list[f_id];
    //             size_t e_num = facet.elist[0];
    //             for(size_t ei = 0; ei < e_num; ++ei)
    //             {
    //                 size_t e_id = facet.elist[1 + ei];
    //                 if(eids_status_[e_id]) continue;
    //                 eids_status_[e_id] = true;
    //                 vorocell_eids_[i].push_back(e_id);
    //                 const auto& ve = ve_list[e_id];
    //                 if(ve.v1 == -1 || ve.v2 == -1) continue;
    //                 if(!vids_status_[ve.v1])
    //                 {
    //                     vorocell_pids_[i].push_back(ve.v1);
    //                     vids_status_[ve.v1] = true;
    //                 }
    //                 if(!vids_status_[ve.v2])
    //                 {
    //                     vorocell_pids_[i].push_back(ve.v2);
    //                     vids_status_[ve.v2] = true;
    //                 }
    //             } 
    //         }
    //     }
        // for(auto vid : vorocell_pids_[i])
        // {
        //     vids_status_[vid] = false;
        // }
        // for(auto eid : vorocell_eids_[i])
        // {
        //     eids_status_[eid] = false;
        // }
        // vorocell_pids_[i] = vids;
        // vorocell_eids_[i] = eids;
    // }
}

void VoronoiGen::GetVertexStar(tetgenmesh::point &p_st, 
            std::set<tetgenmesh::point>& candid_pts, int level = 1)
{
    tetlist_->restart();
    ptlist_->restart();
    tetMesh_.getvertexstar(1, p_st, tetlist_,  ptlist_, NULL);
    candid_pts.insert(p_st);
    for (int i = 0; i < ptlist_->objects; i++) {
        auto pt = * (tetgenmesh::point *) ptlist_-> lookup(i);
        if ((tetMesh_.pointtype(pt) == tetgenmesh::verttype::UNUSEDVERTEX) ||
        (tetMesh_.pointtype(pt) == tetgenmesh::verttype::DUPLICATEDVERTEX) ||
        (tetMesh_.pointtype(pt) == tetgenmesh::verttype::NREGULARVERTEX) ||
        (tetMesh_.pointtype(pt) == tetgenmesh::verttype::DEADVERTEX) ) 
        {
            continue;
        }
        if(point_id_map_.find(pt) == point_id_map_.end()) continue;
        // if(abs(pt[0]) <= 1e-3 && abs(pt[1]) < 1e-3 && abs(pt[2]) < 1e-3)
        // {
        //     // printf("zero pt : %f %f %f \n", pt[0], pt[1], pt[2]);
        //     continue;
        // }
        if(pt != (tetgenmesh::point)(NULL)) candid_pts.insert(pt);
    }
}

void VoronoiGen::BuildTetMeshTetCenterMap()
{
    tetMesh_.tetrahedrons->traversalinit();
    tetgenmesh::triface tetface;
    tetface.tet = tetMesh_.alltetrahedrontraverse();
    size_t tet_count = 0;
    while(tetface.tet != (tetgenmesh::tetrahedron*)(NULL))
    {
        // if(tetMesh_.ishulltet(tetface))
        {
            double cx = 0; 
            double cy = 0;
            double cz = 0; 
            for(size_t i = 0; i < 4; ++i)
            {
                auto cur_pt = (tetgenmesh::point) tetface.tet[4 + i];
                cx += cur_pt[0];
                cy += cur_pt[1];
                cz += cur_pt[2];
            }
            cx /= 3.0;
            cy /= 3.0;
            cz /= 3.0;
            tet_center_pts_.push_back(cx);
            tet_center_pts_.push_back(cy);
            tet_center_pts_.push_back(cz);
            tc_pt_tet_map_[tet_count] = tetface.tet;
            tet_count++;
        }
        tetface.tet = tetMesh_.alltetrahedrontraverse();
    }
    std::string tet_center_path = out_dir_ + filename_ + "/" + "tet_centers";
    writePLYFile(tet_center_path, tet_center_pts_);
}

tetgenmesh::tetrahedron* VoronoiGen::GetClosetTet(double x, double y, double z)
{
    // size_t tet_id = pTree_.SearchNearestPt(x, y, z);
    size_t tet_id = 0;
    return tc_pt_tet_map_[tet_id];
}

void VoronoiGen::BuildPicoTree()
{
    // std::vector<double> tet_center_pts;
    // for(auto pt : tet_center_pts_)
    // {
    //     tet_center_pts.push_back(pt[0]);
    //     tet_center_pts.push_back(pt[1]);
    //     tet_center_pts.push_back(pt[2]);
    // }
    // pTree_.Init(tet_center_pts_);
}

void VoronoiGen::GetVoronoiNeiPts(tetgenmesh::point pt, std::vector<tetgenmesh::point>& candid_pts)
{
    // tetlist_->restart();
    
    tetgenmesh::triface searchtet;
    // auto closet_tet = GetClosetTet(pt[0], pt[1], pt[2]);
    // tetMesh_.decode(*closet_tet, searchtet);
    searchtet.tet = NULL;
    // searchtet.tet = GetClosetTet(pt[0], pt[1], pt[2]);
    // printf(" GetInsertPointBoundry .......... \n");
    // ptlist_->restart();
    // tetMesh_.GetInsertPointBoundry(pt, &search_tet, ptlist_);
    tetgenmesh::insertvertexflags ivf;
    ivf.bowywat = 1; // Use Bowyer-Watson algorithm
    ivf.lawson = 0;
    ivf.validflag = 0;
    // tetgenmesh::triface searchtet;
    // searchtet.tet = NULL;
    ivf.iloc = (int) tetgenmesh::UNKNOWN;
    auto t0 = Clock::now();
    tetgenmesh::face *splitsh = NULL;
    tetgenmesh::face *splitseg = NULL;
    // temp_mesh.setpointtype(newpt, tetgenmesh::UNUSEDVERTEX);
    // printf("start GetCaveBoundryPoint ... \n");
    // tetgenmesh::arraypool *out_cavetetvertlist ;
    // std::vector<tetgenmesh::point> out_cavetetvertlist;
    tetMesh_.GetCaveBoundryPoint(pt, &searchtet, splitsh, splitseg, &ivf, candid_pts);
    // tetMesh_.GetCaveBoundryPointMP(pt, &searchtet, splitsh, splitseg, &ivf, candid_pts);
    // printf(" id out_cavetetvertlist size : %ld \n",  out_cavetetvertlist.size());
    // tetMesh_.GetCaveBoundryPoint(pt, &search_tet, ptlist_);
    // for (int i = 0; i < out_cavetetvertlist->objects; i++) {
    //     // tetgenmesh::point* parypt = (tetgenmesh::point *) tetgenmesh::fastlookup(cavetetvertlist, i);
    //     auto pt = * (tetgenmesh::point *) out_cavetetvertlist-> lookup(i);
    //     if(tetMesh_.pointtype(pt) == tetgenmesh::verttype::UNUSEDVERTEX) continue;
    //     // if ((tetMesh_.pointtype(pt) == tetgenmesh::verttype::UNUSEDVERTEX) ||
    //     // (tetMesh_.pointtype(pt) == tetgenmesh::verttype::DUPLICATEDVERTEX) ||
    //     // (tetMesh_.pointtype(pt) == tetgenmesh::verttype::NREGULARVERTEX) ||
    //     // (tetMesh_.pointtype(pt) == tetgenmesh::verttype::DEADVERTEX) ) 
    //     // {
    //     //     continue;
    //     // }
    //     // if(point_id_map_.find(pt) == point_id_map_.end()) continue;
    //     // if(abs(pt[0]) <= 1e-3 && abs(pt[1]) < 1e-3 && abs(pt[2]) < 1e-3)
    //     // {
    //     //     // printf("zero pt : %f %f %f \n", pt[0], pt[1], pt[2]);
    //     //     continue;
    //     // }
    //     if(pt != (tetgenmesh::point)(NULL)) candid_pts.insert(pt);
    // }
    // delete out_cavetetvertlist;
}

void VoronoiGen::InsertPt(tetgenmesh::point pt)
{
    tetgenmesh::point newpt;
    // auto temp_mesh = tetgenmesh(tetMesh_);
    auto& temp_mesh = tetMesh_;
    printf("deep copy succed! \n");

    temp_mesh.makepoint(&newpt, tetgenmesh::UNUSEDVERTEX);
    newpt[0] = pt[0];
    newpt[1] = pt[1];
    newpt[2] = pt[2];
    tetgenmesh::insertvertexflags ivf;
    ivf.bowywat = 1; // Use Bowyer-Watson algorithm
    ivf.lawson = 0;

    tetgenmesh::triface searchtet;
    searchtet.tet = NULL;
  
    ivf.iloc = (int) tetgenmesh::UNKNOWN;
    auto t0 = Clock::now();
    tetgenmesh::face *splitsh = NULL;
    tetgenmesh::face *splitseg = NULL;
    temp_mesh.setpointtype(newpt, tetgenmesh::VOLVERTEX);

    printf("-----------init tet num: %lu ...\n", tetMesh_.tetrahedrons->items);
    if(temp_mesh.insertpoint(newpt, &searchtet, splitsh, splitseg,  &ivf))
    // if (tetMesh_.insert_vertex_bw(newpt, &searchtet, &ivf))
    {
        std::set<tetgenmesh::point> candidate_pts;
        GetVertexStar(newpt, candidate_pts, 1);
        printf("--------- insert pt  nei points num : %lu \n", candidate_pts.size());
        // temp_mesh.removevertexbyflips();
        // printf("insert successfully ...\n");
        printf("current pt type: %d \n", tetMesh_.pointtype(newpt));
        // printf("current pt num: %lu ...\n", tetMesh_.points->items);
        printf("      Removing a point %d in volume.\n",
             tetMesh_.pointmark(newpt));
             
        printf("-----------insert tet num: %lu ...\n", tetMesh_.tetrahedrons->items);
        
        // tetMesh_.points->dealloc(newpt);
        temp_mesh.b->fliplinklevel = 100000;
        temp_mesh.removevertexbyflips(newpt);
        // temp_mesh.removevertexbyflips(newpt);
        // temp_mesh.removevertexbyflips(newpt);
        candidate_pts.clear();
        GetVertexStar(newpt, candidate_pts, 1);
        printf("after remove insert pt  nei points num : %lu \n", candidate_pts.size());
        // tetMesh_.ed
        // tetMesh_.points->dealloc(newpt);

        // temp_mesh.removevertexbyflips(newpt);
        printf("----------remove tet num : %lu ...\n", tetMesh_.tetrahedrons->items);
        printf("current pt num after remove: %lu ...\n", tetMesh_.points->items);
        // tetMesh_.removevertexbyflips();
        // printf("remove successfully ...\n");
        //  printf("current pt num: %llu ...\n", tetMesh_.points->items);
    }
    auto t1 = Clock::now();
    // double insert_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    // printf(" point insert time : %g \n", insert_time);
    


    // printf("start to insert ........ \n");
    // tetMesh_.insertconstrainedpoints(insertarray, arylen, rejflag);
}

void VoronoiGen::CalcualteNormalWithVIPSS(std::vector<double>& vts, std::vector<double>& normal)
{
    // RBF_API vipss_api;
    vipss_api_.Set_RBF_PARA();
    vipss_api_.run_vipss(vts);

}

std::vector<uint8_t> cal_colors(size_t ptn)
{
    std::vector<uint8_t> colors;
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

void VoronoiGen::BuildPtIdMap()
{
    tetMesh_.points->traversalinit();
    tetgenmesh::point ploop = tetMesh_.pointtraverse();
    size_t pt_num_ = 0;
    points_.clear();
    while(ploop != (tetgenmesh::point)NULL)
    {
        const auto p_type = tetMesh_.pointtype(ploop);
        if ((p_type == tetgenmesh::UNUSEDVERTEX) ) 
        {
            ploop = tetMesh_.pointtraverse();
            continue;
        } 
        points_.push_back(ploop);
        point_id_map_[ploop] = pt_num_;
        ploop = tetMesh_.pointtraverse();
        pt_num_ ++;
    }
    pt_score_mat_.resize(pt_num_, pt_num_);
    pt_adjecent_mat_.resize(pt_num_, pt_num_);
    pt_dist_mat_.resize(pt_num_, pt_num_);
    printf("pt num : %zu \n", pt_num_);
    printf("point_id_map_ size : %zu \n", point_id_map_.size());

}

void VoronoiGen::BuildAdjecentMat()
{
    tetMesh_.points->traversalinit();
    tetgenmesh::point ploop = tetMesh_.pointtraverse();
    // points_.clear()
    while(ploop != (tetgenmesh::point)NULL)
    {
        // printf("new cur pt id 0000001\n");
        if(tetMesh_.pointtype(ploop) == tetgenmesh::UNUSEDVERTEX )
        {
            ploop = tetMesh_.pointtraverse();
            continue;
        }
        // printf("new cur pt id 0000002\n");
        // while(point_id_map_.find(ploop) == point_id_map_.end() && ploop != (tetgenmesh::point)NULL)
        // {
        //     ploop = tetMesh_.pointtraverse();
        // }
        // printf("new cur pt id 0000003\n");
        // points_.push_back(ploop);
        size_t cur_p_id = point_id_map_[ploop]; 
        // printf("cur pt id : %d \n", cur_p_id);
        std::set<tetgenmesh::point> candidate_pts;
        GetVertexStar(ploop, candidate_pts, 1);
        point_cluster_pts_map_[ploop] = candidate_pts;
        for(auto &pt : candidate_pts)
        {
            if(pt == ploop) continue;
            if(point_id_map_.find(pt) == point_id_map_.end()) continue;
            size_t p_id = point_id_map_[pt];
            pt_adjecent_mat_(cur_p_id, p_id) = 1;
        }
        ploop = tetMesh_.pointtraverse();
        // printf("cur pt id 000 : %d \n", cur_p_id);
    }
    // printf("finish loop pt id 000\n");
}

std::vector<double> VoronoiGen::ConvertPtAndNeighborsToVect( std::set<tetgenmesh::point>& candid_pts)
{
    std::vector<double> vts;
    for(auto cur_pt : candid_pts)
    {
        vts.push_back(cur_pt[0]);
        vts.push_back(cur_pt[1]);
        vts.push_back(cur_pt[2]);
    }
    return vts;
}

void VoronoiGen::BuildPtNCluster(P_Set& pset, const std::vector<double>& normals, PtNCluster& pt_cluster)
{
    size_t id = 0;
    if(pset.size() != normals.size() /3) return;
    
    for(auto& pt: pset)
    {
        double nx = normals[3*id];
        double ny = normals[3*id + 1];
        double nz = normals[3*id + 2];
        arma::vec3 cur_n{nx, ny, nz};
        pt_cluster[pt] = cur_n;
        id++;
    }
}


void VoronoiGen::EstimateNormals()
{
    printf("start to estimate normals \n");
    vipss_api_.Set_RBF_PARA();
    vipss_api_.is_surfacing_ = false;

    clock_t tv[7];
    tetMesh_.points->traversalinit();
    tetgenmesh::point ploop = tetMesh_.pointtraverse();
    tv[1] = clock();
    auto t1 = Clock::now();
    size_t id = 0;
    while(ploop != (tetgenmesh::point)NULL)
    {
        candidate_pts_.clear();

        // if ((tetMesh_.pointtype(ploop) == tetgenmesh::verttype::UNUSEDVERTEX) ||
        // (tetMesh_.pointtype(ploop) == tetgenmesh::verttype::DUPLICATEDVERTEX) ||
        // (tetMesh_.pointtype(ploop) == tetgenmesh::verttype::NREGULARVERTEX) ||
        // (tetMesh_.pointtype(ploop) == tetgenmesh::verttype::DEADVERTEX) ||
        // (tetMesh_.pointtype(ploop) == tetgenmesh::verttype::FREEFACETVERTEX) ||
        // (tetMesh_.pointtype(ploop) == tetgenmesh::verttype::VOLVERTEX)) 
        // {
        //     ploop = tetMesh_.pointtraverse();
        //     continue;
        // }

        // printf("%d point type :  %d \n", id, tetMesh_.pointtype(ploop));
        // id ++;
        GetVertexStar(ploop, candidate_pts_, 1);
        // for(auto& new_pt : candidate_pts_)
        // {
        //     if(abs(new_pt[0]) <= 1e-3 && abs(new_pt[1]) < 1e-3 && abs(new_pt[2]) < 1e-3)
        //     {
        //         printf("zero pt : %f %f %f \n", new_pt[0], new_pt[1], new_pt[2]);
        //     }
        // }
        PtCluster new_cluster;
        new_cluster.key_pts_.insert(ploop);
        new_cluster.group_pts_ = candidate_pts_;

        
        point_cluster_map_[ploop] = new_cluster;
        auto vts = ConvertPtAndNeighborsToVect(candidate_pts_);
        vipss_api_.run_vipss(vts);
        
        PtNCluster cur_cluster;
        BuildPtNCluster(candidate_pts_, vipss_api_.normals_, cur_cluster);

        // for(auto new_pt : new_cluster.key_pts)
        // {
        //     printf(" cluster key pt : %f %f %f \n", new_pt[0], new_pt[1], new_pt[2]);
        //     if(cur_cluster.find(new_pt) != cur_cluster.end())
        //     {
        //         printf("group_pts pt : %f %f %f \n", new_pt[0], new_pt[1], new_pt[2]);
        //     }
        // }

        arma::vec3 new_normal;
        if(cur_cluster.find(ploop) != cur_cluster.end()) 
        {
            new_normal = cur_cluster[ploop];
        }
        point_cluster_normal_map_[ploop] = cur_cluster;
        pt_normal_map_[ploop] = new_normal;

        ploop = tetMesh_.pointtraverse();
    }
    tv[2] = clock();
    
    auto t2 = Clock::now();
    double estimate_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
    printf("normal estimation time : %g \n", estimate_time);

    // printf("pt normal map size :  %d \n", pt_normal_map_.size());

}

void VoronoiGen::CalculateScores()
{
    pt_score_mat_.zeros();
    for(auto& ele : point_cluster_normal_map_)
    {
        auto cur_pt = ele.first;
        auto& cluster = ele.second;
        const auto& cur_n = cluster[cur_pt];
        size_t cur_id = point_id_map_[cur_pt];
        for(const auto& pt_ele : cluster)
        {   
            auto& s_pt = pt_ele.first;
            if(cur_pt == s_pt) continue;
            if(pt_normal_map_.find(s_pt) == pt_normal_map_.end()) continue;
            size_t s_id = point_id_map_[s_pt];
            if(pt_score_mat_(cur_id, s_id) != 0.0) continue;
            auto& s_n = pt_ele.second;
            auto& ori_n = pt_normal_map_[s_pt];
            double val = std::min(1.0, abs(arma::dot(ori_n, s_n)));
            double angle = acos (val) * 180.0 / M_PI ;
            if(point_cluster_normal_map_.find(s_pt) == point_cluster_normal_map_.end()) continue;

            if(point_cluster_normal_map_[s_pt].find(cur_pt) == point_cluster_normal_map_[s_pt].end())
            {
                continue;
            }   

            auto cur_n_s = point_cluster_normal_map_[s_pt][cur_pt];
            double val_re = std::min(1.0, abs(arma::dot(cur_n_s, cur_n)));
            double angle_re = acos (val_re) * 180.0 / M_PI ;
            double max_angle = std::max(angle, angle_re);
            // printf("val : %f, val re: %f \n", val, val_re);

            pt_score_mat_(cur_id, s_id) = max_angle;
            pt_score_mat_(s_id, cur_id) = max_angle;
            
            double cur_dist = PtDistance(cur_pt, s_pt);
            cur_dist = std::max(1e-8, cur_dist);
            pt_dist_mat_(cur_id, s_id) = cur_dist;
            pt_dist_mat_(s_id, cur_id) = cur_dist;
        }
    }

    size_t ptn = point_id_map_.size();
    for(auto& ele : point_id_map_)
    {
        auto& cur_pt = ele.first;
        size_t i = ele.second;
        double dist_sum = arma::sum(pt_dist_mat_.row(i));
        // printf("dist sum : %f \n",  dist_sum);
        pt_dist_mat_.row(i) /= dist_sum;
        double cur_score = arma::sum(pt_score_mat_.row(i) % pt_dist_mat_.row(i)) ;
        // arma::vec nozeros = arma::nonzeros(pt_score_mat_.row(i));
        // double cur_score = arma::mean(nozeros);
        // printf("score : %f \n",  arma::sum(pt_score_mat_.row(i))); 
        pt_score_map_[cur_pt] = cur_score;
    }

}

void VoronoiGen::OrientPtNormals()
{
    // std::vector<double> pts;
    // std::vector<double> ptns;
    // for(auto &ele : pt_normal_map_)
    // {
    //     // printf("point : %f %f %f \n", ele.first[0], ele.first[1], ele.first[2]);
    //     // printf("normal : %f %f %f \n", ele.second[0], ele.second[1], ele.second[2]);
    //     if(ele.first == (tetgenmesh::point)NULL) continue;
    //     pts.push_back(ele.first[0]);
    //     pts.push_back(ele.first[1]);
    //     pts.push_back(ele.first[2]);

    //     ptns.push_back(ele.second[0]);
    //     ptns.push_back(ele.second[1]);
    //     ptns.push_back(ele.second[2]);
        
    // }
    // ORIENT::OrientPointNormals(pts, ptns);
}

void VoronoiGen::InitMesh()
{
    
    Tetrahedralize();
    // BuildPtIdMap();
    
    
    tetlist_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::triface), 10);
    ptlist_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::point), 10);

    tetlist2_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::triface), 10);
    ptlist2_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::point), 10);

    tetlist3_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::triface), 10);
    ptlist3_ = new tetgenmesh::arraypool(sizeof(tetgenmesh::point), 10);
}


void VoronoiGen::UpdateVtAndVn()
{
    normals_.clear();
    vertices_.clear();
    // printf("pt normal map size :  %d \n", pt_normal_map_.size());
    for(auto &ele : pt_normal_map_)
    {
        if(ele.first == (tetgenmesh::point)NULL) continue;
        // if(abs(ele.first[0]) < 1e-12 && abs(ele.first[1]) < 1e-12 && abs(ele.first[2]) < 1e-12)
        // {
        //     continue;
        // }
        vertices_.push_back(ele.first[0]);
        vertices_.push_back(ele.first[1]);
        vertices_.push_back(ele.first[2]);

        normals_.push_back(ele.second[0]);
        normals_.push_back(ele.second[1]);
        normals_.push_back(ele.second[2]);
    }
    // ORIENT::OrientPointNormals(vertices_, normals_);
}

void VoronoiGen::SavePtVn(const std::string& path, bool orient_normal)
{
    
    // std::string out_normal_path = "data/vipss_estimated_ns";
    // if(orient_normal)
    // {
    //     ORIENT::OrientPointNormals(vertices_, normals_);
    // }
    writePLYFile_VN(path, vertices_, normals_);
    printf("successfully save points and normals to file : %s", path.c_str());
}


void VoronoiGen::SavePtVnColor(const std::string& path, bool orient_normal)
{
    // std::string out_normal_path = "data/vipss_estimated_ns";
    // if(orient_normal)
    // {
    //     ORIENT::OrientPointNormals(vertices_, normals_);
    // }
    writePLYFile_VN_CO(path, vertices_, normals_, colors_);
    printf("successfully save points and normals to file : %s \n", path.c_str());
}


void VoronoiGen::ScaleNormalByScore()
{

    double max_score = 1e-8;
    // double min_score = 180;
    for(auto& ele : pt_score_map_)
    {
        max_score = max_score > ele.second ? max_score : ele.second;
        // min_score = min_score < ele.second ? min_score : ele.second;
    }
    normals_.clear();
    for(auto &ele : pt_normal_map_)
    {
        if(ele.first == (tetgenmesh::point)NULL) continue;
        double ratio = pt_score_map_[ele.first] / max_score;
        normals_.push_back(ele.second[0] * ratio);
        normals_.push_back(ele.second[1] * ratio);
        normals_.push_back(ele.second[2] * ratio);
    }
}

typedef std::pair<tetgenmesh::point, double> PtScore;

void VoronoiGen::MergeCluster()
{
    std::vector<PtScore> pt_score_pairs;
    for(auto& ele : pt_score_map_)
    {
        PtScore cur_ptscore(ele.first, ele.second);
        pt_score_pairs.push_back(cur_ptscore);
    }
    std::sort(pt_score_pairs.begin(), pt_score_pairs.end(), [](PtScore& a, PtScore& b){ return a.second > b.second;});

    std::set<tetgenmesh::point> visited_pts;
    for(auto & ele: pt_score_pairs)
    {
        // printf(" pt score : %f \n", ele.second );
        if(visited_pts.find(ele.first) != visited_pts.end()) 
        {
            continue;
        }
        visited_pts.insert(ele.first);
        if(ele.second < min_angle_)
        {
            continue;
        }
        if(point_cluster_map_.find(ele.first) != point_cluster_map_.end())
        {
            auto &cluster = point_cluster_map_[ele.first];
            tetgenmesh::point max_pt = tetgenmesh::point(NULL);
            double max_score = 0;
            for(auto& pt : cluster.group_pts_)
            {
                if(pt == ele.first) continue;
                if(pt_score_map_.find(pt) != pt_score_map_.end())
                {
                    if( pt_score_map_[pt] > max_score)
                    {
                        max_score = pt_score_map_[pt];
                        max_pt = pt;
                    }
                }
            }
            if(max_pt != tetgenmesh::point (NULL))
            {
                visited_pts.insert(max_pt);
                auto &max_cluster = point_cluster_map_[max_pt];
                cluster.group_pts_.insert(max_cluster.group_pts_.begin(), max_cluster.group_pts_.end());
                cluster.key_pts_.insert(max_cluster.key_pts_.begin(), max_cluster.key_pts_.end());
                point_cluster_map_[max_pt] = cluster;

                auto vts = ConvertPtAndNeighborsToVect(cluster.group_pts_);
                // printf("vipss vts size : %d \n", vts.size()/3);
                vipss_api_.run_vipss(vts);
   
                // std::string out_path = out_dir_ + filename_ + "_cluster";
                // writePLYFile_VN(out_path, vts, vipss_api_.normals_);

                // printf("vipss normal size : %d \n", vipss_api_.normals_.size()/3);
                // printf("cur group_pts size : %d \n", cluster.group_pts.size());
                
                PtNCluster cur_cluster;
                BuildPtNCluster(cluster.group_pts_, vipss_api_.normals_, cur_cluster);

                // printf("cur cluster size : %d \n", cur_cluster.size());
                // for(auto& k_pt : cluster.group_pts)
                // {
                //     printf("key pt : %f %f %f \n", k_pt[0], k_pt[1], k_pt[2]);
                //     if(cur_cluster.find(k_pt) != cur_cluster.end())
                //     {
                //         printf("contianed pt : %f %f %f \n", k_pt[0], k_pt[1], k_pt[2]);
                //     }
                // }
                for(auto k_pt : cluster.key_pts_)
                {
                    arma::vec3 new_normal;
                    if(cur_cluster.find(k_pt) != cur_cluster.end()) 
                    {
                        new_normal = cur_cluster[k_pt];
                        // printf("new normal : %f %f %f \n", new_normal[0], new_normal[1], new_normal[2]);

                        if(pt_normal_map_.find(k_pt) != pt_normal_map_.end())
                        {
                            pt_normal_map_[k_pt] = new_normal;
                        } 
                    }
                    if(point_cluster_normal_map_.find(k_pt) != point_cluster_normal_map_.end())
                    {
                        point_cluster_normal_map_[k_pt] = cur_cluster;
                    }
                }
            }
        }
    }
}




void VoronoiGen::CalculatePtColors()
{
    colors_.clear();
    double max_score = 0;
    // double min_score = 180;
    for(auto& ele : pt_score_map_)
    {
        max_score = max_score > ele.second ? max_score : ele.second;
        // min_score = min_score < ele.second ? min_score : ele.second;
    }
    max_score = 2.0 * min_angle_; 
    for(auto &ele : pt_normal_map_)
    {

        auto cur_pt = ele.first;
        auto score = pt_score_map_[cur_pt];
        score = std::min(score, max_score);

        uint8_t c_r = uint8_t(score / max_score * 255);
        uint8_t c_g = 0;
        uint8_t c_b = uint8_t( (max_score - score) / max_score * 255); 

        colors_.push_back(c_r);
        colors_.push_back(c_g);
        colors_.push_back(c_b);
    }
}




void VoronoiGen::Run()
{
    InitMesh();
    EstimateNormals();
    CalculateScores();
    UpdateVtAndVn();
    CalculatePtColors();
    std::string out_normal_path = out_dir_ + filename_ + "_origin_ns";
    SavePtVn(out_normal_path, false);
 
    std::string pt_color_path =  out_dir_ +  filename_ + "_color";
    SavePtVnColor(pt_color_path,  false);
    ScaleNormalByScore();
    std::string pt_line_path = out_dir_ +  filename_ + "_line";
    writeObjPtn_line(pt_line_path, vertices_, colors_, normals_);

    MergeCluster();
    CalculateScores();
    UpdateVtAndVn();
    CalculatePtColors();
    std::string out_normal_path1 = out_dir_ + filename_ + "_origin_ns1";
    SavePtVn(out_normal_path1, false);

    std::string pt_color_path1 =  out_dir_ +  filename_ + "_color1";
    SavePtVnColor(pt_color_path1,  false);
    ScaleNormalByScore();
    std::string pt_line_path1 = out_dir_ +  filename_ + "_line1";
    writeObjPtn_line(pt_line_path1, vertices_, colors_, normals_);


}


std::unordered_map<tetgenmesh::point, cluster::PtCluster> VoronoiGen::BuildAllClusters()
{
    tetMesh_.points->traversalinit();
    tetgenmesh::point ploop = tetMesh_.pointtraverse();
    auto t1 = Clock::now();
    size_t id = 0;
    std::unordered_map<tetgenmesh::point, PtCluster> pt_cluster_map;
    while(ploop != (tetgenmesh::point)NULL)
    {
        candidate_pts_.clear();
        GetVertexStar(ploop, candidate_pts_, 1);
        if(candidate_pts_.find(ploop) != candidate_pts_.end()) candidate_pts_.erase(ploop);
        cluster::PtCluster new_cluster;
        new_cluster.key_pts_.insert(ploop);
        new_cluster.group_pts_ = candidate_pts_;
        pt_cluster_map[ploop] = new_cluster;
        ploop = tetMesh_.pointtraverse();
    }
    return pt_cluster_map;
}