#include "voronoi_gen.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "readers.h"
#include "orient_normal.h"
#include <chrono>
#include <math.h>
#include <cmath>

typedef std::chrono::high_resolution_clock Clock;

REAL cps = (REAL) CLOCKS_PER_SEC;
double M_PI  = 2*acos(0.0);

inline double PtDistance(double* p1, double* p2)
{
  return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
              (p2[1] - p1[1]) * (p2[1] - p1[1]) +
              (p2[2] - p1[2]) * (p2[2] - p1[2]));
}

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
        // if(abs(pt[0]) <= 1e-3 && abs(pt[1]) < 1e-3 && abs(pt[2]) < 1e-3)
        // {
        //     // printf("zero pt : %f %f %f \n", pt[0], pt[1], pt[2]);
        //     continue;
        // }
        if(pt != (tetgenmesh::point)(NULL)) candid_pts.insert(pt);
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
    size_t count = 0;
    while(ploop != (tetgenmesh::point)NULL)
    {
        point_id_map_[ploop] = count;
        ploop = tetMesh_.pointtraverse();
        count ++;
    }
    pt_score_mat_.resize(count, count);
    pt_adjecent_mat_.resize(count, count);
    pt_dist_mat_.resize(count, count);
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
        new_cluster.key_pts.insert(ploop);
        new_cluster.group_pts = candidate_pts_;

        
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
            double val = min(1.0, abs(arma::dot(ori_n, s_n)));
            double angle = acos (val) * 180.0 / M_PI ;
            if(point_cluster_normal_map_.find(s_pt) == point_cluster_normal_map_.end()) continue;

            if(point_cluster_normal_map_[s_pt].find(cur_pt) == point_cluster_normal_map_[s_pt].end())
            {
                continue;
            }   

            auto cur_n_s = point_cluster_normal_map_[s_pt][cur_pt];
            double val_re = min(1.0, abs(arma::dot(cur_n_s, cur_n)));
            double angle_re = acos (val_re) * 180.0 / M_PI ;
            double max_angle = max(angle, angle_re);
            // printf("val : %f, val re: %f \n", val, val_re);

            pt_score_mat_(cur_id, s_id) = max_angle;
            pt_score_mat_(s_id, cur_id) = max_angle;
            
            double cur_dist = PtDistance(cur_pt, s_pt);
            cur_dist = max(1e-8, cur_dist);
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
    BuildPtIdMap();
    
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
            for(auto& pt : cluster.group_pts)
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
                cluster.group_pts.insert(max_cluster.group_pts.begin(), max_cluster.group_pts.end());
                cluster.key_pts.insert(max_cluster.key_pts.begin(), max_cluster.key_pts.end());
                point_cluster_map_[max_pt] = cluster;

                auto vts = ConvertPtAndNeighborsToVect(cluster.group_pts);
                // printf("vipss vts size : %d \n", vts.size()/3);
                vipss_api_.run_vipss(vts);
   
                // std::string out_path = out_dir_ + filename_ + "_cluster";
                // writePLYFile_VN(out_path, vts, vipss_api_.normals_);

                // printf("vipss normal size : %d \n", vipss_api_.normals_.size()/3);
                // printf("cur group_pts size : %d \n", cluster.group_pts.size());
                
                PtNCluster cur_cluster;
                BuildPtNCluster(cluster.group_pts, vipss_api_.normals_, cur_cluster);

                // printf("cur cluster size : %d \n", cur_cluster.size());
                // for(auto& k_pt : cluster.group_pts)
                // {
                //     printf("key pt : %f %f %f \n", k_pt[0], k_pt[1], k_pt[2]);
                //     if(cur_cluster.find(k_pt) != cur_cluster.end())
                //     {
                //         printf("contianed pt : %f %f %f \n", k_pt[0], k_pt[1], k_pt[2]);
                //     }
                // }
                for(auto k_pt : cluster.key_pts)
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
        // if(abs(ele.first[0]) < 1e-12 && abs(ele.first[1]) < 1e-12 && abs(ele.first[2]) < 1e-12)
        // {
        //     continue;
        // }
        auto cur_pt = ele.first;
        auto score = pt_score_map_[cur_pt];
        score = min(score, max_score);

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