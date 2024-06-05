#include "local_vipss.hpp"
#include <cmath>
#include <algorithm>
#include "readers.h"  
#include "orient_normal.h"
#include <chrono>
#include <queue>

typedef std::chrono::high_resolution_clock Clock;

void LocalVipss::Init(const std::string & path)
{   
    std::vector<double> in_pts;
    std::vector<double> in_normals;
    readPLYFile(path, in_pts, in_normals);
    printf("load data file : %s \n", path.c_str());

    printf("read point size : %zu \n", in_pts.size());
    auto t0 = Clock::now();
    // printf("read point size : %d \n", in_pts.size());
    // voro_gen_.loadData(path);
    voro_gen_.loadData(in_pts);
    voro_gen_.InitMesh();
    voro_gen_.BuildPtIdMap();
    voro_gen_.BuildAdjecentMat();

    adjacent_mat_ = voro_gen_.pt_adjecent_mat_;
    points_ = voro_gen_.points_;

    printf("adjacent mat rows : %lld, cols : %lld \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    pt_num_ = points_.size();

    cluster_cores_mat_.resize(pt_num_, pt_num_);
    cluster_cores_mat_.eye();
    
    printf("input point size : %zu \n", pt_num_);

    vipss_api_.Set_RBF_PARA();
    vipss_api_.is_surfacing_ = false;
    vipss_api_.outpath_ = out_dir_;
    vipss_api_.user_lambda_ = user_lambda_;

    cluster_normal_x_.resize(pt_num_, pt_num_);
    cluster_normal_y_.resize(pt_num_, pt_num_);
    cluster_normal_z_.resize(pt_num_, pt_num_);

    cluster_scores_mat_.resize(pt_num_, pt_num_);
    cluster_adjacent_mat_.resize(pt_num_, pt_num_);
    auto t1 = Clock::now(); 
    double init_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("finish local vipss initilization : %f ! \n", init_time);
}

inline std::vector<size_t> LocalVipss::GetClusterCoreIds(size_t cluster_id) const
{
    const arma::sp_irowvec& cluster_core_ids = cluster_cores_mat_.row(cluster_id);
    arma::sp_irowvec::const_iterator start = cluster_core_ids.begin();
    arma::sp_irowvec::const_iterator end = cluster_core_ids.end();
    std::vector<size_t> core_ids;
    for ( arma::sp_irowvec::const_iterator i = start; i != end; ++i )
    {
        core_ids.push_back(i.internal_col);
    }
    return core_ids;
}


void LocalVipss::BuidClusterCoresPtIds()
{
    cluster_core_pt_ids_.clear();
    size_t c_num = cluster_cores_mat_.n_rows;
    for(size_t i =0; i < c_num; ++i)
    {
        auto cur_cores = GetClusterCoreIds(i);
        cluster_core_pt_ids_.push_back(cur_cores);
    }
}

void LocalVipss::UpdateClusterCoresPtIds()
{
    size_t c_num = cluster_cores_mat_.n_rows;
    size_t cur_n = cluster_core_pt_ids_.size();
    if(c_num <= cur_n) return; 
    for(size_t i = cur_n; i < c_num; ++i)
    {
        auto p_ids = GetClusterCoreIds(i);
        cluster_core_pt_ids_.push_back(p_ids);
    }
}

inline void LocalVipss::GetInitClusterPtIds(size_t cluster_id, 
        std::vector<double>& pts, std::vector<size_t>& pt_ids)
{
    auto& cluster_pts_map = voro_gen_.point_cluster_pts_map_;
    auto& cluster_pt_id_map = voro_gen_.point_id_map_;
    // auto& points = voro_gen_.points_;
    // printf("cluster id %d, points size %d \n ", cluster_id, points_.size()/3);
    tetgenmesh::point cur_ptr = points_[cluster_id];
    auto cluster_pts = voro_gen_.point_cluster_pts_map_[cur_ptr];

    pt_ids.push_back(cluster_id);
    pts.push_back(cur_ptr[0]);
    pts.push_back(cur_ptr[1]);
    pts.push_back(cur_ptr[2]);

    for(auto& ptr: cluster_pts)
    {
        if(ptr == cur_ptr) continue;
        if(cluster_pt_id_map.find(ptr) == cluster_pt_id_map.end()) continue;
        // if()
        // std::cout << " pt addr : " << ptr << std::endl;
        size_t c_id = cluster_pt_id_map[ptr];
        // printf("cur pt id %d \n", c_id);
        pt_ids.push_back(c_id);
        // printf("cur_ptr %f, %f, %f \n ", ptr[0], ptr[1], ptr[2]);
        pts.push_back(ptr[0]);
        pts.push_back(ptr[1]);
        pts.push_back(ptr[2]);
    }
}

inline std::vector<size_t> LocalVipss::GetClusterPtIds(size_t cluster_id) const
{
    // arma::sp_ivec cluster_core_ids = cluster_cores_mat_.row(cluster_id);
    // arma::sp_ivec cluster_pt_ids = adjacent_mat_ * cluster_core_ids;
    std::vector<size_t> pt_ids;
    std::set<size_t> key_ids;
    const arma::sp_irowvec& cluster_core_ids = cluster_cores_mat_.row(cluster_id);
    arma::sp_irowvec::const_iterator c_start = cluster_core_ids.begin();
    arma::sp_irowvec::const_iterator c_end = cluster_core_ids.end();
    for (auto i = c_start; i != c_end; ++i)
    {
        pt_ids.push_back(i.internal_col);
        key_ids.insert(i.internal_col);
    }
    const arma::sp_irowvec& cluster_pt_ids = cluster_adjacent_pt_mat_.row(cluster_id);
    arma::sp_irowvec::const_iterator start = cluster_pt_ids.begin();
    arma::sp_irowvec::const_iterator end = cluster_pt_ids.end();
    for ( arma::sp_irowvec::const_iterator i = start; i != end; ++i )
    {
        if(key_ids.find(i.internal_col) != key_ids.end()) continue;
        pt_ids.push_back(i.internal_col);
    }
    return pt_ids;
}

inline std::vector<double> LocalVipss::GetClusterVerticesFromIds(const std::vector<size_t>& pt_ids) const
{
    std::vector<double> vertices;
    for(size_t id : pt_ids)
    {
        auto &pt = points_[id];
        vertices.push_back(pt[0]);
        vertices.push_back(pt[1]);
        vertices.push_back(pt[2]);
    }
    return vertices;
}

std::vector<double> LocalVipss::GetClusterVertices(size_t cluster_id) const
{
    auto cluster_pt_ids = GetClusterPtIds(cluster_id);
    return GetClusterVerticesFromIds(cluster_pt_ids);
}


size_t LocalVipss::GetClusterIdFromCorePtId(const size_t pid)
{
    const arma::sp_icolvec& cluster_col = cluster_cores_mat_.col(pid);
    arma::sp_icolvec::const_iterator start = cluster_col.begin();
    return (size_t)start.internal_pos;
}

void LocalVipss::CalculateClusterNormals(size_t cluster_id)
{
    auto cluster_pt_ids = GetClusterPtIds(cluster_id);

    // printf("cluster_pt_ids size : %zu \n", cluster_pt_ids.size());
    auto cluster_vts = GetClusterVerticesFromIds(cluster_pt_ids);
    // printf("cluster pt size : % \n", cluster_vts.size()/3);
    auto t1 = Clock::now();
    // vipss_api_.run_vipss(cluster_vts, cluster_pt_ids.size());
    const arma::sp_irowvec& cluster_core_ids = cluster_cores_mat_.row(cluster_id);
    size_t key_ptn = cluster_core_ids.n_nonzero;
    vipss_api_.run_vipss(cluster_vts, key_ptn);
    auto t2 = Clock::now();
    double vipss_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
    vipss_time_stats_.emplace_back(std::pair<size_t, double>(cluster_pt_ids.size(), vipss_time));

    for(size_t p_id = 0; p_id < cluster_pt_ids.size(); ++p_id)
    {
        size_t col = cluster_pt_ids[p_id];
        cluster_normal_x_(cluster_id, col) = vipss_api_.normals_[3* p_id];
        cluster_normal_y_(cluster_id, col) = vipss_api_.normals_[3* p_id + 1];
        cluster_normal_z_(cluster_id, col) = vipss_api_.normals_[3* p_id + 2];
    }
}


void LocalVipss::InitNormalWithVipss()
{
    size_t cluster_num = cluster_cores_mat_.n_rows;
    printf("cluster num : %zu \n", cluster_num);
    
    vipss_time_stats_.clear();
    double vipss_sum = 0;
    auto t_init = Clock::now();

    for(size_t i =0; i < cluster_num; ++i)
    {
        std::vector<double> vts;
        std::vector<size_t> p_ids;
        GetInitClusterPtIds(i, vts, p_ids);
        // printf("vt size : %d , p ids : %d \n", vts.size(), p_ids.size());
        auto t1 = Clock::now();
        vipss_api_.run_vipss(vts, 1);
        auto t2 = Clock::now();
        double vipss_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
        vipss_sum += vipss_time;
        vipss_time_stats_.emplace_back(std::pair<size_t, double>(p_ids.size(), vipss_time));

        for(size_t p_id = 0; p_id < p_ids.size(); ++p_id)
        {
            size_t col = p_ids[p_id];
            cluster_normal_x_(i, col) = vipss_api_.normals_[3* p_id];
            cluster_normal_y_(i, col) = vipss_api_.normals_[3* p_id + 1];
            cluster_normal_z_(i, col) = vipss_api_.normals_[3* p_id + 2];
        }
    }
    auto t_init2 = Clock::now();
    double all_init_time = std::chrono::nanoseconds(t_init2 - t_init).count()/1e9;
    printf("all init time : %f \n", all_init_time);

    printf("all init vipss time : %f \n", vipss_sum);

    cluster_ptn_vipss_time_stats_.push_back(vipss_time_stats_);
}



inline double LocalVipss::CalculateScores(const arma::mat& a_normals, const arma::mat& b_normals) const
{
    arma::mat dot_mat = a_normals % b_normals;
    arma::colvec dot_sum = arma::sum(dot_mat, 1);
    double min_project_p = dot_sum.min();
    double min_project_n = -1.0 * dot_sum.max();

    min_project_p = min_project_n > min_project_p ? min_project_n : min_project_p;

    // flip_normal_ = false;
    // if(min_project_n > min_project_p)
    // {
    //     // flip_normal_ = true;
    //     min_project_p = min_project_n;
    // }
    min_project_p = std::min(1.0, std::max(-1.0, min_project_p));
    double angle = acos (min_project_p) * Anlge_PI_Rate ;
    return angle;
}

inline bool LocalVipss::IsFlipNormal(const arma::mat& a_normals, const arma::mat& b_normals) const
{
    arma::mat dot_mat = a_normals % b_normals;
    arma::colvec dot_sum = arma::sum(dot_mat, 1);
    double min_project_p = dot_sum.min();
    double min_project_n = -1.0 * dot_sum.max();
    return min_project_n > min_project_p? true : false;

    // flip_normal_ = false;
    // if(min_project_n > min_project_p)
    // {
    //     flip_normal_ = true;
    //     min_project_p = min_project_n;
    // }
    // min_project_p = std::min(1.0, std::max(-1.0, min_project_p));
    // double angle = acos (min_project_p) * 180.0 / M_PI2 ;
    // return angle;
}

inline double LocalVipss::CalculateScores(const std::vector<arma::vec3>& a_normals, const std::vector<arma::vec3>& b_normals) const
{

    double min_project_p = 1.0;
    for(size_t i = 0; i < a_normals.size(); ++i)
    {
        double ab_proj = arma::dot(a_normals[i], b_normals[i]);
        min_project_p = min_project_p < ab_proj ? min_project_p : ab_proj;
    }
    
    double min_project_n = 1.0;
    for(size_t i = 0; i < a_normals.size(); ++i)
    {
        double ab_proj = arma::dot(a_normals[i], - 1.0 * b_normals[i]);
        min_project_n = min_project_n < ab_proj ? min_project_n : ab_proj;
    }
    // min_project_n > min_project_p 
    // flip_normal_ = false;
    // if(min_project_n > min_project_p)
    // {
    //     flip_normal_ = true;
    //     min_project_p = min_project_n;
    // }
    min_project_p = min_project_n > min_project_p ? min_project_n : min_project_p;
    min_project_p = std::max(-1.0, min_project_p);

    double angle = acos (min_project_p) * Anlge_PI_Rate ;
    return angle;
}


inline bool LocalVipss::FlipClusterNormal(size_t c_a, size_t c_b) const
{
    const auto& core_ids_a = cluster_core_pt_ids_[c_a];
    const auto& core_ids_b = cluster_core_pt_ids_[c_b];
    std::vector<size_t> valid_ids;
    for(auto id : core_ids_a)
    {
        if(cluster_normal_x_(c_b, id) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    for(auto id : core_ids_b)
    {
        if(cluster_normal_x_(c_a, id) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    size_t p_size = valid_ids.size();
    arma::mat normal_ma(p_size, 3);
    arma::mat normal_mb(p_size, 3);

    for(size_t i = 0; i < valid_ids.size(); ++i)
    {
        size_t id =  valid_ids[i];
        normal_ma(i, 0) =  cluster_normal_x_(c_a, id);
        normal_ma(i, 1) =  cluster_normal_y_(c_a, id);
        normal_ma(i, 2) =  cluster_normal_z_(c_a, id);

        normal_mb(i, 0) =  cluster_normal_x_(c_b, id);
        normal_mb(i, 1) =  cluster_normal_y_(c_b, id);
        normal_mb(i, 2) =  cluster_normal_z_(c_b, id);
    }
    
    return IsFlipNormal(normal_ma, normal_mb);
}

inline double LocalVipss::CalculateClusterPairScore(size_t c_a, size_t c_b, bool& flip) const
{
    const auto& core_ids_a = cluster_core_pt_ids_[c_a];
    const auto& core_ids_b = cluster_core_pt_ids_[c_b];
    std::vector<size_t> valid_ids;
    for(auto id : core_ids_a)
    {
        if(cluster_normal_x_(c_b, id) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    for(auto id : core_ids_b)
    {
        if(cluster_normal_x_(c_a, id) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    size_t p_size = valid_ids.size();
    arma::mat normal_ma(p_size, 3);
    arma::mat normal_mb(p_size, 3);

    for(size_t i = 0; i < valid_ids.size(); ++i)
    {
        size_t id =  valid_ids[i];
        normal_ma(i, 0) =  cluster_normal_x_(c_a, id);
        normal_ma(i, 1) =  cluster_normal_y_(c_a, id);
        normal_ma(i, 2) =  cluster_normal_z_(c_a, id);

        normal_mb(i, 0) =  cluster_normal_x_(c_b, id);
        normal_mb(i, 1) =  cluster_normal_y_(c_b, id);
        normal_mb(i, 2) =  cluster_normal_z_(c_b, id);
    }

    arma::mat dot_mat = normal_ma % normal_mb;
    arma::colvec dot_sum = arma::sum(dot_mat, 1);
    double min_project_p = dot_sum.min();
    double min_project_n = -1.0 * dot_sum.max();

    // min_project_p = min_project_n > min_project_p ? min_project_n : min_project_p;

    // flip_normal_ = false;
    if(min_project_n > min_project_p)
    {
        flip = true;
        min_project_p = min_project_n;
    }
    min_project_p = std::min(1.0, std::max(-1.0, min_project_p));
    double angle = acos (min_project_p) * Anlge_PI_Rate ;

    return angle;
    
    // return CalculateScores(normal_ma, normal_mb);

    // std::vector<arma::vec3> a_normals;
    // std::vector<arma::vec3> b_normals;

    // for(auto id : core_ids_a)
    // {
    //     if(cluster_normal_x_(c_b, id) != 0)
    //     {
    //         arma::vec3 cur_n_b{cluster_normal_x_(c_b, id), cluster_normal_y_(c_b, id), cluster_normal_z_(c_b, id)};
    //         b_normals.push_back(cur_n_b);
    //         arma::vec3 cur_n_a{cluster_normal_x_(c_a, id), cluster_normal_y_(c_a, id), cluster_normal_z_(c_a, id)};
    //         a_normals.push_back(cur_n_a);
    //     }
    // }

    // for(auto id : core_ids_b)
    // {
    //     if(cluster_normal_x_(c_a, id) != 0)
    //     {
    //         arma::vec3 cur_n_b{cluster_normal_x_(c_b, id), cluster_normal_y_(c_b, id), cluster_normal_z_(c_b, id)};
    //         b_normals.push_back(cur_n_b);
    //         arma::vec3 cur_n_a{cluster_normal_x_(c_a, id), cluster_normal_y_(c_a, id), cluster_normal_z_(c_a, id)};
    //         a_normals.push_back(cur_n_a);
    //     }
    // }

    // return CalculateScores(a_normals, b_normals);   

}

void LocalVipss::BuildClusterAdjacentMat()
{

    // printf("00 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    printf("adjacent no zero count : %llu \n", adjacent_mat_.n_nonzero);
    cluster_adjacent_mat_ = adjacent_mat_ * cluster_cores_mat_;
    cluster_adjacent_pt_mat_ = adjacent_mat_ * cluster_cores_mat_ + cluster_cores_mat_;
    cluster_adjacent_flip_mat_.resize(cluster_cores_mat_.n_rows, cluster_cores_mat_.n_rows); 
    // printf("11 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    // for(size_t i = 0; i < pt_num_; ++i)
    // {
    //     cluster_cores_mat_((i, 0), (i, pt_num_ -1));
    //     cluster_adjacent_mat_.row(i) = adjacent_mat_ * cluster_cores_mat_();
    // }
}

void LocalVipss::InitSingleClusterNeiScores(size_t i)
{
    // const arma::sp_imat::iterator start = cluster_adjacent_mat_.row(i).begin();
    // const arma::sp_irowvec::const_iterator end = cluster_adjacent_mat_.row(i).end();
    // const arma::sp_irowvec& nei = cluster_adjacent_mat_.row(i);
    // const arma::sp_irowvec::const_iterator start = nei.begin();
    // const arma::sp_irowvec::const_iterator end = nei.end();
    auto cur_pt = points_[i];
    auto nei_pts = voro_gen_.point_cluster_pts_map_[cur_pt];
    double sum = 0;
    int count = 0;
    arma::mat cur_i_mat(2, 3);
    arma::mat cur_n_mat(2, 3);
    cur_i_mat(0, 0) = cluster_normal_x_(i, i);
    cur_i_mat(0, 1) = cluster_normal_y_(i, i);
    cur_i_mat(0, 2) = cluster_normal_z_(i, i);
    for(auto n_pt : nei_pts)
    {
        if(n_pt == cur_pt) continue;
        size_t n_pos = voro_gen_.point_id_map_[n_pt];
        if(cluster_scores_mat_(i, n_pos) != 0) 
        {
            continue;
        }
        cur_i_mat(1, 0) = cluster_normal_x_(i, n_pos);
        cur_i_mat(1, 1) = cluster_normal_y_(i, n_pos);
        cur_i_mat(1, 2) = cluster_normal_z_(i, n_pos);

        cur_n_mat(0, 0) = cluster_normal_x_(n_pos, i);
        cur_n_mat(0, 1) = cluster_normal_y_(n_pos, i);
        cur_n_mat(0, 2) = cluster_normal_z_(n_pos, i);

        cur_n_mat(1, 0) = cluster_normal_x_(n_pos, n_pos);
        cur_n_mat(1, 1) = cluster_normal_y_(n_pos, n_pos);
        cur_n_mat(1, 2) = cluster_normal_z_(n_pos, n_pos);

        arma::mat dot_res = cur_i_mat % cur_n_mat;
        arma::vec dot_sum = arma::sum(dot_res, 1); 
        double s1 = dot_sum.min();
        double s2 = -dot_sum.max();
        if(s2 > s1) 
        {
            s1 = s2;
            cluster_adjacent_flip_mat_(i, n_pos) = 1;
            cluster_adjacent_flip_mat_(n_pos, i) = 1;
        }
        double score = std::min(1.0, std::max(-1.0, s1));
        score = acos(score) * Anlge_PI_Rate ;
        // printf("cluster id : %d \n", i);
        cluster_scores_mat_(i, n_pos) = score;
        cluster_scores_mat_(n_pos, i) = score;
    }

    // for(auto iter = start; iter != end; ++ iter)
    // {
    //     size_t n_pos = iter.internal_col;
    //     if(i == n_pos) continue;
    //     count ++;
    //     if(cluster_scores_mat_(i, n_pos) != 0) 
    //     {
    //         sum += cluster_scores_mat_(i, n_pos);
    //         continue;
    //     }
    //     cur_i_mat(1, 0) = cluster_normal_x_(i, n_pos);
    //     cur_i_mat(1, 1) = cluster_normal_y_(i, n_pos);
    //     cur_i_mat(1, 2) = cluster_normal_z_(i, n_pos);

    //     cur_n_mat(0, 0) = cluster_normal_x_(n_pos, i);
    //     cur_n_mat(0, 1) = cluster_normal_y_(n_pos, i);
    //     cur_n_mat(0, 2) = cluster_normal_z_(n_pos, i);

    //     cur_n_mat(1, 0) = cluster_normal_x_(n_pos, n_pos);
    //     cur_n_mat(1, 1) = cluster_normal_y_(n_pos, n_pos);
    //     cur_n_mat(1, 2) = cluster_normal_z_(n_pos, n_pos);

    //     // double d_ii_ni = cluster_normal_x_(i, i) * cluster_normal_x_(n_pos, i) + 
    //     //                 cluster_normal_y_(i, i) * cluster_normal_y_(n_pos, i) +
    //     //                 cluster_normal_z_(i, i) * cluster_normal_z_(n_pos, i);
        
    //     // double d_in_nn = cluster_normal_x_(i, n_pos) * cluster_normal_x_(n_pos, n_pos) + 
    //     //                 cluster_normal_y_(i, n_pos) * cluster_normal_y_(n_pos, n_pos) +
    //     //                 cluster_normal_z_(i, n_pos) * cluster_normal_z_(n_pos, n_pos);
    //     arma::mat dot_res = cur_i_mat % cur_n_mat;
    //     arma::vec dot_sum = arma::sum(dot_res, 1); 
    //     double s1 = dot_sum.min();
    //     double s2 = -dot_sum.max();
    //     if(s2 > s1) 
    //     {
    //         s1 = s2;
    //         cluster_adjacent_flip_mat_(i, n_pos) = 1;
    //         cluster_adjacent_flip_mat_(n_pos, i) = 1;
    //     }
    //     double score = std::min(1.0, std::max(-1.0, s1));
    //     score = acos(score) * Anlge_PI_Rate ;
    //     // printf("cluster id : %d \n", i);
    //     cluster_scores_mat_(i, n_pos) = score;
    //     cluster_scores_mat_(n_pos, i) = score;
    //     sum += score;
    // }
    // if(count > 0)
    // {
    //     cluster_scores_vec_[i] = sum / double(count);
    // }
}


void LocalVipss::CalculateSingleClusterNeiScores(size_t i)
{
    const arma::sp_irowvec& nei = cluster_adjacent_mat_.row(i);
    const arma::sp_irowvec::const_iterator start = nei.begin();
    const arma::sp_irowvec::const_iterator end = nei.end();
    // printf("cluster id : %d \n", i);

    int count = 0;
    double sum = 0;
    for(auto iter = start; iter != end; ++ iter)
    {
        size_t n_pos = iter.internal_col;
        if(i == n_pos) continue;
        count ++;
        if(cluster_scores_mat_(i, n_pos) != 0) 
        {
            sum += cluster_scores_mat_(i, n_pos);
            continue;
        }
        bool flip_normal = false;
        double score = CalculateClusterPairScore(i, n_pos, flip_normal);
        if(flip_normal)
        {
            cluster_adjacent_flip_mat_(i, n_pos) = 1;
            cluster_adjacent_flip_mat_(n_pos, i) = 1;
        }
        // printf("cluster id : %d \n", i);
        cluster_scores_mat_(i, n_pos) = score;
        cluster_scores_mat_(n_pos, i) = score;
        update_score_cluster_ids_.insert(n_pos);
        sum += score;
    }
    // double new_cluster_score = sum / (double(count) + 1e-8);
    // cluster_scores_vec_.push_back(new_cluster_score);
    // printf("--------- CalculateClusterPairScore count : %d \n", count);
}   

void LocalVipss::CalculateClusterNeiScores(bool is_init)
{
    size_t c_num = cluster_adjacent_mat_.n_rows;
    // cluster_scores_vec_.resize(c_num);
    for(size_t i = 0; i < c_num; ++i)
    {
        if(is_init)
        {
            InitSingleClusterNeiScores(i);
        } else {
            CalculateSingleClusterNeiScores(i);
        }
    }
}

void LocalVipss::CalculateClusterScores()
{
    size_t c_num = cluster_scores_mat_.n_rows;
    
    // cluster_id_scores_.clear();
    // cluster_id_scores_.resize(c_num);

    // auto t0 = Clock::now();
    // arma::sp_colvec score_sum = arma::sum(cluster_scores_mat_,1);
    // for(auto id :update_score_cluster_ids_)
    // {
    //     if(cluster_scores_vec_.size() <= id)
    //     printf("id : %d \n", id);
    //     // double sum = arma::sum(cluster_scores_mat_.row(id));
    //     cluster_scores_vec_[id] = score_sum(id) / 
    //         (double(cluster_scores_mat_.row(id).n_nonzero) + 1e-10);

    //     // cluster_scores_vec_[id] = arma::mean(arma::nonzeros(cluster_scores_mat_.row(id)));
    // }

    // auto t1 = Clock::now();
    // double average_score_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    // printf("finish average_score_time time : %f ! \n", average_score_time);
    // atuo t0 = Clock::now();

    cluster_scores_vec_ = arma::sum(cluster_scores_mat_,1);
    // arma::sp_colvec non_zero_count_mat(cluster_scores_vec_.n_rows);
    // arma::sp_imat new_pt_mat = cluster_cores_mat_ + cluster_adjacent_pt_mat_;
    // new_pt_mat /= new_pt_mat;
    arma::sp_mat new_pt_mat = cluster_scores_mat_;
    new_pt_mat.for_each([](arma::sp_mat::elem_type& value){value = 1.0;});
    arma::sp_colvec count_mat =  arma::sum(new_pt_mat, 1);
    arma::vec lower_col(c_num, arma::fill::ones);
    lower_col *= 1e-8;
    count_mat += lower_col;
    // count_mat.save("colsum", arma::csv_ascii);
    cluster_scores_vec_ /= count_mat;
    // t1 = Clock::now();
    // average_score_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    // printf("finish average_score_time time 222------- : %f ! \n", average_score_time);

    // for(size_t i = 0; i < c_num; ++i)
    // {
    //     double score = score_sum(i) / double(cluster_scores_mat_.row(i).n_nonzero);
    //     cluster_id_scores_[i] = std::pair<int, double>(i, score);
    // }


    // for(size_t i = 0; i < c_num; ++i)
    // {
    //     double score = score_sum(i) / double(cluster_scores_mat_.row(i).n_nonzero);
    //     cluster_id_scores_.push_back(std::pair<int, double>(i, score));
    // }

    // std::sort(cluster_id_scores_.begin(), cluster_id_scores_.end(), 
    //         [](std::pair<int, double>& a, std::pair<int, double> &b)
    //         {
    //             return a.second > b.second;
    //         });
}

void LocalVipss::GetClusterCenters()
{
    if(cluster_centers_.empty())
    {
        for(size_t i = 0; i < cluster_cores_mat_.n_rows; ++i)
        {
            tetgenmesh::point new_p = new double[3];
            new_p[0] =0;
            new_p[1] =0;
            new_p[2] =0;
            cluster_centers_.push_back(new_p);
        } 
    } else {
        for(size_t i = 0; i < cluster_cores_mat_.n_rows; ++i)
        {
            tetgenmesh::point& cur_p = cluster_centers_[i];
            cur_p[0] =0;
            cur_p[1] =0;
            cur_p[2] =0;
        } 
    }
    
    for(size_t i = 0; i < cluster_cores_mat_.n_rows; ++i)
    {
        tetgenmesh::point& new_p = cluster_centers_[i];
        const arma::sp_irowvec& cores_pt = cluster_cores_mat_.row(i);
        const arma::sp_irowvec::const_iterator start = cores_pt.begin();
        const arma::sp_irowvec::const_iterator end = cores_pt.end();
        int count = 0;
        for(auto iter = start; iter != end; ++iter)
        {
            size_t pt_id = iter.internal_col;
            new_p[0] += points_[pt_id][0];
            new_p[1] += points_[pt_id][1];
            new_p[2] += points_[pt_id][2];
            count ++;
        }
        if(count > 0)
        {
            new_p[0] /= double(count);
            new_p[1] /= double(count);
            new_p[2] /= double(count);
        }
    }
}

void LocalVipss::MergeClusters()
{
    printf(" start to init merge cluster \n");
    std::set<int> visited_clusters;
    std::vector<size_t> merged_cluster_ids;
    // printf("cluster scores num : %zu \n", cluster_id_scores_.size());
    // auto &cluster_scores = cluster_scores_vec_;
    // std::vector<size_t> indices(cluster_scores_vec_.size());
    // std::iota(indices.begin(), indices.end(), 0);
    // std::sort(indices.begin(), indices.end(), [cluster_scores](size_t a , size_t b)
    // {
    //     return cluster_scores[a] > cluster_scores[b];
    // });
    arma::vec c_scores(cluster_scores_vec_);
    // arma::colvec c_scores = arma::conv_to<arma::colvec>::from(cluster_scores_vec_);
    arma::uvec indices = arma::sort_index(c_scores, "descend");
    printf("cluster scores indices num : %zu \n", indices.size());
    // for(auto &ele : cluster_id_scores_)
    for(auto& id : indices)
    {
        if(visited_clusters.find(id) != visited_clusters.end())
        {
            continue;
        }
        visited_clusters.insert(id);
        if(cluster_scores_vec_[id] > angle_threshold_) 
        {
            const arma::sp_rowvec& c_scores_row = cluster_scores_mat_.row(id);
            const arma::sp_rowvec::const_iterator start = c_scores_row.begin();
            const arma::sp_rowvec::const_iterator end = c_scores_row.end();
            double max_score = 0;
            int max_id = -1;
            // auto cur_center = cluster_centers_[ele.first];
            for(auto iter = start; iter != end; ++iter)
            {
                if(visited_clusters.find(iter.internal_col) == visited_clusters.end())
                {
                    // double cur_dist = PtDistance(cluster_centers_[iter.internal_col], cur_center);
                    // cur_dist = cur_dist > 1e-10? cur_dist: 1e-10;
                    double cur_score = *iter;
                    if(cur_score > max_score)
                    {
                        max_id = iter.internal_col;
                        max_score = cur_score;
                    }
                }
            }

            if(max_id != -1)
            {   
                visited_clusters.insert(max_id);
                merged_cluster_ids.push_back(id);
                merged_cluster_ids.push_back(max_id);

                // printf("before add cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);
                arma::sp_irowvec new_cluster_cores = cluster_cores_mat_.row(id) + cluster_cores_mat_.row(max_id);
                AppendRow(cluster_cores_mat_, new_cluster_cores);

                arma::sp_irowvec adjacent_pts_row = cluster_adjacent_pt_mat_.row(id) + cluster_adjacent_pt_mat_.row(max_id);
                AppendRow(cluster_adjacent_pt_mat_, adjacent_pts_row);
                
                // printf("after add cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);

                arma::sp_irowvec adj_row = cluster_adjacent_mat_.row(id) + cluster_adjacent_mat_.row(max_id);
                AppendRow(cluster_adjacent_mat_, adj_row);
                arma::sp_icolvec adj_col = cluster_adjacent_mat_.col(id) + cluster_adjacent_mat_.col(max_id);
                AppendCol(cluster_adjacent_mat_, adj_col);
            }
        }
        
    }
    
    // cluster_cores_mat_.save(out_dir_ + "c_core.mat", arma::arma_ascii);
    merged_cluster_size_ = merged_cluster_ids.size();
    std::sort(merged_cluster_ids.begin(), merged_cluster_ids.end(), std::greater<>());
    // arma::uvec delete_ids = merged_cluster_ids;
    // cluster_cores_mat_.shed_rows(delete_ids);
    for(size_t id : merged_cluster_ids)
    {
        cluster_cores_mat_.shed_row(id);
        cluster_adjacent_pt_mat_.shed_row(id);

        cluster_scores_mat_.shed_row(id);
        cluster_scores_mat_.shed_col(id);

        cluster_adjacent_mat_.shed_row(id);
        cluster_adjacent_mat_.shed_col(id);

        cluster_normal_x_.shed_row(id);
        cluster_normal_y_.shed_row(id);
        cluster_normal_z_.shed_row(id);

        cluster_adjacent_flip_mat_.shed_row(id);
        cluster_adjacent_flip_mat_.shed_col(id);

        cluster_core_pt_ids_.erase(std::next(cluster_core_pt_ids_.begin(), id));
        // cluster_scores_vec_.erase(std::next(cluster_scores_vec_.begin(), id));
    }
}

void LocalVipss::AppendRow(arma::sp_imat& in_mat,  arma::sp_irowvec& append_row)
{
    size_t rows = in_mat.n_rows;
    size_t cols = in_mat.n_cols;

    arma::sp_imat new_mat(rows+1, cols);
    new_mat.rows(0, rows-1) = in_mat.rows(0, rows-1);
    new_mat.row(rows) = append_row;
    in_mat = new_mat;

}

void LocalVipss::AppendCol(arma::sp_imat& in_mat,  arma::sp_icolvec& append_col)
{
    size_t rows = in_mat.n_rows;
    size_t cols = in_mat.n_cols;

    arma::sp_imat new_mat(rows, cols + 1);
    new_mat.cols(0, cols-1) = in_mat.cols(0, cols-1);
    new_mat.col(cols) = append_col;
    in_mat = new_mat;
}
        
void LocalVipss::UpdateClusterScoreMat()
{
    size_t c_num = cluster_cores_mat_.n_rows;
    size_t s_rows = cluster_scores_mat_.n_rows;
    size_t s_cols = cluster_scores_mat_.n_cols;

    // printf("cluster_scores_mat_ size %d updated to %d \n", s_rows, c_num);

    // printf("00cluster_scores_mat_ rows : %d , cols : %d \n", cluster_scores_mat_.n_rows, cluster_scores_mat_.n_cols);
    // auto t0 = Clock::now();
    arma::sp_mat new_scores_mat(c_num, c_num);
    new_scores_mat(0, 0, arma::size(s_rows -1, s_cols -1)) = cluster_scores_mat_(0, 0, arma::size(s_rows -1, s_cols -1));
    cluster_scores_mat_ = new_scores_mat;

    arma::sp_imat new_cluster_adjacent_flip_mat(c_num, c_num);
    new_cluster_adjacent_flip_mat(0, 0, arma::size(s_rows -1, s_cols -1)) =
             cluster_adjacent_flip_mat_(0, 0, arma::size(s_rows -1, s_cols -1));
    cluster_adjacent_flip_mat_ = new_cluster_adjacent_flip_mat;

    // auto t1 = Clock::now();
    // double update_score_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    // printf("finish update_score_time time : %f ! \n", update_score_time);

    // printf("cluster_scores_mat_ rows : %d , cols : %d \n", cluster_scores_mat_.n_rows, cluster_scores_mat_.n_cols);
    
    // printf("cluster_adjacent_mat_ rows : %d , cols : %d \n", cluster_adjacent_mat_.n_rows, cluster_adjacent_mat_.n_cols);
    
    for(size_t i = s_rows; i < c_num; ++i)
    {
        CalculateSingleClusterNeiScores(i);
    }
    printf("cluster_scores_mat_ size %zu updated to %zu \n", s_rows, c_num);
}


void LocalVipss::UpdateClusterNormals()
{
    size_t c_num = cluster_cores_mat_.n_rows;
    size_t s_rows = cluster_normal_x_.n_rows;
    size_t s_cols = cluster_normal_x_.n_cols;

    // printf("cluster num : %zu , s_rows : %zu, s_cols: %zu \n", c_num, s_rows, s_cols);

    arma::sp_mat new_normal_x(c_num, s_cols);
    arma::sp_mat new_normal_y(c_num, s_cols);
    arma::sp_mat new_normal_z(c_num, s_cols);

    if(s_rows > 0)
    {
        new_normal_x.rows(0, s_rows-1) = cluster_normal_x_.rows(0, s_rows-1);
        new_normal_y.rows(0, s_rows-1) = cluster_normal_y_.rows(0, s_rows-1);
        new_normal_z.rows(0, s_rows-1) = cluster_normal_z_.rows(0, s_rows-1);
    }
    cluster_normal_x_ = new_normal_x;
    cluster_normal_y_ = new_normal_y;
    cluster_normal_z_ = new_normal_z;

    vipss_time_stats_.clear();
    for(size_t i = s_rows; i < c_num; ++i)
    {
        CalculateClusterNormals(i);
    }
    cluster_ptn_vipss_time_stats_.push_back(vipss_time_stats_);
}

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

void LocalVipss::FlipClusterNormalsByScores()
{
    size_t c_num = cluster_cores_mat_.n_rows;

    std::queue<size_t> cluster_queued_ids;
    cluster_queued_ids.push(0);
    std::set<size_t> visited_cluster_ids;
    std::set<size_t> flipped_cluster_ids;
    flipped_cluster_ids.insert(0);
    while(!cluster_queued_ids.empty())
    {
        size_t cur_cid = cluster_queued_ids.front();
        cluster_queued_ids.pop();
        if(visited_cluster_ids.find(cur_cid) != visited_cluster_ids.end()) continue;

        const arma::sp_irowvec& adj_row = cluster_adjacent_mat_.row(cur_cid);
        const arma::sp_irowvec::const_iterator start = adj_row.begin();
        const arma::sp_irowvec::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {   
            size_t n_cid = iter.internal_col;
            if(flipped_cluster_ids.find(n_cid) != flipped_cluster_ids.end()) continue;
            flipped_cluster_ids.insert(n_cid);
            // CalculateClusterPairScore(cur_cid, n_cid);
            // if(cluster_adjacent_flip_mat_(cur_cid, n_cid) == 1)
            if(FlipClusterNormal(cur_cid, n_cid))
            {
                cluster_normal_x_.row(n_cid) *= -1.0;
                cluster_normal_y_.row(n_cid) *= -1.0;
                cluster_normal_z_.row(n_cid) *= -1.0;
            }
            cluster_queued_ids.push(n_cid);
        }
    }
}

void LocalVipss::BuildClusterMST()
{
    std::set<std::string> visited_edge_ids;
    std::vector<C_Edege> tree_edges;

    auto cmp = [](const C_Edege& left, const C_Edege& right) { return left.score_ > right.score_; };
    std::priority_queue<C_Edege, std::vector<C_Edege>, decltype(cmp)> edge_priority_queue(cmp);

    C_Edege st_e(0,0);
    edge_priority_queue.push(st_e);
    std::set<size_t> visited_vids;
    while(!edge_priority_queue.empty())
    {
        C_Edege cur_e = edge_priority_queue.top();
        edge_priority_queue.pop();

        if(visited_vids.find(cur_e.c_b_ ) != visited_vids.end()) continue;
        visited_vids.insert(cur_e.c_b_);
        if(cur_e.c_a_ != cur_e.c_b_)
        {
            tree_edges.push_back(cur_e);
        }

        size_t cur_pid = cur_e.c_b_ ;
        const arma::sp_irowvec& adj_row = cluster_adjacent_mat_.row(cur_pid);
        // printf("row %d contains no zero num : %d \n", cur_pid, adj_row.n_nonzero);
        const arma::sp_irowvec::const_iterator start = adj_row.begin();
        const arma::sp_irowvec::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {
            size_t n_id = iter.internal_col;
            if(visited_vids.find(n_id) != visited_vids.end()) continue;
            if(n_id == cur_pid) continue;
            C_Edege edge(cur_pid, n_id);
            edge.score_ = cluster_scores_mat_(cur_pid, n_id);
            edge_priority_queue.push(edge);
        }
    }
    size_t c_num = cluster_cores_mat_.n_rows;
    cluster_MST_mat_.resize(c_num, c_num);
    for(const auto& edge: tree_edges)
    {
        size_t i = edge.c_a_;
        size_t j = edge.c_b_;
        cluster_MST_mat_(i,j) = 1;
        cluster_MST_mat_(j,i) = 1;
    }
}

void LocalVipss::FlipClusterNormalsByMST()
{
    size_t c_num = cluster_cores_mat_.n_rows;
    std::queue<size_t> cluster_queued_ids;
    cluster_queued_ids.push(0);
    std::set<size_t> visited_cluster_ids;
    std::set<size_t> flipped_cluster_ids;
    flipped_cluster_ids.insert(0);
    while(!cluster_queued_ids.empty())
    {
        size_t cur_cid = cluster_queued_ids.front();
        cluster_queued_ids.pop();
        if(visited_cluster_ids.find(cur_cid) != visited_cluster_ids.end()) continue;

        const arma::sp_irowvec& adj_row = cluster_MST_mat_.row(cur_cid);
        const arma::sp_irowvec::const_iterator start = adj_row.begin();
        const arma::sp_irowvec::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {   
            size_t n_cid = iter.internal_col;
            if(flipped_cluster_ids.find(n_cid) != flipped_cluster_ids.end()) continue;
            flipped_cluster_ids.insert(n_cid);
            // CalculateClusterPairScore(cur_cid, n_cid);
            // FlipClusterNormal
            // bool flip = false;
            // CalculateClusterPairScore(cur_cid, n_cid, flip);
            // bool new_flip = FlipClusterNormal(cur_cid, n_cid);
            // bool pre_flip = cluster_adjacent_flip_mat_(cur_cid, n_cid);
            // printf("flip %d , new flip %d, pre_flip \n", flip, new_flip);
            if(FlipClusterNormal(cur_cid, n_cid))
            // if(cluster_adjacent_flip_mat_(cur_cid, n_cid) > 0)
            {
                cluster_normal_x_.row(n_cid) *= -1.0;
                cluster_normal_y_.row(n_cid) *= -1.0;
                cluster_normal_z_.row(n_cid) *= -1.0;
            }
            cluster_queued_ids.push(n_cid);
        }
    }
}

void LocalVipss::FlipClusterNormalsByMinST()
{
    size_t pt_num = adjacent_mat_.n_rows;
    std::set<std::string> visited_edge_ids;
    std::vector<C_Edege> tree_edges;

    auto cmp = [](C_Edege& left, C_Edege& right) { return left.score_ > right.score_; };
    std::priority_queue<C_Edege, std::vector<C_Edege>, decltype(cmp)> edge_priority_queue(cmp);
 
    C_Edege st_e(0,0);
    edge_priority_queue.push(st_e);
    std::set<size_t> visited_vids;
    while(!edge_priority_queue.empty())
    {
        C_Edege cur_e = edge_priority_queue.top();
        edge_priority_queue.pop();

        if(visited_vids.find(cur_e.c_b_ ) != visited_vids.end()) continue;
        visited_vids.insert(cur_e.c_b_);
        if(cur_e.c_a_ != cur_e.c_b_)
        {
            tree_edges.push_back(cur_e);
        }

        size_t cur_pid = cur_e.c_b_ ;
        const arma::sp_irowvec& adj_row = adjacent_mat_.row(cur_pid);
        // printf("row %d contains no zero num : %d \n", cur_pid, adj_row.n_nonzero);
        const arma::sp_irowvec::const_iterator start = adj_row.begin();
        const arma::sp_irowvec::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {
            size_t n_id = iter.internal_col;
            if(visited_vids.find(n_id) != visited_vids.end()) continue;
            if(n_id == cur_pid) continue;
            C_Edege edge(cur_pid, n_id);
            edge.score_ = PtDistance(points_[cur_pid], points_[n_id]);
            edge_priority_queue.push(edge);
        }
    }

    std::string out_tree_path = out_dir_ + filename_ + "tree.obj";
    std::vector<double> vts ;
    for(auto& p : points_)
    {
        vts.push_back(p[0]);
        vts.push_back(p[1]);
        vts.push_back(p[2]);
    }
    std::vector<unsigned int> edges_out;
    for(auto& e : tree_edges) 
    {
        edges_out.push_back(e.c_a_);
        edges_out.push_back(e.c_b_);
    }

    writeObjFile_line(out_tree_path,vts, edges_out);
}

void LocalVipss::OuputPtN(const std::string& out_path, bool orient_normal)
{
    out_pts_.clear();
    out_normals_.clear();
    std::vector<double>& pts = out_pts_;
    std::vector<double>& normals = out_normals_;
    
    // cluster_cores_mat_.save(out_dir_ + filename_ + "cores_mat_out.mat", arma::raw_ascii );
    // printf("cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);
    for(size_t i = 0; i < cluster_cores_mat_.n_rows; ++i)
    {
        // printf("cluster id : %d \n", i);
        const arma::sp_irowvec& cur_row = cluster_cores_mat_.row(i);
        const arma::sp_irowvec::const_iterator start = cur_row.begin();
        const arma::sp_irowvec::const_iterator end = cur_row.end();
        for(auto iter = start; iter != end; ++iter)
        {
            size_t p_id = iter.internal_col;
            // printf("p_id : %d \n", p_id);
            auto pt = points_[p_id];
            pts.push_back(pt[0]);
            pts.push_back(pt[1]);
            pts.push_back(pt[2]);
            normals.push_back(cluster_normal_x_(i, p_id));
            normals.push_back(cluster_normal_y_(i, p_id));
            normals.push_back(cluster_normal_z_(i, p_id));
        }
    }
    printf("out pt size : %zu \n", pts.size() / 3);

    // if(orient_normal)
    // {
    //     ORIENT::OrientPointNormals(pts, normals);
    // }
    
    writePLYFile_VN(out_path, pts, normals);
}

void LocalVipss::WriteVipssTimeLog()
{
    std::ofstream myfile;
    std::string csv_path = out_dir_ + "_vipss_time.csv";
    myfile.open(csv_path);
    
    double time_sum = 0.0;
    size_t count = 0;
    for(auto& v_times: cluster_ptn_vipss_time_stats_)
    {
        for(auto& ele : v_times)
        {
            time_sum += ele.second;
            myfile << std::to_string(count) << ",";
            myfile << std::to_string(ele.first) << ",";
            myfile << std::to_string(ele.second) << std::endl;
        }
        count ++;
    }
    printf("Vipss total time : %f \n", time_sum);
}

// void LocalVipss::SaveCluster()
// {
//     for(size_t i = 0; i < cluster_cores_mat_.n_rows; ++i)
//     {
//         const arma::sp_irowvec& cur_row = cluster_cores_mat_.row(i);
//         if(cur_row.n_nonzero < 3) continue;
        
//         const arma::sp_irowvec::const_iterator start = cur_row.begin();
//         const arma::sp_irowvec::const_iterator end = cur_row.end();
//         std::set<P3tr> key_pts;
//         for(auto iter = start; iter != end; ++iter)
//         {
//             key_pts.insert(points_[iter.internal_col]);
//         } 
//         std::vector<P3tr> nei_pts;
//         const arma::sp_irowvec& cur_pt_row = cluster_adjacent_pt_mat_.row(i);
//         const arma::sp_irowvec::const_iterator pt_start = cur_pt_row.begin();
//         const arma::sp_irowvec::const_iterator pt_end = cur_pt_row.end();
//         for(auto iter = pt_start; iter != pt_end; ++iter)
//         {
//             auto cur_pt = points_[iter.internal_col];
//             if(key_pts.find(cur_pt) == key_pts.end())
//             {
//                 nei_pts.push_back(cur_pt);
//             }
//         }
//         std::vector<P3tr> key_pts_vec(key_pts.begin(), key_pts.end());
//         std::string out_path = out_dir_ + filename_ + "_cluster_" + std::to_string(i);
//         // SaveClusterPts(out_path, key_pts_vec, nei_pts);
//     }
// }

// void LocalVipss::SaveClusterPts(const std::string& path,
//                             const std::vector<P3tr>& key_pts, 
//                             const std::vector<P3tr>& nei_pts)
// {
//     std::vector<double> pts;
//     std::vector<uint8_t> colors;
//     size_t key_num = key_pts.size();
//     size_t nei_num = nei_pts.size();
//     pts.resize(3*(key_num + nei_num));
//     colors.resize(3*(key_num + nei_num));
//     for(size_t i = 0; i < key_num + nei_num; ++i)
//     {
//         if(i < key_num)
//         {
//             pts[3*i] = key_pts[i][0];
//             pts[3*i + 1] = key_pts[i][1];
//             pts[3*i + 2] = key_pts[i][2];
//             colors[3*i] = 255;
//             colors[3*i + 1] = 0;
//             colors[3*i + 2] = 0;
//         } else {
//             size_t id = i - key_num;
//             pts[3*id] = key_pts[id][0];
//             pts[3*id + 1] = key_pts[id][1];
//             pts[3*id + 2] = key_pts[id][2];
//             colors[3*id] = 0;
//             colors[3*id + 1] = 0;
//             colors[3*id + 2] = 255;
//         } 
//     }
//     writePLYFile_CO(path, pts, colors);
// }



void LocalVipss::Run()
{
    // printf("00000 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);

    double total_time = 0;
    auto t0 = Clock::now();
    BuildClusterAdjacentMat();
    // printf("finish build cluster adjacent mat ! \n");

    // printf("2 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    BuidClusterCoresPtIds();
    auto t1 = Clock::now();
    double build_mat_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("finish init core pt ids time : %f ! \n", build_mat_time);

    InitNormalWithVipss();
    auto t12 = Clock::now();
    double normal_estimate_time = std::chrono::nanoseconds(t12 - t1).count()/1e9;
    printf("finish init cluster normals time : %f ! \n", normal_estimate_time);

    total_time += build_mat_time + normal_estimate_time;
    

    // std::string init_ptn_path = out_dir_ + filename_ + "_init_ptn";
    // OuputPtN(init_ptn_path);

    auto t2 = Clock::now();
    CalculateClusterNeiScores(true);
    auto t3 = Clock::now();
    double scores_time = std::chrono::nanoseconds(t3 - t2).count()/1e9;
    printf("finish init cluster neigh scores time : %f ! \n", scores_time);

    CalculateClusterScores();

    auto t34 = Clock::now();
    double cluster_scores_time = std::chrono::nanoseconds(t34 - t3).count()/1e9;
    printf("finish calculate cluster scores time : %f ! \n", cluster_scores_time);

    total_time += scores_time;
    size_t iter_num = 1;
    while(merged_cluster_size_ > 0)
    {

        auto tt0 = Clock::now();
        // GetClusterCenters();
        MergeClusters();
        auto tt1 = Clock::now();
        double merge_time = std::chrono::nanoseconds(tt1 - tt0).count()/1e9;
        total_time += merge_time;
        printf("MergeClusters time : %f ! \n", merge_time);
        printf("iter %zu merged cluster number : %zu \n", iter_num, merged_cluster_size_);
        if(merged_cluster_size_ == 0)
        {
            break;
        }

        auto tt2 = Clock::now();
        UpdateClusterCoresPtIds();
        UpdateClusterNormals();
        auto tt3 = Clock::now();
        double update_normal_time = std::chrono::nanoseconds(tt3 - tt2).count()/1e9;
        printf("update_normal time : %f ! \n", update_normal_time);
        
        total_time += update_normal_time;
        // printf("cluster_adjacent_mat_ number : %d \n", cluster_adjacent_mat_.n_rows);
        auto tt4 = Clock::now();
        update_score_cluster_ids_.clear();
        UpdateClusterScoreMat();
        
        auto tt5 = Clock::now();
        double update_score_time = std::chrono::nanoseconds(tt5 - tt4).count()/1e9;
        printf("update_score time : %f ! \n", update_score_time);
        total_time += update_score_time;

        CalculateClusterScores();
        auto tt6 = Clock::now();
        double calculate_score_time = std::chrono::nanoseconds(tt6 - tt5).count()/1e9;
        printf("calculate_score_time : %f ! \n", calculate_score_time);

        total_time += calculate_score_time;
        // std::string init_ptn_path_iter = out_dir_ + filename_ + "_ptn" + std::to_string(iter_num);
        // OuputPtN(init_ptn_path_iter);
        // std::string init_ptn_path_iter2 = out_dir_ + filename_ + "_orient_ptn" + std::to_string(iter_num);
        // OuputPtN(init_ptn_path_iter2, true);
        iter_num ++;
        if(iter_num >= max_iter_) break;
    }

    auto ti00 = Clock::now();
    BuildClusterMST();
    auto ti11 = Clock::now();
    double MST_time = std::chrono::nanoseconds(ti11 - ti00).count()/1e9;
    total_time += MST_time;
    printf("normal MST_time time used : %f \n", MST_time);
 
    auto ti0 = Clock::now();
    FlipClusterNormalsByMST();
    // FlipClusterNormalsByScores();
    auto ti1 = Clock::now();
    double flip_time = std::chrono::nanoseconds(ti1 - ti0).count()/1e9;
    total_time += flip_time;
    printf("normal flip time used : %f \n", flip_time);
    // flipClusterNormalsByScores();
    std::string init_ptn_path_iter = out_dir_ + filename_ + "_flipped" + std::to_string(iter_num);
    OuputPtN(init_ptn_path_iter);
    WriteVipssTimeLog();
    printf("total time used : %f \n", total_time);

    if(use_hrbf_surface_)
    {
        vipss_api_.is_surfacing_ = true;
        vipss_api_.run_vipss(out_pts_, out_normals_);
    }
}
