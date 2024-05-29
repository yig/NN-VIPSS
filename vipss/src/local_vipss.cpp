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

    printf("read point size : %d \n", in_pts.size());
    // printf("read point size : %d \n", in_pts.size());
    // voro_gen_.loadData(path);
    voro_gen_.loadData(in_pts);
    voro_gen_.InitMesh();
    voro_gen_.BuildPtIdMap();
    voro_gen_.BuildAdjecentMat();

    adjacent_mat_ = voro_gen_.pt_adjecent_mat_;
    points_ = voro_gen_.points_;

    printf("adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    pt_num_ = points_.size();

    cluster_cores_mat_.resize(pt_num_, pt_num_);
    cluster_cores_mat_.eye();
    
    printf("input point size : %d \n", pt_num_);

    vipss_api_.Set_RBF_PARA();
    vipss_api_.is_surfacing_ = false;

    cluster_normal_x_.resize(pt_num_, pt_num_);
    cluster_normal_y_.resize(pt_num_, pt_num_);
    cluster_normal_z_.resize(pt_num_, pt_num_);

    cluster_scores_mat_.resize(pt_num_, pt_num_);
    cluster_adjacent_mat_.resize(pt_num_, pt_num_);
    printf("finish local vipss initilization ! \n");
}

std::vector<size_t> LocalVipss::GetClusterCoreIds(size_t cluster_id)
{
    arma::sp_irowvec cluster_core_ids = cluster_cores_mat_.row(cluster_id);
    arma::sp_irowvec::const_iterator start = cluster_core_ids.begin();
    arma::sp_irowvec::const_iterator end = cluster_core_ids.end();

    std::vector<size_t> core_ids;
    for ( arma::sp_irowvec::const_iterator i = start; i != end; ++i )
    {
        core_ids.push_back(i.internal_col);
    }
    return core_ids;
}

std::vector<size_t> LocalVipss::GetClusterPtIds(size_t cluster_id)
{
    // arma::sp_ivec cluster_core_ids = cluster_cores_mat_.row(cluster_id);
    

    // arma::sp_ivec cluster_pt_ids = adjacent_mat_ * cluster_core_ids;
    size_t nozero_count = cluster_adjacent_pt_mat_.n_nonzero;

    arma::sp_irowvec cluster_pt_ids = cluster_adjacent_pt_mat_.row(cluster_id);

    // printf("cluster_pt_ids no zero count : %d \n", cluster_pt_ids.n_nonzero);

    arma::sp_irowvec::const_iterator start = cluster_pt_ids.begin();
    arma::sp_irowvec::const_iterator end = cluster_pt_ids.end();
    std::vector<size_t> pt_ids;
    for ( arma::sp_irowvec::const_iterator i = start; i != end; ++i )
    {
        // printf("i pos : %d \n", i.internal_pos);
        pt_ids.push_back(i.internal_col);
    }
    return pt_ids;
}

std::vector<double> LocalVipss::GetClusterVerticesFromIds(const std::vector<size_t>& pt_ids) 
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

std::vector<double> LocalVipss::GetClusterVertices(size_t cluster_id)
{
    auto cluster_pt_ids = GetClusterPtIds(cluster_id);
    return GetClusterVerticesFromIds(cluster_pt_ids);
}


size_t LocalVipss::GetClusterIdFromCorePtId(const size_t pid)
{
    arma::sp_icolvec cluster_col = cluster_cores_mat_.col(pid);
    arma::sp_icolvec::const_iterator start = cluster_col.begin();
    return (size_t)start.internal_pos;
}

void LocalVipss::CalculateClusterNormals(size_t cluster_id)
{
    auto cluster_pt_ids = GetClusterPtIds(cluster_id);
    // printf("cluster_pt_ids size : %d \n", cluster_pt_ids.size());
    auto cluster_vts = GetClusterVerticesFromIds(cluster_pt_ids);
    // printf("cluster pt size : %d \n", cluster_vts.size()/3);
    vipss_api_.run_vipss(cluster_vts);
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
    printf("cluster num : %d \n", cluster_num);
    
    for(size_t i =0; i < cluster_num; ++i)
    {
        CalculateClusterNormals(i);
    }
}

double M_PI2  = 2*acos(0.0);

double LocalVipss::CalculateScores(std::vector<arma::vec3>& a_normals, std::vector<arma::vec3>& b_normals)
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
    flip_normal_ = false;
    if(min_project_n > min_project_p)
    {
        flip_normal_ = true;
        min_project_p = min_project_n;
    }
    // min_project_p = min_project_n > min_project_p ? min_project_n : min_project_p;
    // min_project_p = std::min(1.0, min_project_p);

    double angle = acos (min_project_p) * 180.0 / M_PI2 ;
    return angle;
}

double LocalVipss::CalculateClusterPairScore(size_t c_a, size_t c_b)
{
    auto core_ids_a = GetClusterCoreIds(c_a);
    auto core_ids_b = GetClusterCoreIds(c_b);

    std::vector<arma::vec3> a_normals;
    std::vector<arma::vec3> b_normals;

    for(auto id : core_ids_a)
    {
        if(cluster_normal_x_(c_b, id) != 0)
        {
            arma::vec3 cur_n_b{cluster_normal_x_(c_b, id), cluster_normal_y_(c_b, id), cluster_normal_z_(c_b, id)};
            b_normals.push_back(cur_n_b);
            arma::vec3 cur_n_a{cluster_normal_x_(c_a, id), cluster_normal_y_(c_a, id), cluster_normal_z_(c_a, id)};
            a_normals.push_back(cur_n_a);
        }
    }

    for(auto id : core_ids_b)
    {
        if(cluster_normal_x_(c_a, id) != 0)
        {
            arma::vec3 cur_n_b{cluster_normal_x_(c_b, id), cluster_normal_y_(c_b, id), cluster_normal_z_(c_b, id)};
            b_normals.push_back(cur_n_b);
            arma::vec3 cur_n_a{cluster_normal_x_(c_a, id), cluster_normal_y_(c_a, id), cluster_normal_z_(c_a, id)};
            a_normals.push_back(cur_n_a);
        }
    }

    return CalculateScores(a_normals, b_normals);   

}

void LocalVipss::BuildClusterAdjacentMat()
{

    // printf("00 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    printf("adjacent no zero count : %d \n", adjacent_mat_.n_nonzero);
    cluster_adjacent_mat_ = adjacent_mat_ * cluster_cores_mat_;
    cluster_adjacent_pt_mat_ = adjacent_mat_ * cluster_cores_mat_ + cluster_cores_mat_;
    // printf("11 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    // for(size_t i = 0; i < pt_num_; ++i)
    // {
    //     cluster_cores_mat_((i, 0), (i, pt_num_ -1));
    //     cluster_adjacent_mat_.row(i) = adjacent_mat_ * cluster_cores_mat_();
    // }
}

void LocalVipss::CalculateSingleClusterNeiScores(size_t i)
{
    arma::sp_irowvec nei = cluster_adjacent_mat_.row(i);
    const arma::sp_irowvec::const_iterator start = nei.begin();
    const arma::sp_irowvec::const_iterator end = nei.end();
    // printf("cluster id : %d \n", i);
    for(auto iter = start; iter != end; ++ iter)
    {
        size_t n_pos = iter.internal_col;
        if(i == n_pos) continue;
        if(cluster_scores_mat_(i, n_pos) != 0)  continue;
        double score = CalculateClusterPairScore(i, n_pos);
        // printf("cluster id : %d \n", i);
        cluster_scores_mat_(i, n_pos) = score;
        cluster_scores_mat_(n_pos, i) = score;
    }
}   

void LocalVipss::CalculateClusterNeiScores()
{
    size_t c_num = cluster_adjacent_mat_.n_rows;

    for(size_t i = 0; i < c_num; ++i)
    {
        CalculateSingleClusterNeiScores(i);
    }
}

void LocalVipss::CalculateClusterScores()
{
    
    size_t c_num = cluster_scores_mat_.n_rows;
    

    cluster_id_scores_.clear();
    for(size_t i = 0; i < c_num; ++i)
    {
        arma::sp_rowvec c_mat = cluster_scores_mat_.row(i);
        double score = arma::sum(c_mat) / double(c_mat.n_nonzero);
        cluster_id_scores_.push_back(std::pair<int, double>(i, score));
    }
    std::sort(cluster_id_scores_.begin(), cluster_id_scores_.end(), 
            [](std::pair<int, double>& a, std::pair<int, double> &b)
            {
                return a.second > b.second;
            });

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
        arma::sp_irowvec cores_pt = cluster_cores_mat_.row(i);
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
    std::set<int> visited_clusters;
    std::vector<size_t> merged_cluster_ids;
    printf("cluster scores num : %d \n", cluster_id_scores_.size());
    for(auto &ele : cluster_id_scores_)
    {
        if(visited_clusters.find(ele.first) != visited_clusters.end())
        {
            continue;
        }

        visited_clusters.insert(ele.first);

        if(ele.second > angle_threshold_)
        {
            arma::sp_rowvec c_scores_row =  cluster_scores_mat_.row(ele.first);
            const arma::sp_rowvec::const_iterator start = c_scores_row.begin();
            const arma::sp_rowvec::const_iterator end = c_scores_row.end();
            double max_score = 0;
            int max_id = -1;
            auto cur_center = cluster_centers_[ele.first];
            for(auto iter = start; iter != end; iter++)
            {
                if(visited_clusters.find(iter.internal_col) == visited_clusters.end())
                {
                    double cur_dist = PtDistance(cluster_centers_[iter.internal_col], cur_center);
                    cur_dist = cur_dist > 1e-10? cur_dist: 1e-10;
                    double cur_score = *iter/cur_dist;
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
                merged_cluster_ids.push_back(ele.first);
                merged_cluster_ids.push_back(max_id);

                // printf("before add cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);
                arma::sp_irowvec new_cluster_cores = cluster_cores_mat_.row(ele.first) + cluster_cores_mat_.row(max_id);
                AppendRow(cluster_cores_mat_, new_cluster_cores);

                arma::sp_irowvec adjacent_pts_row = cluster_adjacent_pt_mat_.row(ele.first) + cluster_adjacent_pt_mat_.row(max_id);
                AppendRow(cluster_adjacent_pt_mat_, adjacent_pts_row);
                
                // printf("after add cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);

                arma::sp_irowvec adj_row = cluster_adjacent_mat_.row(ele.first) + cluster_adjacent_mat_.row(max_id);
                AppendRow(cluster_adjacent_mat_, adj_row);
                arma::sp_icolvec adj_col = cluster_adjacent_mat_.col(ele.first) + cluster_adjacent_mat_.col(max_id);
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

    arma::sp_mat new_scores_mat(c_num, c_num);
    new_scores_mat(0, 0, arma::size(s_rows -1, s_cols -1)) = cluster_scores_mat_(0, 0, arma::size(s_rows -1, s_cols -1));
    cluster_scores_mat_ = new_scores_mat;

    // printf("cluster_scores_mat_ rows : %d , cols : %d \n", cluster_scores_mat_.n_rows, cluster_scores_mat_.n_cols);
    
    // printf("cluster_adjacent_mat_ rows : %d , cols : %d \n", cluster_adjacent_mat_.n_rows, cluster_adjacent_mat_.n_cols);
    for(size_t i = s_rows; i < c_num; ++i)
    {
        CalculateSingleClusterNeiScores(i);
    }
    printf("cluster_scores_mat_ size %d updated to %d \n", s_rows, c_num);
}


void LocalVipss::UpdateClusterNormals()
{
    size_t c_num = cluster_cores_mat_.n_rows;
    size_t s_rows = cluster_normal_x_.n_rows;
    size_t s_cols = cluster_normal_x_.n_cols;

    arma::sp_mat new_normal_x(c_num, s_cols);
    arma::sp_mat new_normal_y(c_num, s_cols);
    arma::sp_mat new_normal_z(c_num, s_cols);

    new_normal_x.rows(0, s_rows-1) = cluster_normal_x_.rows(0, s_rows-1);
    new_normal_y.rows(0, s_rows-1) = cluster_normal_y_.rows(0, s_rows-1);
    new_normal_z.rows(0, s_rows-1) = cluster_normal_z_.rows(0, s_rows-1);

    cluster_normal_x_ = new_normal_x;
    cluster_normal_y_ = new_normal_y;
    cluster_normal_z_ = new_normal_z;

    for(size_t i = s_rows; i < c_num; ++i)
    {
        CalculateClusterNormals(i);
    }
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

void LocalVipss::flipClusterNormalsByScores()
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

        arma::sp_irowvec adj_row = cluster_adjacent_mat_.row(cur_cid);
        const arma::sp_irowvec::const_iterator start = adj_row.begin();
        const arma::sp_irowvec::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {   
            size_t n_cid = iter.internal_col;
            if(flipped_cluster_ids.find(n_cid) != flipped_cluster_ids.end()) continue;
            flipped_cluster_ids.insert(n_cid);
            CalculateClusterPairScore(cur_cid, n_cid);
            if(flip_normal_)
            {
                cluster_normal_x_.row(n_cid) *= -1.0;
                cluster_normal_y_.row(n_cid) *= -1.0;
                cluster_normal_z_.row(n_cid) *= -1.0;
            }
            cluster_queued_ids.push(n_cid);
        }
    }
}

void LocalVipss::flipClusterNormalsByMinST()
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
        arma::sp_irowvec adj_row = adjacent_mat_.row(cur_pid);
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
    std::vector<double> pts;
    std::vector<double> normals;
    
    // cluster_cores_mat_.save(out_dir_ + filename_ + "cores_mat_out.mat", arma::raw_ascii );
    // printf("cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);
    for(size_t i = 0; i < cluster_cores_mat_.n_rows; ++i)
    {
        // printf("cluster id : %d \n", i);
        arma::sp_irowvec cur_row = cluster_cores_mat_.row(i);
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
    printf("out pt size : %d \n", pts.size() / 3);

    if(orient_normal)
    {
        ORIENT::OrientPointNormals(pts, normals);
    }
    
    writePLYFile_VN(out_path, pts, normals);
}


void LocalVipss::Run()
{
    // printf("00000 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);

    double total_time = 0;
    auto t0 = Clock::now();
    BuildClusterAdjacentMat();
    // printf("finish build cluster adjacent mat ! \n");

    // printf("2 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    
    InitNormalWithVipss();
    auto t1 = Clock::now();
    double estimate_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("finish init cluster normals time : %f ! \n", estimate_time);

    total_time += estimate_time;
    

    std::string init_ptn_path = out_dir_ + filename_ + "_init_ptn";
    OuputPtN(init_ptn_path);

    auto t2 = Clock::now();
    CalculateClusterNeiScores();
    CalculateClusterScores();
    auto t3 = Clock::now();
    double scores_time = std::chrono::nanoseconds(t3 - t2).count()/1e9;
    printf("finish init cluster scores time : %f ! \n", scores_time);

    total_time += scores_time;

    // printf("0 cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);

    size_t iter_num = 1;
    while(merged_cluster_size_ > 0)
    {
        auto tt0 = Clock::now();
        GetClusterCenters();
        MergeClusters();
        auto tt1 = Clock::now();
        double merge_time = std::chrono::nanoseconds(tt1 - tt0).count()/1e9;
        total_time += merge_time;
        printf("MergeClusters time : %f ! \n", merge_time);
        // printf("merged cluster pair number : %d \n", merged_cluster_size_);
        // printf("1 cluster_cores_mat_ no zero: %d \n", cluster_cores_mat_.n_nonzero);

        auto tt2 = Clock::now();
        UpdateClusterNormals();
        auto tt3 = Clock::now();
        double update_normal_time = std::chrono::nanoseconds(tt3 - tt2).count()/1e9;
        printf("update_normal time : %f ! \n", update_normal_time);
        
        total_time += update_normal_time;
        // printf("cluster_adjacent_mat_ number : %d \n", cluster_adjacent_mat_.n_rows);
        auto tt4 = Clock::now();
        UpdateClusterScoreMat();
        CalculateClusterScores();
        auto tt5 = Clock::now();
        double update_score_time = std::chrono::nanoseconds(tt5 - tt4).count()/1e9;
        printf("update_score time : %f ! \n", update_score_time);

        total_time += update_score_time;

        std::string init_ptn_path_iter = out_dir_ + filename_ + "_ptn" + std::to_string(iter_num);
        OuputPtN(init_ptn_path_iter);
        if(merged_cluster_size_ == 0)
        {
            flipClusterNormalsByScores();
            init_ptn_path_iter = out_dir_ + filename_ + "_flipped" + std::to_string(iter_num);
            OuputPtN(init_ptn_path_iter);
        }


        std::string init_ptn_path_iter2 = out_dir_ + filename_ + "_orient_ptn" + std::to_string(iter_num);
        OuputPtN(init_ptn_path_iter2, true);

        
        printf("iter %d merged cluster number : %d \n", iter_num, merged_cluster_size_);
        iter_num ++;
        // break;
    }

    printf("total time used : %f \n", total_time);
}
