#include "local_vipss.hpp"
#include <cmath>
#include <algorithm>
#include "readers.h"  
#include "orient_normal.h"
#include <chrono>
#include <queue>

typedef std::chrono::high_resolution_clock Clock;

void LocalVipss::TestInsertPt()
{
    // std::string in_pt_path = "/home/jjxia/Documents/projects/VIPSS_LOCAL/data/doghead_100/doghead_100.ply";
    std::string in_pt_path =  "../../data/dragon/dragon_sample1k_noise.ply";
    std::vector<double> vts;
    std::vector<double> vns;
    readPLYFile(in_pt_path, vts, vns);
    // std::vector<tetgenmesh::point> insert_pts;
    auto t0 = Clock::now();
    for(size_t i = 0; i < vts.size() / 3; ++i)
    {
        tetgenmesh::point pt = new double[3];
        pt[0] = vts[3 * i];
        pt[1] = vts[3 * i + 1];
        pt[2] = vts[3 * i + 2];
        // insert_pts_.push_back(pt);
        // auto v_gen = voro_gen_;
        // voro_gen_.InsertPt(pt);
        printf("insert id : %lu \n", i);
    }

    auto t1 = Clock::now();
    double insert_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf(" point insert time total: %g \n", insert_time);

    // for(auto &pt : insert_pts_)
    // {
    //     voro_gen_.InsertPt(pt);
    // }
}

void LocalVipss::TestVoronoiPts()
{
    std::string in_pt_path =  "../../data/doghead_100/doghead_100_sample.ply";
    // std::string in_pt_path =  "../../data/dragon/dragon_sample.ply";
    // std::string in_pt_path =  "../../data/arma_100k/arma_100k_sample.ply";
    std::vector<double> vts;
    std::vector<double> vns;
    readPLYFile(in_pt_path, vts, vns);
    // std::vector<tetgenmesh::point> insert_pts;
    auto t000 = Clock::now();
    // voro_gen_.BuildTetMeshTetCenterMap();
    // voro_gen_.BuildPicoTree();
    auto t001 = Clock::now();
    double build_tree_time = std::chrono::nanoseconds(t001 - t000).count()/1e9;
    printf(" build_tree_time total: %g \n", build_tree_time);
    // std::vector<size_t> ids = {5, 5};
    // for(auto pt : voro_gen_.points_)
    // for(size_t i = 0; i < vts.size()/3; ++i)
    size_t pnei_count = 0;
    auto t0 = Clock::now();
    size_t search_pt_num = vts.size()/3;
// #pragma omp parallel for
    for(size_t i = 0; i < search_pt_num; ++i)
    {
        // count ++;
        // if(i > 5) continue;
        // tetgenmesh::point pt = new double[3];
        tetgenmesh::point pt = &vts.data()[3*i];
        // pt[0] = vts[3*i]; pt[1] = vts[3*i + 1]; pt[2] = vts[3*i + 2];
        // printf("finish get vertex star pts! \n");
        std::vector<tetgenmesh::point> nei_pts;
        voro_gen_.GetVoronoiNeiPts(pt, nei_pts);
        printf(" nei pt size : %ld \n", nei_pts.size());
        // auto t00 = Clock::now();
        
        for(auto nei_pt : nei_pts)
        {
            // double v1 = voro_gen_.CalTruncatedCellVolume(pt, nei_pt);
            // printf("v1 : %f \n", v1);
            // double v2 = voro_gen_.CalTruncatedCellVolumePass(pt, nei_pt);
            // printf("v2 : %f \n", v2);
            // double v3 = voro_gen_.CalUnionCellVolume(pt, nei_pt);
            // printf("v3 : %f \n", v3);
            // break;
        }
        pnei_count += nei_pts.size();
        // auto t11 = Clock::now();
        // double cal_volume_time = std::chrono::nanoseconds(t11 - t00).count()/1e9;
        // printf(" cal_volume_time time average: %g \n", cal_volume_time/double(nei_pts.size()));
        // printf(" cal_volume_time time total: %g \n", cal_volume_time);

        // size_t tet_count = voro_gen_.tetMesh_.tetrahedrons->items;
        // printf(".......... tet_count : %ld \n", tet_count);
        // ofstream file2;
        // std::string my_file2 = "../out/voro_star_" + std::to_string(i) + ".obj";
        // file2.open(my_file2);
        // file2 << "v " << pt[0] << " " <<pt[1] << " " << pt[2] << std::endl;
        // for(auto p2 : nei_pts)
        // {
        //     file2 << "v " << p2[0] << " " <<p2[1] << " " << p2[2] << std::endl;
        // }
        // file2.close();
        // printf("finish get voronoi nei pts! \n");
        // printf("insert id : %lu \n", i);
        // delete []pt;
        // break;
    }
    
    auto t1 = Clock::now();
    double total_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("cal nei voro cell truncated volume time total: %g \n", total_time);
    printf("cal nei voro cell truncated volume per pt time: %g \n", total_time/search_pt_num);
    printf("cal nei voro cell truncated volume per pair time: %g \n", total_time/(pnei_count));

}

void LocalVipss::VisualFuncValues(double (*function)(const R3Pt &in_pt), const VoroPlane& in_plane,
                              const std::string& dist_fuc_color_path)
{
    
    arma::vec3 ax_1 = {1, 1, 0}; 
    ax_1[2] = -(in_plane.nx * ax_1[0] + in_plane.ny * ax_1[1]) / in_plane.nz;
    double len = sqrt(ax_1[0] * ax_1[0] + ax_1[1] * ax_1[1] + ax_1[2] * ax_1[2]);
    ax_1[0] /= len;
    ax_1[1] /= len;
    ax_1[2] /= len;
    arma::vec3 pv{in_plane.nx, in_plane.ny, in_plane.nz};
    arma::vec3 ax_2 = arma::cross(ax_1, pv);

    std::ofstream out_file;
    std::string out_path = "../out/bi_intersect_plane.obj";
    out_file.open(out_path);

    arma::vec3 ori{in_plane.px, in_plane.py, in_plane.pz};
    double step = 0.001;
    int dim = 200;
    std::vector<double> pts;
    std::vector<uint8_t> pts_co;
    //rbf_core_.isHermite = true;
    //cout << "is hermite " << rbf_core_.isHermite << endl;
    
    for (int i = -dim; i < dim + 1; ++i)
    {
        for (int j = -dim; j < dim + 1; ++j)
        {
            arma::vec3 new_pt = ori + i * step * ax_1 + j * step * ax_2;
            
            pts.push_back(new_pt[0]);
            pts.push_back(new_pt[1]);
            pts.push_back(new_pt[2]);
            R3Pt pt(new_pt[0], new_pt[1], new_pt[2]);
            double dist = function(pt);
            double scale = 0.01;
            dist = dist > -scale ? dist : -scale;
            dist = dist < scale ? dist : scale;
            dist = dist / scale;
            int co_val = abs(dist) * 255;
            uint8_t c_val = uint8_t(co_val);
            if (dist >= 0)
            {
                pts_co.push_back(c_val);
                pts_co.push_back(0);
                pts_co.push_back(0);
            }
            else {
                pts_co.push_back(0);
                pts_co.push_back(c_val);
                pts_co.push_back(0);
            }
        }
    }
    printf("pts size : %ld , pts co size : %ld \n", pts.size()/3, pts_co.size()/3);
    writePLYFile_CO(dist_fuc_color_path, pts, pts_co);
}

void LocalVipss::Init(const std::string & path)
{   
    std::vector<double> in_pts;
    std::vector<double> in_normals;
    readPLYFile(path, in_pts, in_normals);
    printf("load data file : %s \n", path.c_str());
    printf("read point size : %lu \n", in_pts.size()/3);
    auto t0 = Clock::now();
    // printf("read point size : %d \n", in_pts.size());
    // voro_gen_.loadData(path);
    voro_gen_.out_dir_ = out_dir_;

    voro_gen_.loadData(in_pts);
    printf("start to init mesh \n");
    voro_gen_.InitMesh();
    // TestInsertPt();
    // TestVoronoiPts();
    // voro_gen_.BuildPtIdMap();
    voro_gen_.BuildAdjecentMat();
    printf("finish build adjecent mat \n");

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
    vipss_api_.user_lambda_ = user_lambda_;
    vipss_api_.n_voxel_line_ = volume_dim_;

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
    arma::sp_imat cluster_core_ids(cluster_cores_mat_.col(cluster_id));
    arma::sp_imat::const_iterator start = cluster_core_ids.begin();
    arma::sp_imat::const_iterator end = cluster_core_ids.end();
    std::vector<size_t> core_ids;
    for ( auto i = start; i != end; ++i )
    {
        core_ids.push_back(i.row());
    }
    return core_ids;
}


void LocalVipss::BuidClusterCoresPtIds()
{
    cluster_core_pt_ids_.clear();
    // size_t c_num = cluster_cores_mat_.n_rows;
    size_t c_num = cluster_cores_mat_.n_cols;
    for(size_t i =0; i < c_num; ++i)
    {
        auto cur_cores = GetClusterCoreIds(i);
        cluster_core_pt_ids_.push_back(cur_cores);
    }
}

void LocalVipss::UpdateClusterCoresPtIds()
{
    size_t c_num = cluster_cores_mat_.n_cols;
    size_t cur_n = cluster_core_pt_ids_.size();
    if(c_num <= cur_n) return; 
    for(size_t i = cur_n; i < c_num; ++i)
    {
        auto p_ids = GetClusterCoreIds(i);
        cluster_core_pt_ids_.push_back(p_ids);
    }
}

void LocalVipss::GetInitClusterPtIds(size_t cluster_id, 
        std::vector<double>& pts, std::vector<size_t>& pt_ids)
{
    auto& cluster_pts_map = voro_gen_.point_cluster_pts_map_;
    auto& cluster_pt_id_map = voro_gen_.point_id_map_;
    // auto& points = voro_gen_.points_;
    // printf("cluster id %d, points size %d \n ", cluster_id, points_.size()/3);
    tetgenmesh::point cur_ptr = points_[cluster_id];
    auto& cluster_pts = voro_gen_.point_cluster_pts_map_[cur_ptr];

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

    std::vector<size_t> pt_ids;
    std::set<size_t> key_ids;
    arma::sp_imat cluster_core_ids(cluster_cores_mat_.col(cluster_id));
    arma::sp_imat::const_iterator c_start = cluster_core_ids.begin();
    arma::sp_imat::const_iterator c_end = cluster_core_ids.end();
    for (auto i = c_start; i != c_end; ++i)
    {
        pt_ids.push_back(i.row());
        key_ids.insert(i.row());
    }
    arma::sp_imat cluster_pt_ids (cluster_adjacent_pt_mat_.col(cluster_id));
    arma::sp_imat::const_iterator start = cluster_pt_ids.begin();
    arma::sp_imat::const_iterator end = cluster_pt_ids.end();
    for ( auto i = start; i != end; ++i )
    {
        if(key_ids.find(i.row()) != key_ids.end()) continue;
        pt_ids.push_back(i.row());
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
    const arma::sp_imat cluster_col(cluster_cores_mat_.col(pid));
    arma::sp_imat::const_iterator start = cluster_col.begin();
    return (size_t)start.row();
}

void LocalVipss::CalculateClusterNormals(size_t cluster_id)
{
    auto cluster_pt_ids = GetClusterPtIds(cluster_id);
    // printf("cluster_pt_ids size : %zu \n", cluster_pt_ids.size());
    auto cluster_vts = GetClusterVerticesFromIds(cluster_pt_ids);
    // printf("cluster pt size : % \n", cluster_vts.size()/3);
    auto t1 = Clock::now();
    // vipss_api_.run_vipss(cluster_vts, cluster_pt_ids.size());
    arma::sp_imat cluster_core_ids(cluster_cores_mat_.col(cluster_id));
    size_t key_ptn = cluster_core_ids.n_nonzero;
    vipss_api_.run_vipss(cluster_vts, key_ptn);
    auto t2 = Clock::now();
    double vipss_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
    vipss_time_stats_.emplace_back(std::pair<size_t, double>(cluster_pt_ids.size(), vipss_time));

    for(size_t p_id = 0; p_id < cluster_pt_ids.size(); ++p_id)
    {
        size_t v_id = cluster_pt_ids[p_id];
        cluster_normal_x_(v_id, cluster_id) = vipss_api_.normals_[3* p_id];
        cluster_normal_y_(v_id, cluster_id) = vipss_api_.normals_[3* p_id + 1];
        cluster_normal_z_(v_id, cluster_id) = vipss_api_.normals_[3* p_id + 2];
    }
}

inline arma::sp_mat BuildClusterPtIdMatrix(const std::vector<size_t>& p_ids, size_t npt)
{
    size_t unit_npt = p_ids.size();
    arma::sp_mat unit_cluster_mat(unit_npt*4, npt * 4);
    for(size_t id = 0; id < unit_npt; ++id)
    {
        size_t pid = p_ids[id];
        unit_cluster_mat(id, pid) = 1.0;
        unit_cluster_mat(id + unit_npt,     pid + npt ) = 1.0;
        unit_cluster_mat(id + unit_npt * 2, pid + npt * 2) = 1.0;
        unit_cluster_mat(id + unit_npt * 3, pid + npt * 3) = 1.0;
    }
    return unit_cluster_mat;
}

void LocalVipss::SampleClusterPts()
{

}

void LocalVipss::AddClusterHMatrix(const std::vector<size_t>& p_ids, const arma::mat& J_m, size_t npt)
{
    size_t unit_npt = p_ids.size();
    // printf("uint pt num : %lu \n", unit_npt);
    auto& final_H = final_H_;
    // #pragma omp parallel shared(final_H)
    // arma::sp_mat final_H(npt * 4, npt * 4); 
    double sum = 0;
    // std::vector< std::unordered_map<int,double>> temp_map_values;
    // temp_map_values.resize(npt * 4);
    // #pragma omp parallel for
    // arma::umat locations(2, unit_npt * 4 * unit_npt * 4);
    // arma::vec values(unit_npt * 4 * unit_npt * 4);

    size_t e_id = 0;
    long col_num = 4 * npt;
    std::vector<Triplet> coefficients;
    for(size_t i = 0; i < unit_npt; ++i)
    {
        size_t pi = p_ids[i];
        for(size_t step = 0; step < 4; ++step)
        {
            size_t row = pi + npt * step;
            // auto& h_row_map = temp_Hmat_values_[row];
            // #pragma omp parallel for
            size_t j_row = i + unit_npt* step;
            for(size_t j = 0; j < unit_npt; ++j)
            {
                size_t pj = p_ids[j];
                // h_row_map[pj]           += J_m(j_row, j);
                // h_row_map[npt + pj]     += J_m(j_row, j + unit_npt);
                // h_row_map[npt * 2 + pj] += J_m(j_row, j + unit_npt * 2);
                // h_row_map[npt * 3 + pj] += J_m(j_row, j + unit_npt * 3);

                // for(size_t k = 0; k < 4; ++k)
                // {
                //     h_ele_triplets_.push_back(std::move(Triplet(row, pj + npt * k, J_m(j_row, j + unit_npt * k))));
                // }
                
                // long id0 = row * col_num + pj;
                // h_pos_value_map_[id0] += J_m(j_row, j);
                // long id1 = row * col_num + npt + pj ;
                // h_pos_value_map_[id1] += J_m(j_row, j + unit_npt);

                // long id2 = row * col_num + npt * 2 + pj ;
                // h_pos_value_map_[id2] += J_m(j_row, j + unit_npt * 2);

                // long id3 = row * col_num + npt * 3 + pj ;
                // h_pos_value_map_[id3] += J_m(j_row, j + unit_npt * 3);

                // for(size_t k = 0; k < 4; ++k)
                // {
                //     locations(0, e_id) = row;
                //     locations(1, e_id) = pj + npt * k;
                //     values(e_id)       = J_m(j_row, j + unit_npt * k);
                //     e_id ++;
                // }
                
                final_H(row, pj)           += J_m(j_row, j);
                final_H(row, npt + pj)     += J_m(j_row, j + unit_npt);
                final_H(row, npt * 2 + pj) += J_m(j_row, j + unit_npt * 2);
                final_H(row, npt * 3 + pj) += J_m(j_row, j + unit_npt * 3);
                
                // final_H(pi + npt * step, pj)           += J_m(i + unit_npt* step, j);
                // final_H(pi + npt * step, npt + pj)     += J_m(i + unit_npt* step, j + unit_npt);
                // final_H(pi + npt * step, npt * 2 + pj) += J_m(i + unit_npt* step, j + unit_npt * 2);
                // final_H(pi + npt * step, npt * 3 + pj) += J_m(i + unit_npt* step, j + unit_npt * 3);
            }
        }
    }

    // SpMat cur_H(4 * npt, 4 * npt);
    // cur_H.setFromTriplets(coefficients.begin(), coefficients.end());
    // final_h_eigen_ += cur_H;
    // auto cur_H = arma::sp_mat(locations, values, 4 * npt, 4 * npt);
    // final_H_ =  cur_H + final_H_;
    // temp_H_vec_.push_back(std::move(final_H));
    // printf("----final h nozero: %llu \n", final_H_.n_nonzero);
    // printf("----J  ele number: %llu \n", J_m.n_elem);
    // printf("----final h sum: %lf \n", arma::accu(final_H_));
    // printf("----J  sum: %lf \n", arma::accu(J_m));
}

void LocalVipss::BuildCurrentH(const arma::sp_mat& unit_cluster_mat, const arma::sp_mat&J_m, double lambda)
{
        // const auto& unit_cluster_mat = cluster_pt_mat_vec[i];
        // const auto& J_m = cluster_J_mat_vec[i];
        // const auto& current_ids = cluster_pt_ids_vec[i];
        if(lambda < 1e-10)
        {            
            arma::sp_mat sp_J(J_m); 
            arma::sp_mat temp_H = unit_cluster_mat.t() * sp_J * unit_cluster_mat;
            temp_H_vec_.push_back(std::move(temp_H));
            // final_H_temp_  += temp_H;
            // auto delt_sum = arma::accu(final_H_temp_ - final_H_);
            // printf("delt_sum : %lf \n", delt_sum);
        } else {
            size_t unit_npt = unit_cluster_mat.n_rows / 4; 
            arma::mat F(4 * unit_npt, 4 * unit_npt);
            arma::mat E(unit_npt, unit_npt);
            E.eye();
            F(0, 0, arma::size(unit_npt, unit_npt)) = E;
            arma::mat K = (F + J_m * (lambda));
            arma::sp_mat sp_K(K);  
            arma::sp_mat temp_H = unit_cluster_mat.t() * sp_K * unit_cluster_mat;
            temp_H_vec_.push_back(std::move(temp_H));
            // final_H_  += temp_H;
        }
}

void LocalVipss::CalHVecSum()
{
    size_t cur_n = temp_H_vec_.size();
    size_t H_n = temp_H_vec_.size();
    size_t cur_loop = 0;
    // printf("npt_ %lu log 2 : %d \n", npt_, (int)log2(npt_));
    size_t loop_level = (size_t)log2(cur_n);
    if(cur_n > pow(2, loop_level)) loop_level ++;
    for(size_t cur_loop = 1; cur_loop <= loop_level; ++cur_loop)
    {
        size_t cur_step = pow(2, cur_loop);
        
        if(cur_n / 2 == 0) break;
        cur_n = cur_n / 2 + cur_n % 2; 
        // printf("current N 0: %llu \n", cur_n);
    // #pragma omp parallel for
        for(size_t m_i = 0; m_i < cur_n; ++m_i)
        {
            if(m_i*cur_step + cur_step/2 < H_n)
            {
                temp_H_vec_[m_i*cur_step] += temp_H_vec_[m_i*cur_step + cur_step/2]; 
            }
        } 
    }
    final_H_ += temp_H_vec_[0];
    temp_H_vec_.clear();
}

void LocalVipss::BuildMatrixH()
{
    auto t00 = Clock::now(); 
    size_t npt = this->points_.size();
    size_t cluster_num = cluster_cores_mat_.n_cols;
    final_H_.resize(4 * npt , 4 * npt);
    // final_H_ = arma::speye<arma::sp_mat>(4 * npt , 4 * npt) * 1e-16; 
    double time_sum = 0;
    temp_Hmat_values_.resize(4 * npt);
    // temp_Hmat_values_.eye();
    final_h_eigen_.resize(4 * npt , 4 * npt);
    // h_ele_triplets_.reserve(cluster_num * 300 * 300);
    for(size_t i =0; i < cluster_num; ++i)
    {
        // std::vector<double> vts;
        // std::vector<size_t> p_ids;
        auto cluster_pt_ids = GetClusterPtIds(i);
        auto cluster_pt_vec = GetClusterVerticesFromIds(cluster_pt_ids);
        auto t3 = Clock::now();
        vipss_api_.build_unit_vipss(cluster_pt_vec);
        time_sum += vipss_api_.rbf_core_.bigM_inv_time;
        auto t4 = Clock::now();
        double build_j_time = std::chrono::nanoseconds(t4 - t3).count()/1e9;
        build_j_time_total_ += build_j_time;
        // printf("build_j_time_total_ : %f \n", build_j_time_total_);
        // cluster_J_mat_vec_.push_back(std::move(vipss_api_.rbf_core_.Minv));
        auto t5 = Clock::now();
        if(user_lambda_ > 1e-10)
        {
            size_t unit_npt = cluster_pt_ids.size(); 
            arma::mat F(4 * unit_npt, 4 * unit_npt);
            arma::mat E(unit_npt, unit_npt);
            E.eye();
            F(0, 0, arma::size(unit_npt, unit_npt)) = E;
            // double cur_lambda = user_lambda_ * double(unit_npt) / double(npt);
            double cur_lambda = user_lambda_;
            arma::mat K = (F + vipss_api_.rbf_core_.Minv * cur_lambda);
            AddClusterHMatrix(cluster_pt_ids, K, npt);
        } else {
            AddClusterHMatrix(cluster_pt_ids, vipss_api_.rbf_core_.Minv, npt);
        }
        // auto unit_pt_mat = BuildClusterPtIdMatrix(p_ids, cluster_num);
        // arma::sp_mat new_J(vipss_api_.rbf_core_.Minv);
        // BuildCurrentH(unit_pt_mat, new_J, user_lambda_);  
        auto t6 = Clock::now();
        double build_h_time = std::chrono::nanoseconds(t6 - t5).count()/1e9;
        build_h_time_total_ += build_h_time;
    }
    auto t_h1 = Clock::now();
    // size_t ele_num = h_pos_value_map_.size();
    // arma::umat locations(2, ele_num);
    // arma::vec values(ele_num);

    // size_t e_id = 0;
    // for(const auto& ele : h_pos_value_map_)
    // {
    //     locations(0,e_id) = ele.first / (4 * npt);
    //     locations(1,e_id) = ele.first % (4 * npt);
    //     values(e_id)      = ele.second;
    //     e_id++;
    // }

    // std::vector<size_t> rows_vec;
    // std::vector<size_t> cols_vec;
    // std::vector<double> val_vec;
    
    // std::vector<Triplet> coefficients;
    // for(size_t i = 0; i < 4 * npt; ++i)
    // {
    //     const auto& row = temp_Hmat_values_[i];
    //     for(const auto& ele : row)
    //     {
    //         coefficients.push_back(Triplet(i, ele.first, ele.second));
    //         // rows_vec.push_back(i);
    //         // cols_vec.push_back(ele.first);
    //         // val_vec.push_back(ele.second);
    //     }
    // }
    // size_t ele_num = coefficients.size();
    // arma::umat locations(2, ele_num);
    // arma::vec values(ele_num);
    // for(size_t i = 0; i < ele_num; ++i)
    // {
    //     locations(0,i) = coefficients[i].row();
    //     locations(1,i) = coefficients[i].col();
    //     values(i)      = coefficients[i].value();
    // } 

    // auto t_h1 = Clock::now();
    // final_H_ = arma::sp_mat(locations, values, 4* npt, 4 * npt);
    // printf("------- X : %llu, %llu \n", X.n_rows, X.n_cols);
    // final_h_eigen_.setFromTriplets(h_ele_triplets_.begin(), h_ele_triplets_.end());
    auto t_h2 = Clock::now();
    double build_h_time = std::chrono::nanoseconds(t_h2 - t_h1).count()/1e9;
    printf("--- build_H_time from map : %f \n", build_h_time);
    build_h_time_total_ += build_h_time;

    auto t11 = Clock::now();
    double build_H_time_total = std::chrono::nanoseconds(t11 - t00).count()/1e9;

    printf("--- build_j_time_total_  time : %f \n", build_j_time_total_);
    printf("--- add_h_time_total_  time : %f \n", build_h_time_total_);
    printf("--- build_H_time_total sum  time : %f \n", build_H_time_total);
    printf("--- big M inv sum  time : %f \n", time_sum);

}

inline std::vector<double> LocalVipss::GetClusterNormalsFromIds
                            (const std::vector<size_t>& pt_ids, const std::vector<double>& all_normals) const
{
    std::vector<double> normals;
    for(size_t id : pt_ids)
    {
        // auto &pt = points_[id];
        normals.push_back(all_normals[3 * id]);
        normals.push_back(all_normals[3 * id + 1]);
        normals.push_back(all_normals[3 * id + 2]);
    }
    return normals;
}

inline std::vector<double> LocalVipss::GetClusterSvalsFromIds(const std::vector<size_t>& pt_ids, 
                        const std::vector<double>& all_svals) const
{
    std::vector<double> svals;
    for(size_t id : pt_ids)
    {
        svals.push_back(all_svals[id]);
    }
    return svals;
}

void LocalVipss::BuildHRBFPerNode()
{
    auto t00 = Clock::now(); 
// if(0)

    size_t npt = this->points_.size();
    size_t cluster_num = cluster_cores_mat_.n_cols;
    double time_sum = 0;
    node_rbf_vec_.resize(cluster_num);
    bool use_partial_vipss = false;
    vipss_api_.user_lambda_ = user_lambda_;

    printf("cluster num : %lu \n", cluster_num);

    for(size_t i =0; i < cluster_num; ++i)
    {
        std::vector<double> cluster_pt_vec;
        std::vector<size_t> cluster_pt_ids;
        GetInitClusterPtIds(i, cluster_pt_vec, cluster_pt_ids);
        std::vector<double> cluster_nl_vec;
        if(use_partial_vipss) 
        {
            // cluster_nl_vec = {normals_[3*i], normals_[3*i + 1], normals_[3*i + 2]};
        } else {
            cluster_nl_vec = GetClusterNormalsFromIds(cluster_pt_ids, normals_);
        }
        auto cluster_sv_vec = GetClusterSvalsFromIds(cluster_pt_ids, s_vals_);

        // std::string cluster_pt_path = out_dir_ + "cluster_pts/" + std::to_string(i) + "_cluster_color"; 
        // // std::cout << " Cluster pts : " << cluster_pt_path << std::endl;
        // output_opt_pts_with_color(cluster_pt_vec, cluster_sv_vec, cluster_pt_path);
        // printf("cluster_sv_vec size : %ld \n", cluster_sv_vec.size());
        node_rbf_vec_[i] = std::make_shared<RBF_Core>();
        vipss_api_.build_cluster_hrbf(cluster_pt_vec, cluster_nl_vec, cluster_sv_vec, node_rbf_vec_[i]);
    }

    auto t11 = Clock::now();
    double build_HRBF_time_total = std::chrono::nanoseconds(t11 - t00).count()/1e9;
    printf("--- build_HRBF_time_total sum  time : %f \n", build_HRBF_time_total);
}

double LocalVipss::NodeDistanceFunction(const tetgenmesh::point nn_pt, const tetgenmesh::point cur_pt)
{
    if(voro_gen_.point_id_map_.find(nn_pt) != voro_gen_.point_id_map_.end())
    {
        size_t pid = voro_gen_.point_id_map_[nn_pt];
        return node_rbf_vec_[pid]->Dist_Function(cur_pt);
    } 
    return 0;
}

void LocalVipss::testNNPtDist()
{
    // printf("start to test nn dist -------- \n");
    // double test_pt[3] = {0.301703, 0.34657, 0.0822563};
    double test_pt[3] =  {-0.143955, 0.149453, -0.273193};
    std::vector<tetgenmesh::point> nei_pts;
    voro_gen_.GetVoronoiNeiPts(test_pt, nei_pts);
    std::vector<double> nn_pts;
    std::vector<double> nn_normals;
    int n_voxel_line_ = 50;
    std::string sphere_path = "/home/jianjun/Documents/sphere.ply";
    std::vector<std::array<double, 3>> vts;
    std::vector<std::vector<size_t>> faces;
    // printf("start to readPlyMesh -------- \n");
    readPlyMesh(sphere_path, vts, faces);
    // printf("finish readPlyMesh -------- \n");
    // readPLYFile()
    size_t nn_num = nei_pts.size();
    // arma::vec dist_vec(nn_num);
    double volume_sum = 0;
    double weight_sum = 0;
    
    for(size_t i = 0; i < nn_num; ++i)
    {
        auto nn_pt = nei_pts[i];
        //  printf("start to cal node distance -------- \n");
        if(nn_pt != (tetgenmesh::point)NULL && voro_gen_.tetMesh_.pointtype(nn_pt) != tetgenmesh::UNUSEDVERTEX)
        {
            double nn_dist = NodeDistanceFunction(nn_pt, test_pt);
            // printf("nn dist %f \n", nn_dist);
            double volume  = voro_gen_.CalTruncatedCellVolumePass(test_pt, nn_pt);
            // double volume2  = voro_gen_.CalTruncatedCellVolume(test_pt, nn_pt);
            // printf("nn dist %f, volume %f  volume2 %f  \n", nn_dist, volume, volume2);
            if(volume < 0.01) continue;
            weight_sum += nn_dist * volume; 
            volume_sum += volume;
            size_t p_id = voro_gen_.point_id_map_[nn_pt];
            std::string out_path = out_dir_ + "/nn_spheres/" + std::to_string(p_id) + "_sphere" + std::to_string(volume) + ".ply";
            std::array<double,3> center = {nn_pt[0], nn_pt[1], nn_pt[2]};
            SaveSphere(out_path, vts, faces, center, volume);
        }
    }

    
    for(auto pt : nei_pts)
    {
        if(voro_gen_.tetMesh_.pointtype(pt) != tetgenmesh::UNUSEDVERTEX)
        {
            nn_pts.push_back(pt[0]);
            nn_pts.push_back(pt[1]);
            nn_pts.push_back(pt[2]);
            size_t p_id = voro_gen_.point_id_map_[pt];

            nn_normals.push_back(normals_[3*p_id]);
            nn_normals.push_back(normals_[3*p_id + 1]);
            nn_normals.push_back(normals_[3*p_id + 2]);

            // size_t p_id = voro_gen_.point_id_map_[pt];
            if(0)
            {
                node_rbf_vec_[p_id]->Surfacing(0,n_voxel_line_);
                node_rbf_vec_[p_id]->Write_Surface(out_dir_ + "cluster_mesh/"+  std::to_string(p_id) + "_surface");
            }
            
        }
        
    }

    std::ofstream out_nn_file;
    std::string nn_pt_path = out_dir_  + "nn_pts";
    writePLYFile_VN(nn_pt_path, nn_pts, nn_normals);
    writeXYZ(nn_pt_path, nn_pts);
}

double LocalVipss::NatureNeighborDistanceFunction(const tetgenmesh::point cur_pt)
{
    std::vector<tetgenmesh::point> nei_pts;
    voro_gen_.GetVoronoiNeiPts(cur_pt, nei_pts);
    // printf("nn num %ld \n", nei_pts.size());
    size_t nn_num = nei_pts.size();
    // arma::vec dist_vec(nn_num);
    double volume_sum = 0;
    double weight_sum = 0;
    // printf("nn num %ld \n", nn_num);

    for(size_t i = 0; i < nn_num; ++i)
    {
        auto nn_pt = nei_pts[i];
        if(nn_pt != (tetgenmesh::point)NULL && voro_gen_.tetMesh_.pointtype(nn_pt) != tetgenmesh::UNUSEDVERTEX)
        {
            double nn_dist = NodeDistanceFunction(nn_pt, cur_pt);
            // printf("nn dist %f \n", nn_dist);
            double volume  = voro_gen_.CalTruncatedCellVolumePass(cur_pt, nn_pt);
            // double volume  = voro_gen_.CalTruncatedCellVolume(cur_pt, nn_pt);
            // double volume  = voro_gen_.CalUnionCellVolume(cur_pt, nn_pt);
            // printf("nn dist %f, volume %f \n", nn_dist, volume);
            weight_sum += nn_dist * volume; 
            volume_sum += volume;
        }
    }

    if(volume_sum != 0)
    {
        double weight_dist = weight_sum / volume_sum;
        return weight_dist;
    }
    return 0;
}


void LocalVipss::InitNormalWithVipss()
{
    size_t cluster_num = cluster_cores_mat_.n_cols;
    // printf("cluster num : %zu \n", cluster_num);
    vipss_time_stats_.clear();
    double vipss_sum = 0;
    auto t_init = Clock::now();
    size_t npt = this->points_.size();

    final_H_.resize(4 * npt , 4 * npt);
    // final_H_temp_.resize(4 * npt , 4 * npt);
    for(size_t i =0; i < cluster_num; ++i)
    {
        std::vector<double> vts;
        std::vector<size_t> p_ids;
        GetInitClusterPtIds(i, vts, p_ids);
        auto t1 = Clock::now();
        size_t unit_npt = p_ids.size(); 
        double cur_lambda = user_lambda_ / double(unit_npt);
        // double cur_lambda = user_lambda_ * double(unit_npt) / double(npt);
        // double cur_lambda = user_lambda_;
        vipss_api_.user_lambda_ = cur_lambda;
        vipss_api_.run_vipss(vts, 1);
        // printf("finish partial vipss 22\n");
        auto t2 = Clock::now();
        double vipss_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
        vipss_sum += vipss_time;
        vipss_time_stats_.emplace_back(std::pair<size_t, double>(p_ids.size(), vipss_time));

        // printf("finish partial vipss 33\n");
        for(size_t p_id = 0; p_id < p_ids.size(); ++p_id)
        {
            size_t v_id = p_ids[p_id];
            cluster_normal_x_(v_id, i) = vipss_api_.normals_[3* p_id];
            cluster_normal_y_(v_id, i) = vipss_api_.normals_[3* p_id + 1];
            cluster_normal_z_(v_id, i) = vipss_api_.normals_[3* p_id + 2];
        }
    }
    auto t_init2 = Clock::now();
    double all_init_time = std::chrono::nanoseconds(t_init2 - t_init).count()/1e9;
    printf("all init time : %f \n", all_init_time);
    printf("all init vipss time : %f \n", vipss_sum);
    

    // cluster_ptn_vipss_time_stats_.push_back(vipss_time_stats_);
}

static LocalVipss* local_vipss_ptr;
int LocalVipss::DistCallNum = 0;
double LocalVipss::DistCallTime = 0.0;

double LocalVipss::NNDistFunction(const R3Pt &in_pt)
{
    DistCallNum ++;
    auto t0 = Clock::now();
    double new_pt[3] = {in_pt[0], in_pt[1], in_pt[2]};
    double dist = local_vipss_ptr->NatureNeighborDistanceFunction(&(new_pt[0]));
    auto t1 = Clock::now();
    double dist_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    DistCallTime += dist_time;
    return dist;
}

void LocalVipss::SetThis()
{
    local_vipss_ptr = this;
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
        if(cluster_normal_x_(id, c_b) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    for(auto id : core_ids_b)
    {
        if(cluster_normal_x_(id, c_a) != 0)
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
        normal_ma(i, 0) =  cluster_normal_x_(id, c_a);
        normal_ma(i, 1) =  cluster_normal_y_(id, c_a);
        normal_ma(i, 2) =  cluster_normal_z_(id, c_a);

        normal_mb(i, 0) =  cluster_normal_x_(id, c_b);
        normal_mb(i, 1) =  cluster_normal_y_(id, c_b);
        normal_mb(i, 2) =  cluster_normal_z_(id, c_b);
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
        if(cluster_normal_x_(id, c_b) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    for(auto id : core_ids_b)
    {
        if(cluster_normal_x_(id, c_a) != 0)
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
        normal_ma(i, 0) =  cluster_normal_x_(id, c_a);
        normal_ma(i, 1) =  cluster_normal_y_(id, c_a);
        normal_ma(i, 2) =  cluster_normal_z_(id, c_a);

        normal_mb(i, 0) =  cluster_normal_x_(id, c_b);
        normal_mb(i, 1) =  cluster_normal_y_(id, c_b);
        normal_mb(i, 2) =  cluster_normal_z_(id, c_b);
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
}

void LocalVipss::BuildClusterAdjacentMat()
{
    // printf("00 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    printf("adjacent no zero count : %llu \n", adjacent_mat_.n_nonzero);
    cluster_adjacent_mat_ = adjacent_mat_ * cluster_cores_mat_;
    cluster_adjacent_pt_mat_ = adjacent_mat_ * cluster_cores_mat_ + cluster_cores_mat_;
   
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
        if(cluster_scores_mat_(n_pos, i) != 0) 
        {
            continue;
        }
        cur_i_mat(1, 0) = cluster_normal_x_(n_pos, i);
        cur_i_mat(1, 1) = cluster_normal_y_(n_pos, i);
        cur_i_mat(1, 2) = cluster_normal_z_(n_pos, i);

        cur_n_mat(0, 0) = cluster_normal_x_(i, n_pos);
        cur_n_mat(0, 1) = cluster_normal_y_(i, n_pos);
        cur_n_mat(0, 2) = cluster_normal_z_(i, n_pos);

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
            // cluster_adjacent_flip_mat_(i, n_pos) = 1;
            // cluster_adjacent_flip_mat_(n_pos, i) = 1;
        }
        double score = std::min(1.0, std::max(-1.0, s1));
        score = acos(score) * Anlge_PI_Rate ;
        // printf("cluster id : %d \n", i);
        cluster_scores_mat_(n_pos, i) = score;
        cluster_scores_mat_(i, n_pos) = score;
    }
}


void LocalVipss::CalculateSingleClusterNeiScores(size_t i)
{
    arma::sp_imat nei(cluster_adjacent_mat_.col(i));
    const arma::sp_imat::const_iterator start = nei.begin();
    const arma::sp_imat::const_iterator end = nei.end();
    // printf("cluster id : %d \n", i);

    int count = 0;
    double sum = 0;
    for(auto iter = start; iter != end; ++ iter)
    {
        size_t n_pos = iter.row();
        if(i == n_pos) continue;
        count ++;
        if(cluster_scores_mat_(i, n_pos) != 0) 
        {
            continue;
        }
        bool flip_normal = false;
        double score = CalculateClusterPairScore(i, n_pos, flip_normal);
        // if(flip_normal)
        // {
        //     cluster_adjacent_flip_mat_(i, n_pos) = 1;
        //     cluster_adjacent_flip_mat_(n_pos, i) = 1;
        // }
        // printf("cluster id : %d \n", i);
        cluster_scores_mat_(i, n_pos) = score;
        cluster_scores_mat_(n_pos, i) = score;
        update_score_cluster_ids_.insert(n_pos);
        // sum += score;
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
    size_t c_num = cluster_scores_mat_.n_cols;

    // arma::sp_vec scores_vec(arma::sum(cluster_scores_mat_,0));
    // cluster_scores_vec_ = arma::sp_colvec(scores_vec);
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
}

void LocalVipss::GetClusterCenters()
{
    if(cluster_centers_.empty())
    {
        for(size_t i = 0; i < cluster_cores_mat_.n_cols; ++i)
        {
            tetgenmesh::point new_p = new double[3];
            new_p[0] =0;
            new_p[1] =0;
            new_p[2] =0;
            cluster_centers_.push_back(new_p);
        } 
    } else {
        for(size_t i = 0; i < cluster_cores_mat_.n_cols; ++i)
        {
            tetgenmesh::point& cur_p = cluster_centers_[i];
            cur_p[0] =0;
            cur_p[1] =0;
            cur_p[2] =0;
        } 
    }
    
    for(size_t i = 0; i < cluster_cores_mat_.n_cols; ++i)
    {
        tetgenmesh::point& new_p = cluster_centers_[i];
        const arma::sp_imat cores_pt(cluster_cores_mat_.col(i));
        const arma::sp_imat::const_iterator start = cores_pt.begin();
        const arma::sp_imat::const_iterator end = cores_pt.end();
        int count = 0;
        for(auto iter = start; iter != end; ++iter)
        {
            size_t pt_id = iter.row();
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
    std::vector<arma::uword> merged_cluster_ids;
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
    printf("cluster scores indices num : %llu \n", indices.size());
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
            const arma::sp_mat c_scores_vec(cluster_scores_mat_.col(id));
            const arma::sp_mat::const_iterator start = c_scores_vec.begin();
            const arma::sp_mat::const_iterator end = c_scores_vec.end();
            double max_score = 0;
            int max_id = -1;
            // auto cur_center = cluster_centers_[ele.first];
            for(auto iter = start; iter != end; ++iter)
            {
                if(visited_clusters.find(iter.row()) == visited_clusters.end())
                {
                    // double cur_dist = PtDistance(cluster_centers_[iter.internal_col], cur_center);
                    // cur_dist = cur_dist > 1e-10? cur_dist: 1e-10;
                    double cur_score = *iter;
                    if(cur_score > max_score && cluster_scores_vec_[iter.row()] > angle_threshold_)
                    {
                        max_id = iter.row();
                        max_score = cur_score;
                    }
                }
            }

            if(max_id != -1)
            {   
                visited_clusters.insert(max_id);
                merged_cluster_ids.push_back(id);
                merged_cluster_ids.push_back(max_id);

                // arma::sp_icolvec new_cluster_cores = cluster_cores_mat_.col(id) + cluster_cores_mat_.col(max_id);
                // AppendCol(cluster_cores_mat_, new_cluster_cores);

                // arma::sp_icolvec adjacent_pts_col = cluster_adjacent_pt_mat_.col(id) + cluster_adjacent_pt_mat_.col(max_id);
                // AppendCol(cluster_adjacent_pt_mat_, adjacent_pts_col);
                
                // arma::sp_irowvec adj_row = cluster_adjacent_mat_.row(id) + cluster_adjacent_mat_.row(max_id);
                // AppendRow(cluster_adjacent_mat_, adj_row);
                // arma::sp_icolvec adj_col = cluster_adjacent_mat_.col(id) + cluster_adjacent_mat_.col(max_id);
                // AppendCol(cluster_adjacent_mat_, adj_col);
            }
        }
    }

    auto t00 = Clock::now();
    size_t merged_pair_num = merged_cluster_ids.size() / 2;
    if(merged_pair_num > 0)
    {
        size_t c_num = cluster_cores_mat_.n_cols;
        size_t pt_num = cluster_cores_mat_.n_rows;

        arma::sp_imat new_cluster_cores_mat(pt_num, c_num + merged_pair_num); 
        new_cluster_cores_mat.cols(0, c_num-1) = cluster_cores_mat_.cols(0, c_num-1);
        cluster_cores_mat_ = new_cluster_cores_mat;    
        arma::sp_imat new_cluster_adjacent_pt_mat(pt_num, c_num + merged_pair_num);
        new_cluster_adjacent_pt_mat.cols(0, c_num-1) = cluster_adjacent_pt_mat_.cols(0, c_num-1);
        cluster_adjacent_pt_mat_ = new_cluster_adjacent_pt_mat;

        arma::sp_imat new_cluster_adjacent_mat(c_num + merged_pair_num, c_num + merged_pair_num);
        new_cluster_adjacent_mat(0, 0, arma::size(c_num, c_num)) = 
            cluster_adjacent_mat_(0, 0, arma::size(c_num, c_num));
        cluster_adjacent_mat_ = new_cluster_adjacent_mat;

        for(size_t i = 0; i < merged_pair_num; ++i)
        {
            size_t ca = merged_cluster_ids[2 * i];
            size_t cb = merged_cluster_ids[2 * i + 1];
            cluster_cores_mat_.col(c_num + i) = cluster_cores_mat_.col(ca) + cluster_cores_mat_.col(cb);
            cluster_adjacent_pt_mat_.col(c_num + i) = cluster_adjacent_pt_mat_.col(ca) + cluster_adjacent_pt_mat_.col(cb);
            cluster_adjacent_mat_.row(c_num + i) = cluster_adjacent_mat_.row(ca) + cluster_adjacent_mat_.row(cb);
            cluster_adjacent_mat_.col(c_num + i) = cluster_adjacent_mat_.col(ca) + cluster_adjacent_mat_.col(cb);
        // cluster_adjacent_mat_.insert
        }
    }
    auto t01 = Clock::now();
    double append_data_time = std::chrono::nanoseconds(t01 - t00).count()/1e9;
    printf("------finish append_data_time time : %f ! \n", append_data_time);
        
    merged_cluster_size_ = merged_cluster_ids.size();

    // auto t00v1 = Clock::now();
    std::sort(merged_cluster_ids.begin(), merged_cluster_ids.end(), std::greater<>());
    // auto t00v2 = Clock::now();
    // double scores_sort_time = std::chrono::nanoseconds(t00v2 - t00v1).count()/1e9;
    // printf("------finish scores_sort_time time : %f ! \n", scores_sort_time);

    auto t000 = Clock::now();
    // arma::uvec delete_ids(merged_cluster_ids);
    // // (merged_cluster_ids.begin(), merged_cluster_ids.end());
    // // arma::uvec delete_ids = merged_cluster_ids;
    // cluster_normal_x_.shed_cols(delete_ids);
    // cluster_normal_y_.shed_cols(delete_ids);
    // cluster_normal_z_.shed_cols(delete_ids);
    // cluster_cores_mat_.shed_cols(merged_cluster_ids);

    // ShedCols(cluster_cores_mat_, merged_cluster_ids);
    // ShedCols(cluster_adjacent_pt_mat_, merged_cluster_ids);
    // ShedCols(cluster_adjacent_mat_, merged_cluster_ids);

    // ShedCols(cluster_scores_mat_, merged_cluster_ids);
    // ShedCols(cluster_normal_x_, merged_cluster_ids);
    // ShedCols(cluster_normal_y_, merged_cluster_ids);
    // ShedCols(cluster_normal_z_, merged_cluster_ids);
    
    for(auto id : merged_cluster_ids)
    {
        cluster_cores_mat_.shed_col(id);
        cluster_adjacent_pt_mat_.shed_col(id);
        // cluster_scores_mat_.shed_row(id);
        cluster_scores_mat_.shed_col(id);
        // cluster_adjacent_mat_.shed_row(id);
        cluster_adjacent_mat_.shed_col(id);
        cluster_normal_x_.shed_col(id);
        cluster_normal_y_.shed_col(id);
        cluster_normal_z_.shed_col(id);
        cluster_core_pt_ids_.erase(std::next(cluster_core_pt_ids_.begin(), id));
    }

    
    // for(auto id : merged_cluster_ids)
    // {
        // cluster_core_pt_ids_.erase(std::next(cluster_core_pt_ids_.begin(), id));
    // }
    

    auto t00t1 = Clock::now();
    cluster_scores_mat_ = cluster_scores_mat_.t();
    cluster_adjacent_mat_ = cluster_adjacent_mat_.t();
    
    for(auto id : merged_cluster_ids)
    {
        cluster_scores_mat_.shed_col(id);
        cluster_adjacent_mat_.shed_col(id);
    }
    auto t00t2 = Clock::now();
    double update_mat_trans_time = std::chrono::nanoseconds(t00t2 - t00t1).count()/1e9;
    printf("------finish update update_mat_trans_time time : %f ! \n", update_mat_trans_time);
    // cluster_scores_mat_ = cluster_scores_mat_.t();
    // cluster_adjacent_mat_ = cluster_adjacent_mat_.t();

    auto t001 = Clock::now();
    double update_score_time = std::chrono::nanoseconds(t001 - t000).count()/1e9;
    printf("------finish update shed rows and cols time : %f ! \n", update_score_time);

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

inline void LocalVipss::ShedCols(arma::sp_imat& in_mat, const std::vector<arma::uword>& delete_ids)
{
    for(auto id : delete_ids)
    {
        in_mat.shed_col(id);
    }
}

inline void LocalVipss::ShedCols(arma::sp_mat& in_mat, const std::vector<arma::uword>& delete_ids)
{
    for(auto id : delete_ids)
    {
        in_mat.shed_col(id);
    }
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
    size_t c_num = cluster_cores_mat_.n_cols;
    size_t s_rows = cluster_scores_mat_.n_rows;
    size_t s_cols = cluster_scores_mat_.n_cols;

    // auto t0 = Clock::now();
    arma::sp_mat new_scores_mat(c_num, c_num);
    new_scores_mat(0, 0, arma::size(s_rows, s_cols)) = cluster_scores_mat_(0, 0, arma::size(s_rows, s_cols));
    cluster_scores_mat_ = new_scores_mat;

    // arma::sp_imat new_cluster_adjacent_flip_mat(c_num, c_num);
    // new_cluster_adjacent_flip_mat(0, 0, arma::size(s_rows -1, s_cols -1)) =
    //          cluster_adjacent_flip_mat_(0, 0, arma::size(s_rows -1, s_cols -1));
    // cluster_adjacent_flip_mat_ = new_cluster_adjacent_flip_mat;

    // auto t1 = Clock::now();
    // double update_score_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    // printf("finish update_score_time time : %f ! \n", update_score_time);

    
    // printf("cluster_adjacent_mat_ rows : %d , cols : %d \n", cluster_adjacent_mat_.n_rows, cluster_adjacent_mat_.n_cols);
    
    for(size_t i = s_rows; i < c_num; ++i)
    {
        CalculateSingleClusterNeiScores(i);
    }
    printf("cluster_scores_mat_ size %zu updated to %zu \n", s_rows, c_num);
}


void LocalVipss::UpdateClusterNormals()
{
    size_t c_num = cluster_cores_mat_.n_cols;
    size_t s_rows = cluster_normal_x_.n_rows;
    size_t s_cols = cluster_normal_x_.n_cols;

    // printf("cluster num : %zu , s_rows : %zu, s_cols: %zu \n", c_num, s_rows, s_cols);

    arma::sp_mat new_normal_x(s_rows, c_num);
    arma::sp_mat new_normal_y(s_rows, c_num);
    arma::sp_mat new_normal_z(s_rows, c_num);

    if(s_cols > 0)
    {
        new_normal_x.cols(0, s_cols-1) = cluster_normal_x_.cols(0, s_cols-1);
        new_normal_y.cols(0, s_cols-1) = cluster_normal_y_.cols(0, s_cols-1);
        new_normal_z.cols(0, s_cols-1) = cluster_normal_z_.cols(0, s_cols-1);
    }
    // printf("00 cluster num : %zu , s_rows : %zu, s_cols: %zu \n", c_num, s_rows, s_cols);

    cluster_normal_x_ = new_normal_x;
    cluster_normal_y_ = new_normal_y;
    cluster_normal_z_ = new_normal_z;

    vipss_time_stats_.clear();
    for(size_t i = s_cols; i < c_num; ++i)
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
    size_t c_num = cluster_cores_mat_.n_cols;

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
                cluster_normal_x_.col(n_cid) *= -1.0;
                cluster_normal_y_.col(n_cid) *= -1.0;
                cluster_normal_z_.col(n_cid) *= -1.0;
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
        const arma::sp_imat& adj_row = cluster_adjacent_mat_.col(cur_pid);
        // printf("row %d contains no zero num : %d \n", cur_pid, adj_row.n_nonzero);
        const arma::sp_imat::const_iterator start = adj_row.begin();
        const arma::sp_imat::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {
            size_t n_id = iter.row();
            if(visited_vids.find(n_id) != visited_vids.end()) continue;
            if(n_id == cur_pid) continue;
            C_Edege edge(cur_pid, n_id);
            edge.score_ = cluster_scores_mat_(n_id, cur_pid);
            edge_priority_queue.push(edge);
        }
    }
    size_t c_num = cluster_cores_mat_.n_cols;
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
    size_t c_num = cluster_cores_mat_.n_cols;
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

        const arma::sp_imat& adj_row = cluster_MST_mat_.col(cur_cid);
        const arma::sp_imat::const_iterator start = adj_row.begin();
        const arma::sp_imat::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {   
            size_t n_cid = iter.row();
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
                cluster_normal_x_.col(n_cid) *= -1.0;
                cluster_normal_y_.col(n_cid) *= -1.0;
                cluster_normal_z_.col(n_cid) *= -1.0;
            }
            cluster_queued_ids.push(n_cid);
        }
    }
}

// void LocalVipss::FlipClusterNormalsByMinST()
// {
//     size_t pt_num = adjacent_mat_.n_rows;
//     std::set<std::string> visited_edge_ids;
//     std::vector<C_Edege> tree_edges;

//     auto cmp = [](C_Edege& left, C_Edege& right) { return left.score_ > right.score_; };
//     std::priority_queue<C_Edege, std::vector<C_Edege>, decltype(cmp)> edge_priority_queue(cmp);
 
//     C_Edege st_e(0,0);
//     edge_priority_queue.push(st_e);
//     std::set<size_t> visited_vids;
//     while(!edge_priority_queue.empty())
//     {
//         C_Edege cur_e = edge_priority_queue.top();
//         edge_priority_queue.pop();

//         if(visited_vids.find(cur_e.c_b_ ) != visited_vids.end()) continue;
//         visited_vids.insert(cur_e.c_b_);
//         if(cur_e.c_a_ != cur_e.c_b_)
//         {
//             tree_edges.push_back(cur_e);
//         }

//         size_t cur_pid = cur_e.c_b_ ;
//         const arma::sp_irowvec& adj_row = adjacent_mat_.row(cur_pid);
//         // printf("row %d contains no zero num : %d \n", cur_pid, adj_row.n_nonzero);
//         const arma::sp_irowvec::const_iterator start = adj_row.begin();
//         const arma::sp_irowvec::const_iterator end = adj_row.end();
//         for(auto iter = start; iter != end; ++iter)
//         {
//             size_t n_id = iter.internal_col;
//             if(visited_vids.find(n_id) != visited_vids.end()) continue;
//             if(n_id == cur_pid) continue;
//             C_Edege edge(cur_pid, n_id);
//             edge.score_ = PtDistance(points_[cur_pid], points_[n_id]);
//             edge_priority_queue.push(edge);
//         }
//     }

//     std::string out_tree_path = out_dir_ + filename_ + "tree.obj";
//     std::vector<double> vts ;
//     for(auto& p : points_)
//     {
//         vts.push_back(p[0]);
//         vts.push_back(p[1]);
//         vts.push_back(p[2]);
//     }
//     std::vector<unsigned int> edges_out;
//     for(auto& e : tree_edges) 
//     {
//         edges_out.push_back(e.c_a_);
//         edges_out.push_back(e.c_b_);
//     }

//     writeObjFile_line(out_tree_path,vts, edges_out);
// }

void LocalVipss::OuputPtN(const std::string& out_path, bool orient_normal)
{
    out_pts_.clear();
    out_normals_.clear();
    std::vector<double>& pts = out_pts_;
    std::vector<double>& normals = out_normals_;

    size_t npt = points_.size();
    normals.resize(npt*3);

    for(size_t i = 0; i < npt; ++i)
    {
        auto& pt = points_[i];
        pts.push_back(pt[0]);
        pts.push_back(pt[1]);
        pts.push_back(pt[2]);
    }

    for(size_t i = 0; i < cluster_cores_mat_.n_cols; ++i)
    {
        const arma::sp_imat new_row(cluster_cores_mat_.col(i));
        const arma::sp_imat::const_iterator start = new_row.begin();
        const arma::sp_imat::const_iterator end = new_row.end();
        for(auto iter = start; iter != end; ++iter)
        {
            size_t p_id = iter.row();     
            normals[p_id *3]     = cluster_normal_x_(p_id, i);
            normals[p_id *3 + 1] = cluster_normal_y_(p_id, i);
            normals[p_id *3 + 2] = cluster_normal_z_(p_id, i);
        }
    }
    

    // for(size_t i = 0; i < cluster_cores_mat_.n_cols; ++i)
    // {
    //     // printf("cluster id : %d \n", i);
    //     // const arma::sp_icolvec cur_row = cluster_cores_mat_.col(i);
    //     const arma::sp_imat new_row(cluster_cores_mat_.col(i));
    //     const arma::sp_imat::const_iterator start = new_row.begin();
    //     const arma::sp_imat::const_iterator end = new_row.end();
    //     for(auto iter = start; iter != end; ++iter)
    //     {
    //         size_t p_id = iter.row();
            
    //         // printf("p_id : %d \n", p_id);
    //         // printf("non zero count : %d \n", new_row.n_nonzero);
    //         auto pt = points_[p_id];
    //         pts.push_back(pt[0]);
    //         pts.push_back(pt[1]);
    //         pts.push_back(pt[2]);
    //         normals.push_back(cluster_normal_x_(p_id, i));
    //         normals.push_back(cluster_normal_y_(p_id, i));
    //         normals.push_back(cluster_normal_z_(p_id, i));
    //     }
    // }
    // printf("out pt size : %zu \n", pts.size() / 3);

    // if(orient_normal)
    // {
    //     ORIENT::OrientPointNormals(pts, normals);
    // }
    
    // writePLYFile_VN(out_path, pts, normals);
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

void LocalVipss::SaveCluster()
{
    std::vector<std::vector<P3tr>> key_pts_vec;
    for(size_t i = 0; i < cluster_cores_mat_.n_cols; ++i)
    {
        const arma::sp_imat& cur_row(cluster_cores_mat_.col(i));
        if(cur_row.n_nonzero < 3) continue;
        
        const arma::sp_imat::const_iterator start = cur_row.begin();
        const arma::sp_imat::const_iterator end = cur_row.end();
        std::set<P3tr> key_pts;
        for(auto iter = start; iter != end; ++iter)
        {
            key_pts.insert(points_[iter.row()]);
        } 
        
        key_pts_vec.push_back(std::vector<P3tr>(key_pts.begin(), key_pts.end()) );
        // std::vector<P3tr> nei_pts;
        // const arma::sp_irowvec& cur_pt_row = cluster_adjacent_pt_mat_.row(i);
        // const arma::sp_irowvec::const_iterator pt_start = cur_pt_row.begin();
        // const arma::sp_irowvec::const_iterator pt_end = cur_pt_row.end();
        // for(auto iter = pt_start; iter != pt_end; ++iter)
        // {
        //     auto cur_pt = points_[iter.internal_col];
        //     if(key_pts.find(cur_pt) == key_pts.end())
        //     {
        //         nei_pts.push_back(cur_pt);
        //     }
        // }
        // std::vector<P3tr> key_pts_vec(key_pts.begin(), key_pts.end());
        // SaveClusterPts(out_path, key_pts_vec, nei_pts);
    }
    std::sort(key_pts_vec.begin(), key_pts_vec.end(), [](const auto& a, const auto& b)
    {
        return a.size() > b.size();
    });
    std::vector<std::vector<P3tr>> top_k_key_pts_vec;
    for(size_t i = 0; i < top_k_ && i < key_pts_vec.size(); ++i)
    {
        top_k_key_pts_vec.push_back(key_pts_vec[i]);
    }
    std::string out_path = out_dir_ + filename_ + "_cluster_core_pts";
    SaveClusterCorePts(out_path, top_k_key_pts_vec);

}



void LocalVipss::SaveClusterPts(const std::string& path,
                            const std::vector<P3tr>& key_pts, 
                            const std::vector<P3tr>& nei_pts)
{
    std::vector<double> pts;
    std::vector<uint8_t> colors;
    size_t key_num = key_pts.size();
    size_t nei_num = nei_pts.size();
    pts.resize(3*(key_num + nei_num));
    colors.resize(3*(key_num + nei_num));
    for(size_t i = 0; i < key_num + nei_num; ++i)
    {
        if(i < key_num)
        {
            pts[3*i] = key_pts[i][0];
            pts[3*i + 1] = key_pts[i][1];
            pts[3*i + 2] = key_pts[i][2];
            colors[3*i] = 255;
            colors[3*i + 1] = 0;
            colors[3*i + 2] = 0;
        } else {
            size_t id = i - key_num;
            pts[3*id] = nei_pts[id][0];
            pts[3*id + 1] = nei_pts[id][1];
            pts[3*id + 2] = nei_pts[id][2];
            colors[3*id] = 0;
            colors[3*id + 1] = 0;
            colors[3*id + 2] = 255;
        } 
    }
    writePLYFile_CO(path, pts, colors);
}

void LocalVipss::SaveClusterCorePts(const std::string& path,
                            const std::vector<std::vector<P3tr>>& key_pts_vec)
{
    std::vector<double> pts;
    std::vector<uint8_t> colors;
    int k = key_pts_vec.size();
    float step = 255.0 / k * 2.0;
    std::vector<int> red_steps;
    std::vector<int> blue_steps;
    for(size_t i = 0; i < k; ++i)
    {
        uint8_t red = i < k / 2.0 ? uint8_t(std::abs(255 - step * i)) : 0;
        uint8_t green = uint8_t(150 - std::abs(float(128 - step *i)));
        uint8_t blue = i > k / 2.0 ? uint8_t(255 - step * (k-i-1 )) :  0;
        
        const auto& key_pts = key_pts_vec[i];
        for(auto pt : key_pts)
        {
            pts.push_back(pt[0]);
            pts.push_back(pt[1]);
            pts.push_back(pt[2]);
            colors.push_back(red);
            colors.push_back(green);
            colors.push_back(blue);
        }
    }
    writePLYFile_CO(path, pts, colors);  
}




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

    std::string init_ptn_path = out_dir_ + filename_ + "_init_ptn";
    OuputPtN(init_ptn_path);

    InitNormalWithVipss();
    auto t12 = Clock::now();
    double normal_estimate_time = std::chrono::nanoseconds(t12 - t1).count()/1e9;
    printf("finish init cluster normals time : %f ! \n", normal_estimate_time);

    total_time += build_mat_time + normal_estimate_time;
    

    double sum_score_time = 0;

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
    total_time += cluster_scores_time;

    sum_score_time += scores_time;
    sum_score_time += cluster_scores_time;

    size_t iter_num = 1;

    double sum_merge_time = 0;
    while(merged_cluster_size_ > 0)
    {
        iter_num ++;
        if(iter_num > max_iter_) 
        {
            iter_num = max_iter_;
            break;
        }
        auto tt0 = Clock::now();
        // GetClusterCenters();
        MergeClusters();
        auto tt1 = Clock::now();
        double merge_time = std::chrono::nanoseconds(tt1 - tt0).count()/1e9;
        total_time += merge_time;
        sum_merge_time += merge_time;
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

        if(cluster_cores_mat_.n_cols <= 1) break;
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
        sum_score_time += update_score_time + calculate_score_time;
        // std::string init_ptn_path_iter = out_dir_ + filename_ + "_ptn" + std::to_string(iter_num);
        // OuputPtN(init_ptn_path_iter);
        // std::string init_ptn_path_iter2 = out_dir_ + filename_ + "_orient_ptn" + std::to_string(iter_num);
        // OuputPtN(init_ptn_path_iter2, true);
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
    SaveCluster();
    WriteVipssTimeLog();
    printf("total time used : %f \n", total_time);
    printf("total merge time used : %f \n", sum_merge_time);
    printf("total score time used : %f \n", sum_score_time);
    if(use_hrbf_surface_)
    {
        vipss_api_.is_surfacing_ = true;
        vipss_api_.run_vipss(out_pts_, out_normals_);
    }
}


void LocalVipss::InitNormals()
{
    double total_time = 0;
    auto t0 = Clock::now();
    BuildClusterAdjacentMat();
    // printf("finish build cluster adjacent mat ! \n");

    // printf("2 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    BuidClusterCoresPtIds();
    auto t1 = Clock::now();
    double build_mat_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("finish init core pt ids time : %f ! \n", build_mat_time);
    // std::string init_ptn_path = out_dir_ + filename_ + "_init_ptn";
    // OuputPtN(init_ptn_path);

    InitNormalWithVipss();
    auto t12 = Clock::now();
    double normal_estimate_time = std::chrono::nanoseconds(t12 - t1).count()/1e9;
    printf("finish init cluster normals time : %f ! \n", normal_estimate_time);

    total_time += build_mat_time + normal_estimate_time;
    double sum_score_time = 0;

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
    total_time += cluster_scores_time;
    sum_score_time += scores_time;
    sum_score_time += cluster_scores_time;

    size_t iter_num = 1;
    auto ti00 = Clock::now();
    BuildClusterMST();
    auto ti11 = Clock::now();
    double MST_time = std::chrono::nanoseconds(ti11 - ti00).count()/1e9;
    total_time += MST_time;
    printf("normal MST_time time used : %f \n", MST_time); 
    auto ti0 = Clock::now();
    FlipClusterNormalsByMST();
    auto ti1 = Clock::now();
    double flip_time = std::chrono::nanoseconds(ti1 - ti0).count()/1e9;
    total_time += flip_time;
    printf("normal flip time used : %f \n", flip_time);

    std::string init_ptn_path_iter = out_dir_ + filename_ + "_flipped" + std::to_string(iter_num);
    OuputPtN(init_ptn_path_iter);
    // SaveCluster();
    // WriteVipssTimeLog();
    printf("total time used : %f \n", total_time);
    printf("total score time used : %f \n", sum_score_time);
    // if(use_hrbf_surface_)
    // {
    //     vipss_api_.is_surfacing_ = true;
    //     vipss_api_.run_vipss(out_pts_, out_normals_);
    // }
}

void LocalVipss::InitNormalsWithMerge()
{
    // printf("00000 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    double total_time = 0;
    auto t0 = Clock::now();
    BuildClusterAdjacentMat();
    // printf("2 adjacent mat rows : %d, cols : %d \n", adjacent_mat_.n_rows, adjacent_mat_.n_cols);
    BuidClusterCoresPtIds();
    auto t1 = Clock::now();
    double build_mat_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("finish init core pt ids time : %f ! \n", build_mat_time);

    std::string init_ptn_path = out_dir_ + filename_ + "_init_ptn";
    OuputPtN(init_ptn_path);

    InitNormalWithVipss();
    auto t12 = Clock::now();
    double normal_estimate_time = std::chrono::nanoseconds(t12 - t1).count()/1e9;
    printf("finish init cluster normals time : %f ! \n", normal_estimate_time);
    total_time += build_mat_time + normal_estimate_time;
    double sum_score_time = 0;

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
    total_time += cluster_scores_time;
    sum_score_time += scores_time;
    sum_score_time += cluster_scores_time;

    size_t iter_num = 1;

    double sum_merge_time = 0;
    while(merged_cluster_size_ > 0)
    {
        iter_num ++;
        if(iter_num > max_iter_) 
        {
            iter_num = max_iter_;
            break;
        }
        auto tt0 = Clock::now();
        // GetClusterCenters();
        MergeClusters();
        auto tt1 = Clock::now();
        double merge_time = std::chrono::nanoseconds(tt1 - tt0).count()/1e9;
        total_time += merge_time;
        sum_merge_time += merge_time;
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

        if(cluster_cores_mat_.n_cols <= 1) break;
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
        sum_score_time += update_score_time + calculate_score_time;
        // std::string init_ptn_path_iter = out_dir_ + filename_ + "_ptn" + std::to_string(iter_num);
        // OuputPtN(init_ptn_path_iter);
        // std::string init_ptn_path_iter2 = out_dir_ + filename_ + "_orient_ptn" + std::to_string(iter_num);
        // OuputPtN(init_ptn_path_iter2, true);
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
    SaveCluster();
    // WriteVipssTimeLog();
    printf("total time used : %f \n", total_time);
    printf("total merge time used : %f \n", sum_merge_time);
    printf("total score time used : %f \n", sum_score_time);
    if(use_hrbf_surface_)
    {
        vipss_api_.is_surfacing_ = true;
        vipss_api_.run_vipss(out_pts_, out_normals_);
    }
}