#include "local_vipss.hpp"
#include <cmath>
#include <algorithm>
#include "readers.h"  
#include <chrono>
#include <queue>
#include <omp.h>
#include <random>
#include "stats.h"
#include "kernel.h"

typedef std::chrono::high_resolution_clock Clock;

double LocalVipss::search_nn_time_sum_ = 0;
double LocalVipss::pass_time_sum_ = 0;
int LocalVipss::ave_voxel_nn_pt_num_ = 0;

std::vector<tetgenmesh::point> LocalVipss::points_;
// std::vector<std::vector<size_t>> LocalVipss::cluster_all_pt_ids_;


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
        printf("insert id : %ld \n", i);
    }

    auto t1 = Clock::now();
    double insert_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf(" point insert time total: %g \n", insert_time);

    // for(auto &pt : insert_pts_)
    // {
    //     voro_gen_.InsertPt(pt);
    // }
}

void LocalVipss::PtPCA(std::vector<double>& pts)
{
    // std::vect<double> pts;
    size_t ptn = points_.size();
    arma::mat pts_mat(ptn, 3);
    for(size_t i = 0; i < ptn; ++i)
    {
        pts_mat(i, 0) = points_[i][0];
        pts_mat(i, 1) = points_[i][1];
        pts_mat(i, 2) = points_[i][2];
    }
    arma::mat coeff = arma::princomp(pts_mat);
    std::cout << " pca coeffs:  " << coeff.n_rows << " " << coeff.n_cols << std::endl;
}

void LocalVipss::GroupPtsWithVolume()
{

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

void LocalVipss::Init(const std::string & path, const std::string& ext)
{   
    std::vector<double> in_pts;
    std::vector<double> in_normals;
    if(ext == ".ply") 
    {
        readPLYFile(path, in_pts, in_normals);
    } else {
        readXYZ(path, in_pts);
    }
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
    arma::sp_umat cluster_core_ids(cluster_cores_mat_.col(cluster_id));
    arma::sp_umat::const_iterator start = cluster_core_ids.begin();
    arma::sp_umat::const_iterator end = cluster_core_ids.end();
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

void LocalVipss::AddClusterHMatrix(const std::vector<size_t>& p_ids, const arma::mat& J_m, size_t npt)
{
    size_t unit_npt = p_ids.size();
    size_t unit_key_npt = unit_npt;
    std::vector<Triplet> coefficients;
    for(size_t i = 0; i < unit_npt; ++i)
    {
        size_t pi = p_ids[i];
        if(user_lambda_ > 1e-12)
        {
            for(size_t j = 0; j < unit_npt; ++j)
            {
                h_ele_triplets_.push_back(std::move(Triplet(pi, p_ids[j], J_m(i, j))));
            }
            for(size_t j = 0; j < unit_key_npt; ++j)
            {
                h_ele_triplets_.push_back(std::move(Triplet(pi, p_ids[j] + npt, J_m(i, j + unit_npt))));
                h_ele_triplets_.push_back(std::move(Triplet(pi, p_ids[j] + npt*2, J_m(i, j + unit_npt + unit_key_npt)))); 
                h_ele_triplets_.push_back(std::move(Triplet(pi, p_ids[j] + npt*3, J_m(i, j + unit_npt + unit_key_npt* 2))));   
            }
            for(size_t j = 0; j < unit_key_npt; ++j)
            {
                h_ele_triplets_.push_back(std::move(Triplet(p_ids[j] + npt,   pi, J_m(i, j + unit_npt))));
                h_ele_triplets_.push_back(std::move(Triplet(p_ids[j] + npt*2, pi, J_m(i, j + unit_npt + unit_key_npt)))); 
                h_ele_triplets_.push_back(std::move(Triplet(p_ids[j] + npt*3, pi, J_m(i, j + unit_npt + unit_key_npt* 2))));   
            }
        }
        
        if(i < unit_key_npt)
        {
            for(size_t step = 0; step < 3; ++step)
            {
                size_t row_i = npt + npt * step + pi;
                size_t rowv_i = i + unit_npt + unit_key_npt * step;
                for(size_t j = 0; j < unit_key_npt; ++j)
                {
                    h_ele_triplets_.push_back(std::move(Triplet(row_i, p_ids[j] + npt,     J_m( rowv_i, j + unit_npt))));
                    h_ele_triplets_.push_back(std::move(Triplet(row_i, p_ids[j] + npt * 2, J_m( rowv_i, j + unit_npt + unit_key_npt))));
                    h_ele_triplets_.push_back(std::move(Triplet(row_i, p_ids[j] + npt * 3, J_m( rowv_i, j + unit_npt + unit_key_npt * 2))));
                }
            }
        }

        // for(size_t step = 0; step < 4; ++step)
        // {
        //     size_t row = pi + npt * step;
        //     size_t j_row = i + unit_npt* step;
        //     for(size_t j = 0; j < unit_npt; ++j)
        //     {
        //         size_t pj = p_ids[j];
        //         // h_row_map[pj]           += J_m(j_row, j);
        //         // h_row_map[npt + pj]     += J_m(j_row, j + unit_npt);
        //         // h_row_map[npt * 2 + pj] += J_m(j_row, j + unit_npt * 2);
        //         // h_row_map[npt * 3 + pj] += J_m(j_row, j + unit_npt * 3);

                    
        //         for(size_t k = 0; k < 4; ++k)
        //         {
        //             h_ele_triplets_.push_back(std::move(Triplet(row, pj + npt * k, J_m(j_row, j + unit_npt * k))));
        //         }    
        //     }
        // }
    }
}


void LocalVipss::AddClusterHMatrix(const std::vector<size_t>& p_ids, const arma::mat& J_m, size_t npt, std::vector<Triplet>& ele_vect )
{
    size_t unit_npt = p_ids.size();
    size_t unit_key_npt = unit_npt;
    std::vector<Triplet> coefficients;
    for(size_t i = 0; i < unit_npt; ++i)
    {
        size_t pi = p_ids[i];
        if(user_lambda_ > 1e-12)
        {
            for(size_t j = 0; j < unit_npt; ++j)
            {
                ele_vect.push_back(std::move(Triplet(pi, p_ids[j], J_m(i, j))));
            }
            for(size_t j = 0; j < unit_key_npt; ++j)
            {
                ele_vect.push_back(std::move(Triplet(pi, p_ids[j] + npt, J_m(i, j + unit_npt))));
                ele_vect.push_back(std::move(Triplet(pi, p_ids[j] + npt*2, J_m(i, j + unit_npt + unit_key_npt)))); 
                ele_vect.push_back(std::move(Triplet(pi, p_ids[j] + npt*3, J_m(i, j + unit_npt + unit_key_npt* 2))));   
            }
            for(size_t j = 0; j < unit_key_npt; ++j)
            {
                ele_vect.push_back(std::move(Triplet(p_ids[j] + npt,   pi, J_m(i, j + unit_npt))));
                ele_vect.push_back(std::move(Triplet(p_ids[j] + npt*2, pi, J_m(i, j + unit_npt + unit_key_npt)))); 
                ele_vect.push_back(std::move(Triplet(p_ids[j] + npt*3, pi, J_m(i, j + unit_npt + unit_key_npt* 2))));   
            }
        }
        
        if(i < unit_key_npt)
        {
            for(size_t step = 0; step < 3; ++step)
            {
                size_t row_i = npt + npt * step + pi;
                size_t rowv_i = i + unit_npt + unit_key_npt * step;
                for(size_t j = 0; j < unit_key_npt; ++j)
                {
                    ele_vect.push_back(std::move(Triplet(row_i, p_ids[j] + npt,     J_m( rowv_i, j + unit_npt))));
                    ele_vect.push_back(std::move(Triplet(row_i, p_ids[j] + npt * 2, J_m( rowv_i, j + unit_npt + unit_key_npt))));
                    ele_vect.push_back(std::move(Triplet(row_i, p_ids[j] + npt * 3, J_m( rowv_i, j + unit_npt + unit_key_npt * 2))));
                }
            }
        }
    }
}


void LocalVipss::AddClusterHMatrix(const std::vector<size_t>& p_ids, const arma::mat& J_m, size_t npt, 
                                    std::vector<Triplet>::iterator& ele_iter )
{
    size_t unit_npt = p_ids.size();
    size_t unit_key_npt = unit_npt;
    std::vector<Triplet> coefficients;
    for(size_t i = 0; i < unit_npt; ++i)
    {
        size_t pi = p_ids[i];
        if(user_lambda_ > 1e-12)
        {
            for(size_t j = 0; j < unit_npt; ++j)
            {
                *(ele_iter ++) =  Triplet(pi, p_ids[j], J_m(i, j));
            }
            for(size_t j = 0; j < unit_key_npt; ++j)
            {
                *(ele_iter ++) = Triplet(pi, p_ids[j] + npt, J_m(i, j + unit_npt));
                *(ele_iter ++) = Triplet(pi, p_ids[j] + npt*2, J_m(i, j + unit_npt + unit_key_npt)); 
                *(ele_iter ++) = Triplet(pi, p_ids[j] + npt*3, J_m(i, j + unit_npt + unit_key_npt* 2));   
            }
            for(size_t j = 0; j < unit_key_npt; ++j)
            {
                *(ele_iter ++) = Triplet(p_ids[j] + npt,   pi, J_m(i, j + unit_npt));
                *(ele_iter ++) = Triplet(p_ids[j] + npt*2, pi, J_m(i, j + unit_npt + unit_key_npt)); 
                *(ele_iter ++) = Triplet(p_ids[j] + npt*3, pi, J_m(i, j + unit_npt + unit_key_npt* 2));   
            }
        }
        
        if(i < unit_key_npt)
        {
            for(size_t step = 0; step < 3; ++step)
            {
                size_t row_i = npt + npt * step + pi;
                size_t rowv_i = i + unit_npt + unit_key_npt * step;
                for(size_t j = 0; j < unit_key_npt; ++j)
                {
                    *(ele_iter ++) = Triplet(row_i, p_ids[j] + npt,     J_m( rowv_i, j + unit_npt));
                    *(ele_iter ++) = Triplet(row_i, p_ids[j] + npt * 2, J_m( rowv_i, j + unit_npt + unit_key_npt));
                    *(ele_iter ++) = Triplet(row_i, p_ids[j] + npt * 3, J_m( rowv_i, j + unit_npt + unit_key_npt * 2));
                }
            }
        }
    }
}


void LocalVipss::BuildMatrixH()
{
    auto t00 = Clock::now(); 
    size_t npt = this->points_.size();
    size_t cluster_num = cluster_cores_mat_.n_cols;
    final_h_eigen_.resize(4 * npt , 4 * npt);
    double add_ele_to_vector_time = 0;

    // int all_ele_num = arma::dot(VoronoiGen::cluster_size_vec_ , VoronoiGen::cluster_size_vec_) * 16;
    arma::ivec cluster_j_size_vec = VoronoiGen::cluster_size_vec_  % VoronoiGen::cluster_size_vec_ * 16; 
    int all_ele_num = arma::accu(cluster_j_size_vec);
    h_ele_triplets_.resize(all_ele_num);
    printf("average cluster J ele num : %d \n", int(arma::mean(cluster_j_size_vec)));

    arma::ivec acc_j_size_vec = arma::cumsum(cluster_j_size_vec);
    auto iter = h_ele_triplets_.begin();

    //std::vector<std::vector<Triplet>> ele_vector(cluster_num);
    auto t5 = Clock::now();
//#pragma omp parallel for shared(points_, VoronoiGen::cluster_init_pids_) 
    for(int i =0; i < cluster_num; ++i)
    {
        const auto& cluster_pt_ids = VoronoiGen::cluster_init_pids_[i];
        auto Minv =  VIPSSKernel::BuildHrbfMat(points_, cluster_pt_ids);
        size_t unit_npt = cluster_pt_ids.size();
        size_t j_ele_num = unit_npt * unit_npt * 16;

        // cur_eles.reserve(j_ele_num);
        auto cur_iter = iter + acc_j_size_vec[i];
        if(user_lambda_ > 1e-10)
        {
            arma::mat F(4 * unit_npt, 4 * unit_npt);
            arma::mat E(unit_npt, unit_npt);
            E.eye();
            F(0, 0, arma::size(unit_npt, unit_npt)) = E;
            double cur_lambda = user_lambda_;
            arma::mat K = (F + Minv * cur_lambda);
            AddClusterHMatrix(cluster_pt_ids, K, npt, cur_iter);
        } else {
            AddClusterHMatrix(cluster_pt_ids, Minv, npt, cur_iter);
        }
    }
    auto t6 = Clock::now();
    double add_time = std::chrono::nanoseconds(t6 - t5).count()/1e9;
    add_ele_to_vector_time += add_time;
    
    auto t_h1 = Clock::now();
    auto t_h111 = Clock::now();
    double push_to_vec_time = std::chrono::nanoseconds(t_h111 - t_h1).count() / 1e9;
    printf("--- push ele to vector  time : %f \n", push_to_vec_time);

    final_h_eigen_.setFromTriplets(h_ele_triplets_.begin(), h_ele_triplets_.end());
    auto t_h2 = Clock::now();
    double build_h_from_tris_time = std::chrono::nanoseconds(t_h2 - t_h1).count()/1e9;
    
    auto t11 = Clock::now();
    double build_H_time_total = std::chrono::nanoseconds(t11 - t00).count()/1e9;
    printf("--- build vipss j total  time : %f \n", build_j_time_total_);
    printf("--- add  J matrix to triplet vector time : %f \n", add_ele_to_vector_time);
    printf("--- build final_h  from triplets vector time : %f \n", build_h_from_tris_time);
    printf("--- build_H_time_total sum  time : %f \n", build_H_time_total);

    G_VP_stats.build_H_total_time_ += build_H_time_total;
    G_VP_stats.cal_cluster_J_total_time_ += build_j_time_total_;
    G_VP_stats.add_J_ele_to_triplet_vector_time_ += add_ele_to_vector_time;
    G_VP_stats.build_eigen_final_h_from_tris_time_ += build_h_from_tris_time;

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
std::vector<std::shared_ptr<RBF_Core>> LocalVipss::node_rbf_vec_;
VoronoiGen LocalVipss::voro_gen_;

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
        std::vector<double> cluster_pt_vec = VoronoiGen::cluster_init_pts_[i];
        const std::vector<size_t>& cluster_pt_ids = VoronoiGen::cluster_init_pids_[i];
        std::vector<double> cluster_nl_vec;
        if(use_partial_vipss) 
        {
            cluster_nl_vec = {normals_[3*i], normals_[3*i + 1], normals_[3*i + 2]};
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


void LocalVipss::BuildHRBFPerCluster()
{
    auto t00 = Clock::now(); 
    size_t npt = this->points_.size();
    size_t cluster_num = points_.size();
    double time_sum = 0;
    node_rbf_vec_.resize(cluster_num);
    cluster_volume_vals_.resize(cluster_num);
    bool use_partial_vipss = false;
    vipss_api_.user_lambda_ = user_lambda_;

    // printf("cluster num : %lu \n", cluster_num);
    size_t valid_cluster_num = arma::accu(cluster_valid_sign_vec_);
    // valid_cluster_dist_map_.resize(valid_cluster_num);
    int valid_pt_count =0;
    int valid_core_pt_count =0;
    for(size_t i =0; i < cluster_num; ++i)
    {

        if(!cluster_valid_sign_vec_[i]) continue;
        const auto& c_pids = cluster_pt_ids_[i];
        // printf("---- cur id : %ld group pt number : %ld \n", i, c_pids.size());
        std::vector<double> cluster_pt_vec;
        std::vector<double> cluster_nl_vec;
        std::vector<double> cluster_sv_vec;
        for(const auto pid : c_pids)
        {
            cluster_pt_vec.push_back(points_[pid][0]);
            cluster_pt_vec.push_back(points_[pid][1]);
            cluster_pt_vec.push_back(points_[pid][2]);
            cluster_nl_vec.push_back(normals_[3*pid]);
            cluster_nl_vec.push_back(normals_[3*pid + 1]);
            cluster_nl_vec.push_back(normals_[3*pid + 2]);
            cluster_sv_vec.push_back(s_vals_[pid]);
            // cluster_id_map_[pid] = i;
        }
        const auto& core_ids = cluster_core_pt_ids_vec_[i];
        valid_pt_count += c_pids.size();
        valid_core_pt_count += core_ids.size();
        for(const auto id : core_ids)
        {
            cluster_id_map_[id] = i;
        }
        node_rbf_vec_[i] = std::make_shared<RBF_Core>();
        vipss_api_.build_cluster_hrbf(cluster_pt_vec, cluster_nl_vec, cluster_sv_vec, node_rbf_vec_[i]);
    }

    auto t11 = Clock::now();
    double build_HRBF_time_total = std::chrono::nanoseconds(t11 - t00).count()/1e9;
    printf("--- build cluster HRBF total  time : %f \n", build_HRBF_time_total);
    printf("valid core pt count : %d \n", valid_pt_count);
    printf("average cluster num : %d \n", valid_pt_count/int(valid_cluster_num));
    printf("average cluster core pt num : %d \n", valid_core_pt_count/int(valid_cluster_num));
}

double LocalVipss::NodeDistanceFunction(const tetgenmesh::point nn_pt, const tetgenmesh::point cur_pt) const
{
    if(voro_gen_.point_id_map_.find(nn_pt) != voro_gen_.point_id_map_.end())
    {
        size_t pid = voro_gen_.point_id_map_[nn_pt];
        return node_rbf_vec_[pid]->Dist_Function(cur_pt);
    } 
    return 0;
}

double LocalVipss::NatureNeighborDistanceFunctionOMP(const tetgenmesh::point cur_pt) const
{
    std::vector<tetgenmesh::point> nei_pts;
    auto tn0 = Clock::now();
    voro_gen_.GetVoronoiNeiPts(cur_pt, nei_pts);
    auto tn1 = Clock::now();
    double search_nn_time = std::chrono::nanoseconds(tn1 - tn0).count()/1e9;
    search_nn_time_sum_ += search_nn_time;
    // printf("nn num %ld \n", nei_pts.size());
    size_t nn_num = nei_pts.size();
    int i;  
    auto t0 = Clock::now();
    arma::vec nn_dist_vec_(nn_num);
    arma::vec nn_volume_vec_(nn_num);
    ave_voxel_nn_pt_num_ += nn_num;
    const std::vector<double*>& all_pts = points_;

#pragma omp parallel for shared(node_rbf_vec_, voro_gen_, VoronoiGen::point_id_map_, nei_pts, cur_pt, nn_dist_vec_, nn_volume_vec_) private(i)
    for( i = 0; i < nn_num; ++i)
    {
        auto nn_pt = nei_pts[i];
        if(VoronoiGen::point_id_map_.find(nn_pt) != VoronoiGen::point_id_map_.end())
        {
            size_t pid = VoronoiGen::point_id_map_[nn_pt];
            const arma::vec& a = node_rbf_vec_[pid]->a;
            const arma::vec& b = node_rbf_vec_[pid]->b;
            const std::vector<size_t>& cluster_pids = VoronoiGen::cluster_init_pids_[pid];
            nn_dist_vec_[i] = HRBF_Dist_Alone(cur_pt,  a, b, cluster_pids, all_pts);
            // nn_dist_vec_[i] = 1;
            // nn_dist_vec_[i] = node_rbf_vec_[pid]->Dist_Function(cur_pt);
            int thread_id = omp_get_thread_num();
            // nn_volume_vec_[i] = 1.0;
            nn_volume_vec_[i] = voro_gen_.CalTruncatedCellVolumePassOMP(cur_pt, nn_pt, thread_id); 
        }
    }
    auto t1 = Clock::now();
    double pass_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    pass_time_sum_ += pass_time;

    double volume_sum = arma::accu(nn_volume_vec_);
    if(volume_sum > 1e-20)
    {
        return arma::dot(nn_dist_vec_, nn_volume_vec_) / volume_sum;
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
    // cluster_all_pt_ids_.resize(cluster_num);
    node_rbf_vec_.resize(cluster_num);
    int ele_num = arma::accu(VoronoiGen::cluster_size_vec_);
    std::vector<Triplet> cluster_normals_xele(ele_num);
    std::vector<Triplet> cluster_normals_yele(ele_num);
    std::vector<Triplet> cluster_normals_zele(ele_num);

    // std::vector<Triplet> cluster_normals_xele;
    // std::vector<Triplet> cluster_normals_yele;
    // std::vector<Triplet> cluster_normals_zele;

    cluster_normal_x_.resize(npt, npt);
    cluster_normal_y_.resize(npt, npt);
    cluster_normal_z_.resize(npt, npt);
    
#pragma omp parallel for shared(VoronoiGen::cluster_init_pts_, VoronoiGen::cluster_init_pids_, VoronoiGen::cluster_accum_size_vec_)
    for(int i =0; i < int(cluster_num); ++i)
    {
        auto vts = VoronoiGen::cluster_init_pts_[i];
        const auto& p_ids = VoronoiGen::cluster_init_pids_[i];
        // cluster_all_pt_ids_[i] = p_ids;
        // auto t1 = Clock::now();
        size_t unit_npt = p_ids.size(); 
        double cur_lambda = user_lambda_ / double(unit_npt);
        // if(!node_rbf_vec_[i])  
        
        // node_rbf_vec_[i] = std::make_shared<RBF_Core>();
        std::shared_ptr rbf_temp_ptr = std::make_shared<RBF_Core>();
        InitNormalPartialVipss(vts, 1,rbf_temp_ptr, cur_lambda);
        // auto t2 = Clock::now();
        // double vipss_time = std::chrono::nanoseconds(t2 - t1).count()/1e9;
   
        auto iterx = cluster_normals_xele.begin() + VoronoiGen::cluster_accum_size_vec_[i];
        auto itery = cluster_normals_yele.begin() + VoronoiGen::cluster_accum_size_vec_[i];
        auto iterz = cluster_normals_zele.begin() + VoronoiGen::cluster_accum_size_vec_[i];
        
        for(size_t p_id = 0; p_id < p_ids.size(); ++p_id)
        {
            size_t v_id = p_ids[p_id];
            // cluster_normals_xele.push_back(Triplet(v_id, i, node_rbf_vec_[i]->out_normals_[3* p_id]));
            // cluster_normals_yele.push_back(Triplet(v_id, i, node_rbf_vec_[i]->out_normals_[3* p_id + 1]));
            // cluster_normals_zele.push_back(Triplet(v_id, i, node_rbf_vec_[i]->out_normals_[3* p_id + 2]));

            *(iterx + p_id) = Triplet(v_id, i, rbf_temp_ptr->out_normals_[3* p_id]);
            *(itery + p_id) = Triplet(v_id, i, rbf_temp_ptr->out_normals_[3* p_id + 1]);
            *(iterz + p_id) = Triplet(v_id, i, rbf_temp_ptr->out_normals_[3* p_id + 2]);
        }
        // printf("Init vipss cluster id  : %d, cluster pt num : %d \n", i, int(p_ids.size()));
    }

    cluster_normal_x_.setFromTriplets(cluster_normals_xele.begin(), cluster_normals_xele.end());
    cluster_normal_y_.setFromTriplets(cluster_normals_yele.begin(), cluster_normals_yele.end());
    cluster_normal_z_.setFromTriplets(cluster_normals_zele.begin(), cluster_normals_zele.end());

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
    // double dist = local_vipss_ptr->NatureNeighborDistanceFunction(&(new_pt[0]));
    double dist = local_vipss_ptr->NatureNeighborDistanceFunctionOMP(&(new_pt[0]));
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
        if(cluster_normal_x_.coeff(id, c_b) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    for(auto id : core_ids_b)
    {
        if(cluster_normal_x_.coeff(id, c_a) != 0)
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
        normal_ma(i, 0) =  cluster_normal_x_.coeff(id, c_a);
        normal_ma(i, 1) =  cluster_normal_y_.coeff(id, c_a);
        normal_ma(i, 2) =  cluster_normal_z_.coeff(id, c_a);

        normal_mb(i, 0) =  cluster_normal_x_.coeff(id, c_b);
        normal_mb(i, 1) =  cluster_normal_y_.coeff(id, c_b);
        normal_mb(i, 2) =  cluster_normal_z_.coeff(id, c_b);
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
        if(cluster_normal_x_.coeff(id, c_b) != 0)
        {
            valid_ids.push_back(id); 
        }
    }
    for(auto id : core_ids_b)
    {
        if(cluster_normal_x_.coeff(id, c_a) != 0)
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
        normal_ma(i, 0) =  cluster_normal_x_.coeff(id, c_a);
        normal_ma(i, 1) =  cluster_normal_y_.coeff(id, c_a);
        normal_ma(i, 2) =  cluster_normal_z_.coeff(id, c_a);

        normal_mb(i, 0) =  cluster_normal_x_.coeff(id, c_b);
        normal_mb(i, 1) =  cluster_normal_y_.coeff(id, c_b);
        normal_mb(i, 2) =  cluster_normal_z_.coeff(id, c_b);
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
    // cluster_adjacent_pt_mat_ = adjacent_mat_ * cluster_cores_mat_ + cluster_cores_mat_;   
}

void LocalVipss::CalculateClusterNeiScores(bool is_init)
{
    size_t c_num = points_.size();
    // cluster_scores_vec_.resize(c_num);
    int score_ele_num = arma::accu(VoronoiGen::cluster_size_vec_) - c_num;
    std::vector<Triplet> score_eles(score_ele_num);

#pragma omp parallel for
    for(int i = 0; i < c_num; ++i)
    {
        // InitSingleClusterNeiScores(i);
        const auto& nei_pt_ids = VoronoiGen::cluster_init_pids_[i];
        arma::mat cur_i_mat(2, 3);
        arma::mat cur_n_mat(2, 3);
        cur_i_mat(0, 0) = cluster_normal_x_.coeff(i, i);
        cur_i_mat(0, 1) = cluster_normal_y_.coeff(i, i);
        cur_i_mat(0, 2) = cluster_normal_z_.coeff(i, i);
        auto ele_iter =  score_eles.begin() + VoronoiGen::cluster_accum_size_vec_[i] - i; 
        // for(auto n_pt : nei_pts)
        for(auto n_pid : nei_pt_ids)
        {
            if(n_pid == i) continue;
            cur_i_mat(1, 0) = cluster_normal_x_.coeff(n_pid, i);
            cur_i_mat(1, 1) = cluster_normal_y_.coeff(n_pid, i);
            cur_i_mat(1, 2) = cluster_normal_z_.coeff(n_pid, i);

            cur_n_mat(0, 0) = cluster_normal_x_.coeff(i, n_pid);
            cur_n_mat(0, 1) = cluster_normal_y_.coeff(i, n_pid);
            cur_n_mat(0, 2) = cluster_normal_z_.coeff(i, n_pid);

            cur_n_mat(1, 0) = cluster_normal_x_.coeff(n_pid, n_pid);
            cur_n_mat(1, 1) = cluster_normal_y_.coeff(n_pid, n_pid);
            cur_n_mat(1, 2) = cluster_normal_z_.coeff(n_pid, n_pid);

            arma::mat dot_res = cur_i_mat % cur_n_mat;
            arma::vec dot_sum = arma::sum(dot_res, 1); 
            double s1 = dot_sum.min();
            double s2 = -dot_sum.max();
            if(s2 > s1) 
            {
                s1 = s2;
            }
            double score = std::min(1.0, std::max(-1.0, s1));
            score = acos(score) * Anlge_PI_Rate ;
            *(ele_iter ++) = Triplet(n_pid, i, score);
            // printf("cluster id : %d \n", i);   
        // #pragma omp critical
            // cluster_scores_mat_(n_pid, i) = score;
        // #pragma omp critical
            // cluster_scores_mat_(i, n_pid) = score;
        }
        
    }
    cluster_scores_mat_.setFromTriplets(score_eles.begin(), score_eles.end());
}

void LocalVipss::GroupClustersWithDegree()
{
    auto t0 = Clock::now();
    size_t cluster_size = cluster_adjacent_mat_.n_cols;

    cluster_valid_sign_vec_.ones(cluster_size);
    cluster_core_pt_nums_.ones(cluster_size);
    
    // OptimizeAdjacentMat();
    // ConvertToEigenSparseMat();
    // std::cout << " start to init adjacent data " << std::endl;
    InitAdjacentData();
    for(size_t i = 0; i < max_group_iter_; ++i)
    {
        // std::cout << " iter :: " << i << std::endl;
        // CalClusterDegrees();
        // MergeHRBFClustersWithDegree();
        // MergeHRBFClustersWithEigen();
        MergeHRBFClustersWithMap();
    }
    auto t1 = Clock::now();
    double merge_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("------finish cluster merge time : %f ! \n", merge_time);
}

void LocalVipss::InitAdjacentData()
{
    int cols = cluster_adjacent_mat_.n_cols;
    cluster_adjacent_ids_.resize(cols);
    cluster_id_map_.resize(cols);
    cluster_core_pt_ids_vec_.resize(cols);
    cluster_pt_ids_.resize(cols);
    std::cout << " start to get col vec by id "<< std::endl;
    for(size_t i = 0; i < cols; ++i)
    {
        // const arma::sp_ucolvec& col_vec = cluster_adjacent_mat_.col(i);
        // arma::sp_umat col_vec(cluster_adjacent_mat_opt_.col(i));
        arma::sp_umat col_vec(cluster_adjacent_mat_.col(i));
        // std::cout << " get col vec by id "<< std::endl;
        arma::sp_umat::const_iterator start  = col_vec.begin();
        arma::sp_umat::const_iterator end    = col_vec.end();

        // std::cout << " get col vec by id 00"<< std::endl;
        for(auto iter = start; iter != end; ++iter)
        {
            cluster_adjacent_ids_[i].insert(size_t(iter.row()));
            cluster_pt_ids_[i].insert(size_t(iter.row()));
        }
        // std::cout << " get col vec by id 1"<< std::endl;
        cluster_pt_ids_[i].insert(i);
        cluster_id_map_[i] = i;
        cluster_core_pt_ids_vec_[i].push_back(i);
    }
    // cluster_core_pt_nums_.ones(cols); 
    
}


void LocalVipss::MergeHRBFClustersWithMap()
{
    // auto t0 = Clock::now();
    size_t cluster_num = cluster_adjacent_mat_.n_cols;
    int max_val = std::numeric_limits<int>::max();
    for(size_t i = 0; i < cluster_num; ++i)
    {
        if(cluster_valid_sign_vec_[i])
        {
            std::unordered_set<size_t> update_ids;
            for(const auto id :  cluster_adjacent_ids_[i])
            {
                update_ids.insert(cluster_id_map_[id]);
            }
            cluster_adjacent_ids_[i] = update_ids;
        } 
    }

    // std::cout << " update cluster_adjacent_ids_  " << std::endl;

    std::vector<arma::uword> merged_cluster_ids;
    arma::uvec indices = arma::sort_index(cluster_core_pt_nums_, "ascend");
    std::vector<bool> visited_c_vec(indices.n_elem, false);

    for(auto& id : indices)
    {
        if(visited_c_vec[id]) continue;
        visited_c_vec[id] = true;
        if(!cluster_valid_sign_vec_[id]) continue;
        
        {
            int cand_id = -1;
            int max_share_adj_pt_num = 0;
            int cur_core_pt_num = cluster_core_pt_nums_[id];
            auto& cur_pids = cluster_pt_ids_[id];
            for(const auto nei_id :  cluster_adjacent_ids_[id])
            {
                if(visited_c_vec[nei_id]) continue;
                // if(nid == id) continue;
                if(!cluster_valid_sign_vec_[nei_id]) continue;
                if(cluster_core_pt_nums_[id] + cluster_core_pt_nums_[nei_id] > max_group_pt_num_) continue;
                const auto n_pids = cluster_pt_ids_[nei_id];
                int share_pt_num = 0;
                for(const auto npid : n_pids)
                {
                    if(cur_pids.find(npid) != cur_pids.end()) share_pt_num ++;
                }
                if(share_pt_num > max_share_adj_pt_num )
                {
                    cand_id = nei_id;
                    max_share_adj_pt_num = share_pt_num;
                }
            }
            //  std::cout << " finish get candid ids  " << std::endl;
            if(cand_id != -1)
            {   
                visited_c_vec[cand_id] = true; 
                // merged_cluster_ids.push_back(id);
                // merged_cluster_ids.push_back(cand_id);
                size_t ca = id; 
                size_t cb = cand_id;
                cluster_valid_sign_vec_[cb] = 0;
                cluster_id_map_[cb] = ca;
                cluster_core_pt_nums_[ca] += cluster_core_pt_nums_[cb]; 
                cluster_core_pt_nums_[cb] = std::numeric_limits<size_t>::max();
                const auto n_pids = cluster_pt_ids_[cb];
                for(const auto pid : n_pids)
                {
                    cur_pids.insert(pid);
                }
                for(const auto c_id : cluster_adjacent_ids_[cb])
                {
                    cluster_adjacent_ids_[ca].insert(c_id);
                }
                for(const auto core_id : cluster_core_pt_ids_vec_[cand_id])
                {
                    cluster_core_pt_ids_vec_[ca].push_back(core_id);
                }
            }
        }
    }
}
        
void LocalVipss::SaveGroupPtsWithColor(const std::string& path)
{
    size_t group_size = cluster_cores_mat_.n_cols;
    std::ofstream pt_file;
    pt_file.open(path);
    size_t pt_count = 0;
    std::cout << " -----save cluster group size :  " << arma::accu(cluster_valid_sign_vec_) << std::endl;
    for(size_t i = 0; i < group_size; ++i)
    {
        if(!cluster_valid_sign_vec_[i]) continue;
        int base = 1000000;
        double r = ((double) (rand() % base) / double(base)); r = std::min(sqrt(r) , 1.0);
        double g = ((double) (rand() % base) / double(base)); g = std::max(sqrt(g) , 0.0);
        double b = ((double) (rand() % base) / double(base)); b = std::max(sqrt(b) , 0.0);
        // const arma::sp_umat cluster_col(cluster_cores_mat_.col(i));
        // arma::sp_umat::const_iterator start = cluster_col.begin();
        // arma::sp_umat::const_iterator end = cluster_col.end();
        const auto& cur_pids = cluster_core_pt_ids_vec_[i];
        // printf("cur cluster core pt num : %ld \n", cur_pids.size());
        for(auto pid : cur_pids)
        {
            auto cur_pt = points_[pid];
            pt_file << "v " << cur_pt[0] << " " << cur_pt[1] << " " << cur_pt[2] ;
            pt_file << " " << r << " " << g << " " << b << std::endl;
        }
    }
    pt_file.close();
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
        const arma::sp_umat& adj_row = cluster_adjacent_mat_.col(cur_pid);
        // printf("row %d contains no zero num : %d \n", cur_pid, adj_row.n_nonzero);
        const arma::sp_umat::const_iterator start = adj_row.begin();
        const arma::sp_umat::const_iterator end = adj_row.end();
        for(auto iter = start; iter != end; ++iter)
        {
            size_t n_id = iter.row();
            if(visited_vids.find(n_id) != visited_vids.end()) continue;
            if(n_id == cur_pid) continue;
            C_Edege edge(cur_pid, n_id);
            edge.score_ = cluster_scores_mat_.coeff(n_id, cur_pid);
            edge_priority_queue.push(edge);
        }
    }
    size_t c_num = cluster_cores_mat_.n_cols;
    cluster_MST_mat_.resize(c_num, c_num);
    std::vector<TripletInt> edge_eles(tree_edges.size() *2);
    auto e_iter = edge_eles.begin();
    for(const auto& edge: tree_edges)
    {
        size_t i = edge.c_a_;
        size_t j = edge.c_b_;
        *(e_iter ++) = TripletInt(i,j,1);
        *(e_iter ++) = TripletInt(j,i,1);
        // cluster_MST_mat_(i,j) = 1;
        // cluster_MST_mat_(j,i) = 1;
    }
    cluster_MST_mat_.setFromTriplets(edge_eles.begin(), edge_eles.end());
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
        for(SpiMat::InnerIterator iter(cluster_MST_mat_, cur_cid); iter ; ++iter)
        {   
            size_t n_cid = iter.row();
            if(flipped_cluster_ids.find(n_cid) != flipped_cluster_ids.end()) continue;
            flipped_cluster_ids.insert(n_cid);
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
        const arma::sp_umat new_row(cluster_cores_mat_.col(i));
        const arma::sp_umat::const_iterator start = new_row.begin();
        const arma::sp_umat::const_iterator end = new_row.end();
        for(auto iter = start; iter != end; ++iter)
        {
            size_t p_id = iter.row();     
            normals[p_id *3]     = cluster_normal_x_.coeff(p_id, i);
            normals[p_id *3 + 1] = cluster_normal_y_.coeff(p_id, i);
            normals[p_id *3 + 2] = cluster_normal_z_.coeff(p_id, i);
        }
    }
    
    // for(size_t i = 0; i < cluster_cores_mat_.n_cols; ++i)
    // {
    //     // printf("cluster id : %d \n", i);
    //     // const arma::sp_ucolvec cur_row = cluster_cores_mat_.col(i);
    //     const arma::sp_umat new_row(cluster_cores_mat_.col(i));
    //     const arma::sp_umat::const_iterator start = new_row.begin();
    //     const arma::sp_umat::const_iterator end = new_row.end();
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
        const arma::sp_umat& cur_row(cluster_cores_mat_.col(i));
        if(cur_row.n_nonzero < 3) continue;
        
        const arma::sp_umat::const_iterator start = cur_row.begin();
        const arma::sp_umat::const_iterator end = cur_row.end();
        std::set<P3tr> key_pts;
        for(auto iter = start; iter != end; ++iter)
        {
            key_pts.insert(points_[iter.row()]);
        } 
        key_pts_vec.push_back(std::vector<P3tr>(key_pts.begin(), key_pts.end()) );
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
    // SaveClusterCorePts(out_path, top_k_key_pts_vec);

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

// void LocalVipss::SaveClusterCorePts(const std::string& path,
//                             const std::vector<std::vector<P3tr>>& key_pts_vec)
// {
//     std::vector<double> pts;
//     std::vector<uint8_t> colors;
//     int k = key_pts_vec.size();
//     float step = 255.0 / k * 2.0;
//     std::vector<int> red_steps;
//     std::vector<int> blue_steps;
//     for(size_t i = 0; i < k; ++i)
//     {
//         uint8_t red = i < k / 2.0 ? uint8_t(std::abs(255 - step * i)) : 0;
//         uint8_t green = uint8_t(150 - std::abs(float(128 - step *i)));
//         uint8_t blue = i > k / 2.0 ? uint8_t(255 - step * (k-i-1 )) :  0;
        
//         const auto& key_pts = key_pts_vec[i];
//         for(auto pt : key_pts)
//         {
//             pts.push_back(pt[0]);
//             pts.push_back(pt[1]);
//             pts.push_back(pt[2]);
//             colors.push_back(red);
//             colors.push_back(green);
//             colors.push_back(blue);
//         }
//     }
//     writePLYFile_CO(path, pts, colors);  
// }


void LocalVipss::OptimizeAdjacentMat()
{   
    size_t pt_num = cluster_adjacent_mat_.n_cols;
    cluster_adjacent_mat_opt_.resize(pt_num, pt_num);
    std::vector<std::pair<size_t, size_t> > edge_ids;
    for(size_t i = 0; i < pt_num; ++i)
    {
        const arma::sp_umat& adj_col = cluster_adjacent_mat_.col(i);
        const arma::sp_umat::const_iterator start = adj_col.begin();
        const arma::sp_umat::const_iterator end = adj_col.end();
        arma::vec3 p_normal = {out_normals_[3 *i],out_normals_[3 *i + 1],out_normals_[3 *i + 2]};
        std::vector< std::pair<size_t, size_t> > cur_edges;
        std::vector<double> dist_vec;
        arma::vec3 pt = {out_pts_[3*i], out_pts_[3*i + 1], out_pts_[3*i +2]};
        for(auto iter = start; iter != end; ++iter)
        {
            size_t n_id = iter.row();
            arma::vec3 n_normal = {out_normals_[3 *n_id],out_normals_[3 *n_id + 1],out_normals_[3 *n_id + 2]};
            double pro = arma::dot(p_normal, n_normal);
            // if(pro > -0.5)
            {
                arma::vec3 n_pt = {out_pts_[3*n_id], out_pts_[3*n_id + 1], out_pts_[3*n_id +2]};
                
                cur_edges.push_back(std::pair(i, n_id)); 
                double dist = arma::norm(pt- n_pt);
                dist_vec.push_back(dist);
            }
        }
        std::vector<int> v(dist_vec.size());
        for(int j =0; j < int(v.size()); ++j)
        {
            v[j] = j;
        }
        std::sort(v.begin(), v.end(), [&dist_vec](int a, int b){return dist_vec[a] <  dist_vec[b];});
        // int n_num = std::min(int(8), int(dist_vec.size()));
        for(size_t j = 0; j < 8 && j < dist_vec.size() ; ++j)
        {
            const auto& cur_e = cur_edges[v[j]];
            edge_ids.push_back(cur_e);
            cluster_adjacent_mat_opt_(cur_e.second, cur_e.first) = 1;
            cluster_adjacent_mat_opt_(cur_e.first, cur_e.second) = 1;
        }
    }
    // std::cout << "------- cluster_adjacent_mat_opt_  no zero count : " << cluster_adjacent_mat_opt_.n_nonzero << std::endl;
    // std::cout << "------- cluster_adjacent_mat_  no zero count : " << cluster_adjacent_mat_.n_nonzero << std::endl;

    // std::ofstream g_file;
    // std::string graph_path = out_dir_ + "/" + "opt_graph.obj";
    // std::cout << "~~~~~~~ out opt graph path : " << graph_path << std::endl;
    // g_file.open(graph_path);
    // for(size_t i  =0; i < pt_num; ++i)
    // {
    //     g_file << "v " << out_pts_[3*i] << " " << out_pts_[3*i + 1] <<" "<< out_pts_[3*i+ 2] << std::endl;
    //     g_file << "vn "  << out_normals_[3*i] << " " << out_normals_[3*i + 1] <<" "<< out_normals_[3*i+ 2] << std::endl;
    // }
    // for(auto e_id : edge_ids)
    // {
    //     g_file << "l " << e_id.first + 1 << " " << e_id.second + 1 << std::endl; 
    // }
    // g_file.close();
}


void LocalVipss::InitNormals()
{
    double total_time = 0;
    auto t0 = Clock::now();
    BuildClusterAdjacentMat();
    BuidClusterCoresPtIds();
    auto t1 = Clock::now();
    double build_mat_time = std::chrono::nanoseconds(t1 - t0).count()/1e9;
    printf("finish init adj mat and core pt ids time : %f ! \n", build_mat_time);

    InitNormalWithVipss();
    auto t12 = Clock::now();
    double normal_estimate_time = std::chrono::nanoseconds(t12 - t1).count()/1e9;
    printf("finish init cluster normals time : %f ! \n", normal_estimate_time);

    total_time += build_mat_time + normal_estimate_time;
    G_VP_stats.init_cluster_normal_time_ += (build_mat_time + normal_estimate_time);
    double sum_score_time = 0;

    auto t2 = Clock::now();
    CalculateClusterNeiScores(true);
    auto t3 = Clock::now();
    double scores_time = std::chrono::nanoseconds(t3 - t2).count()/1e9;
    printf("finish init cluster neigh scores time : %f ! \n", scores_time);
    // CalculateClusterScores();
    auto t34 = Clock::now();
    double cluster_scores_time = std::chrono::nanoseconds(t34 - t3).count()/1e9;
    printf("finish calculate cluster scores time : %f ! \n", cluster_scores_time);
    total_time += scores_time;
    total_time += cluster_scores_time;
    sum_score_time += scores_time;
    sum_score_time += cluster_scores_time;
    G_VP_stats.cal_cluster_neigbor_scores_time_ += (scores_time + cluster_scores_time);

    size_t iter_num = 1;
    auto ti00 = Clock::now();
    BuildClusterMST();
    auto ti11 = Clock::now();
    double MST_time = std::chrono::nanoseconds(ti11 - ti00).count()/1e9;
    total_time += MST_time;
    G_VP_stats.build_normal_MST_time_ += MST_time;
    printf("normal MST_time time used : %f \n", MST_time); 
    auto ti0 = Clock::now();
    FlipClusterNormalsByMST();
    auto ti1 = Clock::now();
    double flip_time = std::chrono::nanoseconds(ti1 - ti0).count()/1e9;
    G_VP_stats.normal_flip_time_ += flip_time;
    total_time += flip_time;
    printf("normal flip time used : %f \n", flip_time);

    std::string init_ptn_path_iter = out_dir_ + filename_ + "_flipped" + std::to_string(iter_num);
    OuputPtN(init_ptn_path_iter);
    auto finat_t = Clock::now();
    double init_total_time = std::chrono::nanoseconds(finat_t - t0).count()/1e9;

    printf("total time used : %f \n", init_total_time);
    G_VP_stats.init_normal_total_time_ += init_total_time;
}


void LocalVipss::ClearPartialMemory()
{
    adjacent_mat_.clear();
    cluster_adjacent_mat_.clear();
    cluster_MST_mat_.resize(0,0);
    cluster_scores_mat_.resize(1,1);
    cluster_normal_x_.resize(0,0);
    cluster_normal_y_.resize(0,0);
    cluster_normal_z_.resize(0,0);
}