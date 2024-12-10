#include <chrono>
#include "vipss_unit.hpp"
#include "stats.h"
#include "adgrid.h"
#include "SimpleOctree.h"
#include "test_timing.h"

typedef std::chrono::high_resolution_clock Clock;

VP_STATS G_VP_stats;
int VIPSSUnit::opt_func_count_g = 0;
double VIPSSUnit::opt_func_time_g = 0;
Eigen::VectorXd arma_x_opt_g;
Eigen::VectorXd res_vec_g;
Eigen::VectorXd scsc_vec;

void VIPSSUnit::InitPtNormalWithLocalVipss()
{
    std::string& data_dir =  data_dir_;
    local_vipss_.filename_ = file_name_;
    // l_vp.filename_ = "planck";
    local_vipss_.out_dir_ = out_dir_;
    // std::cout << " out dir --------------- " << local_vipss_.out_dir_ << std::endl;
    std::string path = data_dir + local_vipss_.filename_ + "/" + local_vipss_.filename_ + ".ply";
    // local_vipss_.angle_threshold_ = 30;
    // local_vipss_.user_lambda_ = user_lambda_;
    local_vipss_.user_lambda_ = user_lambda_;
    // local_vipss_.max_iter_ = 30;
    local_vipss_.use_hrbf_surface_ = false;
    local_vipss_.angle_threshold_ = merge_angle_;
    auto t1 = Clock::now();
    // local_vipss_.volume_dim_ = 100;
    local_vipss_.Init(input_data_path_, input_data_ext_);
    local_vipss_.InitNormals();
    auto t2 = Clock::now();
    G_VP_stats.init_normal_total_time_ = std::chrono::nanoseconds(t2 - t1).count() / 1e9;
    G_VP_stats.tetgen_triangulation_time_ = local_vipss_.tet_gen_triangulation_time_;
    npt_ = local_vipss_.points_.size();
    initnormals_ = local_vipss_.out_normals_;
    // printf("adaptive grid generation time : %f ! \n", adgrid_gen_time);
    auto t3 = Clock::now();
    local_vipss_.ClearPartialMemory();
    local_vipss_.BuildMatrixH();
    if(user_lambda_ < 1e-10)
    {
         local_vipss_.final_h_eigen_ = local_vipss_.final_h_eigen_.block(npt_, npt_, 3 * npt_, 3 * npt_);
    } 
    auto t4 = Clock::now();
    G_VP_stats.build_H_total_time_ = std::chrono::nanoseconds(t4 - t3).count() / 1e9;
    // double construct_Hmat_time = std::chrono::nanoseconds(t4 - t3).count() / 1e9;
    
    
    // printf("unit vipss J mat init time : %f \n", local_vipss_.vipss_api_.u_v_time);
}

// double acc_tim1;
// static int countopt = 0;


double optfunc_unit_vipss_simple_eigen(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    // printf("input point number : %llu \n", n);
    // Eigen::VectorXd arma_x(n*3);

    // arma_x_opt_g;
    // res_vec_g;

    //(  sin(a)cos(b), sin(a)sin(b), cos(a)  )  a =>[0, pi], b => [-pi, pi];
    // std::vector<double>sina_cosa_sinb_cosb(n * 4);
    #pragma omp parallel for 
    for(int i=0;i<n;++i){
        size_t ind = i*4;
        scsc_vec[ind] = sin(x[i*2]);
        scsc_vec[ind+1] = cos(x[i*2]);
        scsc_vec[ind+2] = sin(x[i*2+1]);
        scsc_vec[ind+3] = cos(x[i*2+1]);
        if(drbf->axi_plane_ == AXI_PlANE::XYZ)
        {
            arma_x_opt_g(i) = scsc_vec[ind] * scsc_vec[ind + 3];
            arma_x_opt_g(i+n) = scsc_vec[ind] * scsc_vec[ind + 2];
            arma_x_opt_g(i+n*2) = scsc_vec[ind + 1];
        } else {
            arma_x_opt_g(i) = scsc_vec[ind] * scsc_vec[ind + 3];
            arma_x_opt_g(i+n) = scsc_vec[ind + 1];
            arma_x_opt_g(i+n*2) = scsc_vec[ind] * scsc_vec[ind + 2];
        }
    }
    // for(int i=0;i<n;++i){
        
    // }

    Eigen::VectorXd a2 = (arma_x_opt_g.transpose() * drbf->local_vipss_.final_h_eigen_).transpose();
    if (!grad.empty()) {
        grad.resize(n*2);
        #pragma omp parallel for 
        for(int i=0;i<n;++i){
            // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
            size_t ind = i*4;
            
            if(drbf->axi_plane_ == AXI_PlANE::XYZ)
            {
                grad[i*2] = a2(i) * scsc_vec[ind + 1] * scsc_vec[ind + 3] + a2(i+n) * scsc_vec[ind + 1] * scsc_vec[ind + 2] - a2(i+n*2) * scsc_vec[ind];
                grad[i*2+1] = -a2(i) * scsc_vec[ind] * scsc_vec[ind + 2] + a2(i+n) * scsc_vec[ind] * scsc_vec[ind + 3];
            } else {
                grad[i*2] = a2(i) * scsc_vec[ind + 1] * scsc_vec[ind + 3] + a2(i+n*2) * scsc_vec[ind + 1] * scsc_vec[ind + 2] - a2(i+n) * scsc_vec[ind];
                grad[i*2+1] = -a2(i) * scsc_vec[ind] * scsc_vec[ind+ 2] + a2(i+n*2) * scsc_vec[ind] * scsc_vec[ind + 3];
            }
            
        }
    }
    double re =  arma_x_opt_g.dot( a2 );
    // printf("residual val : %f \n", re);
    // printf("Final_H_ non zero  : %d \n", drbf->Final_H_.n_nonzero);
    // countopt++;
    // acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);
    drbf->countopt_ ++;
    VIPSSUnit::opt_func_count_g ++;
    return re;
}

double optfunc_unit_vipss_direct_simple_eigen(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    Eigen::VectorXd arma_x(n*3);
    for(int i=0;i<n;++i){
        // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        arma_x(i)     = x[3*i];
        arma_x(i+n)   = x[3*i + 1];
        arma_x(i+n*2) = x[3*i + 2];
    }

    Eigen::VectorXd a2 = (arma_x.transpose() * drbf->local_vipss_.final_h_eigen_).transpose();
    if (!grad.empty()) {
        // std::cout << " grad size " <<grad.size() << std::endl;
        grad.resize(n*3);
        for(size_t i=0;i<n;++i){
            grad[i*3]   = a2(i) * 2;
            grad[i*3+1] = a2(n+i) * 2;
            grad[i*3+2] = a2(2*n+i) * 2;
        }
    }
    double re = arma_x.dot(a2);
    double alpha = 1000.0;

    for(size_t id =0; id < n; ++id)
    {
        double cur_re = x[3*id] * x[3*id] + x[3*id + 1] * x[3*id + 1] + x[3*id + 2] * x[3*id + 2] - 1;
        if(!grad.empty()) 
        {
            grad[3* id]     += alpha * 2 * x[3*id] * cur_re; 
            grad[3* id + 1] += alpha * 2 * x[3*id + 1] * cur_re; 
            grad[3* id + 2] += alpha * 2 * x[3*id + 2] * cur_re; 
        }
        re += alpha* cur_re * cur_re;
    }
    // printf("res val : %f \n", re);
    VIPSSUnit::opt_func_count_g ++;
    return re;
}

double optfunc_unit_vipss_direct_eigen(const std::vector<double>& x, std::vector<double>& grad, void* fdata) {

    auto t1 = Clock::now();
    VIPSSUnit::opt_func_count_g++;
    VIPSSUnit* drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    size_t u_size = 4;
    // Eigen::VectorXd arma_x(n * u_size);
    for (int i = 0; i < n; ++i) {
        // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        arma_x_opt_g(i) = x[u_size * i];
        arma_x_opt_g(i + n) = x[u_size * i + 1];
        arma_x_opt_g(i + n * 2) = x[u_size * i + 2];
        arma_x_opt_g(i + n * 3) = x[u_size * i + 3];
    }
    //Eigen::VectorXd a2 = (arma_x_opt_g.transpose() * drbf->local_vipss_.final_h_eigen_).transpose();

   
    int id;
    //auto x_t = arma_x_opt_g.transpose();
    double alpha = 1.0;
    const auto& final_h = drbf->local_vipss_.final_h_eigen_;
    if (grad.empty()) grad.resize(4 * n);
    #pragma omp parallel for shared(x, alpha, grad, final_h, arma_x_opt_g, res_vec_g) private(id)
    for (id = 0; id < n; ++id)
    {
        res_vec_g[id] = 0;
        for (int j = 0; j < 4; ++j)
        {
            double val = arma_x_opt_g .transpose() * final_h.col(id + j *n);
            res_vec_g[id] += arma_x_opt_g[id + j *n] * val;
            grad[id * 4 + j] = val * 2;
        }
        double cur_re = x[4 * id + 1] * x[4 * id + 1] + x[4 * id + 2] * x[4 * id + 2]
            + x[4 * id + 3] * x[4 * id + 3] - 1;
        grad[4 * id + 1] += alpha * 2 * x[4 * id + 1] * cur_re;
        grad[4 * id + 2] += alpha * 2 * x[4 * id + 2] * cur_re;
        grad[4 * id + 3] += alpha * 2 * x[4 * id + 3] * cur_re;
        res_vec_g[id] += alpha * cur_re * cur_re;
    }
    double re = res_vec_g.sum();

    auto t2 = Clock::now();
    double opt_time = std::chrono::nanoseconds(t2 - t1).count()/ 1e9;
    VIPSSUnit::opt_func_time_g += opt_time;

    // printf("opt fun call time accu: %f \n", VIPSSUnit::opt_func_time_g);

    return re;
}

double optfunc_unit_vipss(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    size_t n = drbf->npt_;
    Eigen::VectorXd arma_x(n*4);
    std::vector<double>sina_cosa_sinb_cosb(n * 4);
    for(size_t i=0;i<n;++i){
        size_t ind = i*4;
        sina_cosa_sinb_cosb[ind]   = sin(x[i*3 + 1]);
        sina_cosa_sinb_cosb[ind+1] = cos(x[i*3 + 1]);
        sina_cosa_sinb_cosb[ind+2] = sin(x[i*3 + 2]);
        sina_cosa_sinb_cosb[ind+3] = cos(x[i*3 + 2]);
    }
    for(int i=0;i<n;++i){
        auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        arma_x(i) = x[i *3];
        arma_x(i + n) = p_scsc[0] * p_scsc[3];
        arma_x(i + n * 2) = p_scsc[0] * p_scsc[2];
        arma_x(i + n * 3) = p_scsc[1];
    }

    Eigen::VectorXd a2 =  (arma_x.transpose() * drbf ->local_vipss_.final_h_eigen_).transpose();
    if (!grad.empty()) {
        grad.resize(n*3);
        for(size_t i=0;i<n;++i){
            auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
            grad[i*3] = 2 * a2[i];
            grad[i*3 + 1] = 2 * (a2(n + i) * p_scsc[1] * p_scsc[3] + a2(i+ n * 2) * p_scsc[1] * p_scsc[2] - a2(i+n*3) * p_scsc[0]);
            grad[i*3 + 2] = 2* (-a2(i + n) * p_scsc[0] * p_scsc[2] + a2(i+n* 2) * p_scsc[0] * p_scsc[3]);
        }
    }
    double re =  arma_x.dot( a2 );
    // printf("residual val : %f \n", re);
    // printf("Final_H_ non zero  : %d \n", drbf->Final_H_.n_nonzero);
    // countopt++;
    // acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);
    drbf->countopt_ ++;
    VIPSSUnit::opt_func_count_g ++;
    
    return re;
}

void VIPSSUnit::OptUnitVipssNormalSimple(){

    printf("start to call solver ! \n");
    solver_.solveval.resize(npt_ * 2);

    for(size_t i=0;i<npt_;++i){
        double *veccc = initnormals_.data()+i*3;
        {
            if(axi_plane_ == AXI_PlANE::XYZ)
            {
                solver_.solveval[i*2] = atan2(sqrt(veccc[0]*veccc[0]+veccc[1]*veccc[1]),veccc[2] );
                solver_.solveval[i*2 + 1] = atan2( veccc[1], veccc[0]);
            } else {
                solver_.solveval[i*2] = atan2(sqrt(veccc[0]*veccc[0]+veccc[2]*veccc[2]),veccc[1] );
                solver_.solveval[i*2 + 1] = atan2( veccc[2], veccc[0]);
            }
            // solver_.solveval[i*2] = atan2(sqrt(veccc[1]*veccc[1]+veccc[2]*veccc[2]),veccc[0] );
            // solver_.solveval[i*2 + 1] = atan2( veccc[2], veccc[1]);
        }
    }
    arma_x_opt_g.resize(npt_ * 3);
    res_vec_g.resize(npt_);
    scsc_vec.resize(npt_ * 4);
    // printf("finish init solver ! \n");
    if(1){ 
        std::vector<double>upper(npt_*2);
        std::vector<double>lower(npt_*2);
        for(int i=0;i<npt_;++i){
            upper[i*2] = 1 * M_PI_;
            upper[i*2 + 1] = 1 * M_PI_;

            lower[i*2] = -1 * M_PI_;
            lower[i*2 + 1] = -1 * M_PI_;
        }
        // countopt = 0;
        // acc_time = 0;
        //LocalIterativeSolver(sol,kk==0?normals:newnormals,300,1e-7);
        // printf("start the solver ! \n");
        Solver::nloptwrapper(lower,upper,optfunc_unit_vipss_simple_eigen,this,opt_tor_, max_opt_iter_,solver_);
        // callfunc_time = acc_time;
        // solve_time = sol.time;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;
    }
    newnormals_.resize(npt_*3);
    s_func_vals_.resize(npt_, 0); 
    // printf("-----------newnormal size : %lu \n", newnormals_.size());
    arma::vec y(npt_ + 3 * npt_);
    for(size_t i=0;i<npt_;++i)y(i) = 0;
    for(size_t i=0;i<npt_;++i){
        double a = solver_.solveval[i*2], b = solver_.solveval[i*2+1];

        if(axi_plane_ == AXI_PlANE::XYZ)
        {
            newnormals_[i*3]   = sin(a) * cos(b);
            newnormals_[i*3+1] = sin(a) * sin(b);
            newnormals_[i*3+2] = cos(a);
        } else {
            newnormals_[i*3]   = sin(a) * cos(b);
            newnormals_[i*3+1] = cos(a);
            newnormals_[i*3+2] = sin(a) * sin(b);
        }
        
        // MyUtility::normalize(newnormals.data()+i*3);
    }
    // Set_RBFCoef(y);
    //sol.energy = arma::dot(a,M*a);
    // if(open_debug_log)
    // cout<<"Opt_Hermite_PredictNormal_UnitNormal"<<endl;
    return;
}

void VIPSSUnit::OptUnitVipssNormalDirectSimple(){

    printf("start to call solver ! \n");
    solver_.solveval.resize(npt_ * 3);

    for(size_t i=0;i<npt_;++i){
        double *veccc = initnormals_.data()+i*3;
        {
            solver_.solveval[i*3] =veccc[0];
            solver_.solveval[i*3 + 1] = veccc[1];
            solver_.solveval[i*3 + 2] = veccc[2];
        }
    }
    std::vector<double>upper(npt_*3);
    std::vector<double>lower(npt_*3);
    for(int i=0;i<npt_;++i){
        upper[i*3] = 1;
        upper[i*3 + 1] = 1;
        upper[i*3 + 2] = 1;

        lower[i*3] = -1.0;
        lower[i*3 + 1] = -1.0;
        lower[i*3 + 2] = -1.0;
    }

    Solver::nloptwrapperDirect(lower,upper,optfunc_unit_vipss_direct_simple_eigen,
                this,opt_tor_, max_opt_iter_,solver_);
    // Solver::nloptwrapper(lower,upper,optfunc_unit_vipss_simple,this,1e-7,3000,solver_);
    
    newnormals_.resize(npt_*3);
    s_func_vals_.resize(npt_, 0);
    arma::vec y(npt_ + 3 * npt_);
    for(size_t i=0;i<npt_;++i)y(i) = 0;
    for(size_t i=0;i<npt_;++i){
        newnormals_[i*3]   = solver_.solveval[i*3];
        newnormals_[i*3+1] = solver_.solveval[i*3 + 1];
        newnormals_[i*3+2] = solver_.solveval[i*3 + 2];
        // MyUtility::normalize(newnormals.data()+i*3);
    }
    return;
}


void VIPSSUnit::OptUnitVipssNormalDirect(){

    printf("start to call solver ! \n");

    size_t u_size = 4;
    solver_.solveval.resize(npt_ * u_size);

    for(size_t i=0;i<npt_;++i){
        double *veccc = initnormals_.data()+i*3;
        {
            solver_.solveval[i*u_size]     = 0;
            solver_.solveval[i*u_size + 1] = veccc[0];
            solver_.solveval[i*u_size + 2] = veccc[1];
            solver_.solveval[i*u_size + 3] = veccc[2];
        }
    }
    std::vector<double>upper(npt_ * u_size);
    std::vector<double>lower(npt_ * u_size);
    for(int i=0;i<npt_;++i){
        for(size_t j = 0; j < u_size; ++j)
        {
            upper[i*u_size + j] = 1;
            lower[i*u_size + j] = -1.0;
        }
    }
    arma_x_opt_g.resize(npt_ * 4);
    res_vec_g.resize(npt_);
    Solver::nloptwrapperDirect(lower,upper,optfunc_unit_vipss_direct_eigen,
                this, opt_tor_, max_opt_iter_, solver_);
    
    // Solver::nloptwrapper(lower,upper,optfunc_unit_vipss_simple,this,1e-7,3000,solver_);
    newnormals_.resize(npt_*3);
    s_func_vals_.resize(npt_);
    // arma::vec y(npt_ + 3 * npt_);
    // for(size_t i=0;i<npt_;++i)y(i) = 0;
    for(size_t i=0;i<npt_;++i){
        s_func_vals_[i]    = solver_.solveval[i*u_size];
        newnormals_[i*3]   = solver_.solveval[i*u_size + 1];
        newnormals_[i*3+1] = solver_.solveval[i*u_size + 2];
        newnormals_[i*3+2] = solver_.solveval[i*u_size + 3];
        // MyUtility::normalize(newnormals.data()+i*3);
    }
    return;
}


void VIPSSUnit::OptUnitVipssNormal(){

    printf("start to call solver ! \n");
    solver_.solveval.resize(npt_ * 3);

    for(size_t i=0;i<npt_;++i){
        double *veccc = initnormals_.data()+i*3;
        {
            solver_.solveval[i*3] = 0.0;
            solver_.solveval[i*3 + 1] = atan2(sqrt(veccc[0]*veccc[0]+veccc[1]*veccc[1]),veccc[2] );
            solver_.solveval[i*3 + 2] = atan2( veccc[1], veccc[0]   );
        }
    }
    arma_x_opt_g.resize(npt_ * 4);
    res_vec_g.resize(npt_);
    // printf("finish init solver ! \n");
    if(1){
        std::vector<double>upper(npt_*3);
        std::vector<double>lower(npt_*3);
        for(int i=0;i<npt_;++i){
            upper[i*3 ] = 1.0;
            upper[i*3 + 1] = 1 * M_PI_;
            upper[i*3 + 2] = 1 * M_PI_;

            lower[i*3] = -1.0;
            lower[i*3 + 1] = -1.0 * M_PI_;
            lower[i*3 + 2] = -1 * M_PI_;
        }
        // printf("start the solver ! \n");
        Solver::nloptwrapper(lower,upper,optfunc_unit_vipss,this,opt_tor_, max_opt_iter_,solver_);
    }
    newnormals_.resize(npt_*3);
    s_func_vals_.resize(npt_);
    for(size_t i=0;i<npt_;++i) s_func_vals_[i] = solver_.solveval[i*3];;
    for(size_t i=0;i<npt_;++i){
        double a = solver_.solveval[i*3 + 1], b = solver_.solveval[i*3+2];
        newnormals_[i*3]   = sin(a) * cos(b);
        newnormals_[i*3+1] = sin(a) * sin(b);
        newnormals_[i*3+2] = cos(a);
        // MyUtility::normalize(newnormals.data()+i*3);
    }
    // Set_RBFCoef(y);
    //sol.energy = arma::dot(a,M*a);
    // if(open_debug_log)
    // cout<<"Opt_Hermite_PredictNormal_UnitNormal"<<endl;
    return;
}

void VIPSSUnit::BuildNNHRBFFunctions()
{
    auto t000 = Clock::now();
    std::vector<std::array<double,3>> octree_sample_pts;
    if(make_nn_const_neighbor_num_)
    {
        SimOctree::SimpleOctree octree;
        // std::cout << " start to init octree " << std::endl;
        octree.InitOctTree(local_vipss_.origin_in_pts_, 4);
        std::cout << " insert octree center pts num : " << octree.octree_centers_.size() << std::endl; 
        octree_sample_pts = octree.octree_centers_;
    }

    local_vipss_.voro_gen_.GenerateVoroData();
    local_vipss_.voro_gen_.SetInsertBoundaryPtsToUnused();
    local_vipss_.voro_gen_.BuildTetMeshTetCenterMap();
    local_vipss_.voro_gen_.BuildPicoTree();
    // auto boundary_pts = local_vipss_.voro_gen_.insert_boundary_pts_;
    local_vipss_.voro_gen_.insert_boundary_pts_.clear();
    auto t001 = Clock::now();
    G_VP_stats.generate_voro_data_time_ = std::chrono::nanoseconds(t001 - t000).count() / 1e9;
    // G_VP_stats.generate_voro_data_time_ = generate_voro_data_time
    
    // auto t002 = Clock::now();
    // double generate_voro_data_time = std::chrono::nanoseconds(t002 - t001).count() / 1e9;
    local_vipss_.normals_ = newnormals_;
    local_vipss_.s_vals_ = s_func_vals_;
    local_vipss_.user_lambda_ = user_lambda_;
    local_vipss_.BuildHRBFPerNode();
    local_vipss_.SetThis(); 
    local_vipss_.dummy_sign_ = LocalVipss::NNDistFunction(R3Pt(local_vipss_.voro_gen_.dummy_sign_pt_[0],
                                                               local_vipss_.voro_gen_.dummy_sign_pt_[1], 
                                                               local_vipss_.voro_gen_.dummy_sign_pt_[2]));    
    if(abs(local_vipss_.dummy_sign_ ) > 1e-18)
    {
        local_vipss_.dummy_sign_ = local_vipss_.dummy_sign_ / abs(local_vipss_.dummy_sign_);
    }
    std::cout << " ********** dummy pt sign val : " << local_vipss_.dummy_sign_ << std::endl;

    
    if(make_nn_const_neighbor_num_)
    {
        std::vector<double> insert_pt_func_vals;
        std::vector<double> insert_pt_func_gradients;
        std::vector<std::array<double,3>> valid_pts;
        std::vector<double> dummy_dist_vals; 
        auto t001 = Clock::now();
        // std::string octree_sample_path = out_dir_ + file_name_ + "octree_sample.xyz";
        // std::ofstream octree_file(octree_sample_path);
    if(0)
    {
        for(auto pt : octree_sample_pts)
        {
            double gradient[3];
            double dist_val = local_vipss_.NatureNeighborGradientOMP(&pt[0], gradient);
            local_vipss_.s_vals_.push_back(dist_val);
            local_vipss_.normals_.push_back(-1.0 * gradient[0] );
            local_vipss_.normals_.push_back(-1.0 * gradient[1] );
            local_vipss_.normals_.push_back(-1.0 * gradient[2] );
            // octree_file << pt[0] << " " << pt[1] << " " << pt[2] << " ";
            // octree_file << gradient[0] << " " << gradient[1] << " " << gradient[2] << std::endl;
        }
    }
        auto t0022 = Clock::now();
        G_VP_stats.octree_pt_gradient_cal_time_ = std::chrono::nanoseconds(t0022 - t001).count() / 1e9;
        std::cout << "evaluate octree sample time : " << G_VP_stats.octree_pt_gradient_cal_time_ << std::endl;

        G_VP_stats.octree_dummy_pt_num_ = octree_sample_pts.size();
        local_vipss_.voro_gen_.InsertPts(octree_sample_pts);
        local_vipss_.voro_gen_.BuildTetMeshTetCenterMap();
        local_vipss_.voro_gen_.BuildPicoTree();
        local_vipss_.voro_gen_.voronoi_data_.clean_memory();
        local_vipss_.voro_gen_.tetMesh_.generate_voronoi_cell(&(local_vipss_.voro_gen_.voronoi_data_));
        local_vipss_.voro_gen_.SetInsertBoundaryPtsToUnused();
        auto t0033 = Clock::now();
        double rebuild_pico_and_voro_time = std::chrono::nanoseconds(t0033 - t0022).count() / 1e9;
        std::cout << "rebuild_pico_and_voro_time : " << rebuild_pico_and_voro_time << std::endl;

    if(0)
    {
        const auto& insert_pts = local_vipss_.voro_gen_.insert_boundary_pts_;
        int input_pt_size = local_vipss_.points_.size();
        // std::vector<tetgenmesh::point> octree_insert_pts;
        for(auto pt : insert_pts)
        {
            int new_pid = local_vipss_.points_.size();
            local_vipss_.points_.push_back(pt);
            local_vipss_.voro_gen_.point_id_map_[pt] = new_pid;
        }
        auto all_valid_pt_size = local_vipss_.points_.size();
        std::cout << "octree_insert_pts size : " << insert_pts.size() << std::endl;
        int cluster_sum = 0;
        local_vipss_.voro_gen_.cluster_init_pids_.clear();
        local_vipss_.voro_gen_.cluster_init_pids_.resize(all_valid_pt_size);
        // local_vipss_.node_rbf_vec_.clear();
        local_vipss_.node_rbf_vec_.resize(all_valid_pt_size);

        auto t0041 = Clock::now();
        VoronoiGen::cluster_init_pts_.resize(all_valid_pt_size);
        std::vector<std::vector<double>> cluster_sv_vecs(all_valid_pt_size); 
        std::vector<std::vector<double>> cluster_pt_vecs(all_valid_pt_size);
        std::vector<std::vector<double>> cluster_nl_vecs(all_valid_pt_size);

        int max_cluster_size = 0;
        for(int i =0; i < all_valid_pt_size; ++i)
        {
            auto cur_pt = local_vipss_.points_[i];
            std::set<tetgenmesh::point> candidate_pts;
            local_vipss_.voro_gen_.GetVertexStar(cur_pt, candidate_pts, 1);
            std::vector<size_t> cluster_pt_ids;
            // cluster_pt_ids.push_back(local_vipss_.voro_gen_.point_id_map_[cur_pt]);
            // std::vector<double> cluster_nl_vec;
            auto cur_pid = local_vipss_.voro_gen_.point_id_map_[cur_pt];
            std::vector<double> cluster_nl_vec;
            std::vector<double> cluster_sv_vec; 
            std::vector<double> cluster_pt_vec;
            for(auto nn_pt : candidate_pts)
            {
                // if( nn_pt == cur_pt) continue;
                auto pid = local_vipss_.voro_gen_.point_id_map_[nn_pt];
                cluster_pt_ids.push_back(pid);
                cluster_pt_vec.push_back(nn_pt[0]);
                cluster_pt_vec.push_back(nn_pt[1]);
                cluster_pt_vec.push_back(nn_pt[2]);
                // cluster_sv_vec.push_back(local_vipss_.s_vals_[pid]);
                // cluster_nl_vec.push_back(local_vipss_.normals_[3*pid]);
                // cluster_nl_vec.push_back(local_vipss_.normals_[3*pid + 1]);
                // cluster_nl_vec.push_back(local_vipss_.normals_[3*pid + 2]);
            }
            max_cluster_size = max_cluster_size > cluster_pt_ids.size() ? max_cluster_size : cluster_pt_ids.size();
            cluster_sum += cluster_pt_ids.size();
            local_vipss_.voro_gen_.cluster_init_pids_[i] = cluster_pt_ids;
            VoronoiGen::cluster_init_pts_[i] = cluster_pt_vec;
            // cluster_sv_vecs[i] = cluster_sv_vec;
            // cluster_pt_vecs[i] = cluster_pt_vec; 
            // cluster_nl_vecs[i] = cluster_nl_vec;
        }
        std::cout << " max_cluster_size : " << max_cluster_size << std::endl;
        auto t0042 = Clock::now();
        double get_pt_neigbor_time = std::chrono::nanoseconds(t0042 - t0041).count() / 1e9;
        std::cout << "get_pt_neigbor_time : " << get_pt_neigbor_time << std::endl;
        local_vipss_.BuildHRBFPerNode();
        local_vipss_.SetThis();

        auto t0043 = Clock::now();
        double rebuild_neigbor_time = std::chrono::nanoseconds(t0043 - t0042).count() / 1e9;
        std::cout << "rebuild_NN HRBF time : " << rebuild_neigbor_time << std::endl;

    if(0)
    {
        #pragma omp parallel for 
        for(int i =0; i < all_valid_pt_size; ++i)
        {
            // auto cur_pt = local_vipss_.points_;
            // const auto& cluster_pt_ids = local_vipss_.voro_gen_.cluster_init_pids_[i]; 
            // int cluster_size = cluster_pt_ids.size();
            auto& cluster_nl_vec = cluster_nl_vecs[i];
            auto& cluster_sv_vec = cluster_sv_vecs[i]; 
            auto& cluster_pt_vec = cluster_pt_vecs[i];
            // std::vector<double> cluster_nl_vec (cluster_size * 3);
            // std::vector<double> cluster_sv_vec(cluster_size); 
            // std::vector<double> cluster_pt_vec(cluster_size * 3);
            // for(int j = 0; j < cluster_pt_ids.size(); ++j)
            {
                // auto pid = cluster_pt_ids[j];
                // cluster_pt_vec[3*j]     = local_vipss_.points_[pid][0];
                // cluster_pt_vec[3*j + 1] = local_vipss_.points_[pid][1];
                // cluster_pt_vec[3*j + 2] = local_vipss_.points_[pid][2];
                // cluster_sv_vec[j]       = local_vipss_.s_vals_[pid];
                // cluster_nl_vec[3*j]     = local_vipss_.normals_[3*pid];
                // cluster_nl_vec[3*j + 1] = local_vipss_.normals_[3*pid + 1];
                // cluster_nl_vec[3*j + 2] = local_vipss_.normals_[3*pid + 2];
                // cluster_pt_vec.push_back(local_vipss_.points_[pid][0]);
                // cluster_pt_vec.push_back(local_vipss_.points_[pid][1]);
                // cluster_pt_vec.push_back(local_vipss_.points_[pid][2]);
                // cluster_sv_vec.push_back(local_vipss_.s_vals_[pid]);
                // cluster_nl_vec.push_back(local_vipss_.normals_[3*pid]);
                // cluster_nl_vec.push_back(local_vipss_.normals_[3*pid + 1]);
                // cluster_nl_vec.push_back(local_vipss_.normals_[3*pid + 2]);
            }
            local_vipss_.node_rbf_vec_[i] = std::make_shared<RBF_Core>();
            local_vipss_.vipss_api_.build_cluster_hrbf(cluster_pt_vec, cluster_nl_vec, cluster_sv_vec, local_vipss_.node_rbf_vec_[i]);
        }
    }
    
        int ave_cluster_size = int(double(cluster_sum) / double(local_vipss_.points_.size()));
        std::cout << " ------ ave cluster size " << ave_cluster_size << std::endl;
    }

    }
    auto t003 = Clock::now();
    // double build_nn_rbf_time  
    // G_VP_stats.build_nn_rbf_time_ = std::chrono::nanoseconds(t003 - t001).count() / 1e9;
    G_VP_stats.build_nn_rbf_time_ = std::chrono::nanoseconds(t003 - t000).count() / 1e9;
    G_VP_stats.average_cluster_size_ = local_vipss_.voro_gen_.average_neighbor_num_;
    G_VP_stats.pt_num_ = local_vipss_.points_.size();

     std::cout << " ------ build HRBF time all :  " << G_VP_stats.build_nn_rbf_time_ << std::endl;

    // local_vipss_.voro_gen_.SetInsertBoundaryPtsToUnused();

}

void VIPSSUnit::ReconSurface()
{
    // BuildNNHRBFFunctions();
    printf(" start ReconSurface \n");
    // local_vipss_.TestVoronoiPts();
    // local_vipss_.voro_gen_.BuildTetMeshTetCenterMap();
    // printf("init pico tree time  : %f ! \n", init_pico_tree_time);
    // bool use_nn_interpolation = true;
    // auto t002 = Clock::now();

    if(LOCAL_HRBF_NN == hrbf_type_)
    {
        auto t000 = Clock::now();
        // double build_HRBF_time = std::chrono::nanoseconds(t003 - t002).count() / 1e9;
        // G_VP_stats.build_per_cluster_hrbf_total_time_ += build_HRBF_time;     
        // printf(" start set local vipss static ptr \n");
        // local_vipss_.SetThis();
        size_t n_voxels_1d = volume_dim_;
        Surfacer sf;
        auto surf_time = sf.Surfacing_Implicit(local_vipss_.out_pts_, n_voxels_1d, false, LocalVipss::NNDistFunction);
        auto t004 = Clock::now();
        sf.WriteSurface(finalMesh_v_,finalMesh_fv_);
        // std::string out_path = data_dir_ + "/" + file_name_  + "/nn_surface";
        auto t005 = Clock::now();
        double surface_file_save_time = std::chrono::nanoseconds(t005 - t004).count() / 1e9;
        double total_surface_time = std::chrono::nanoseconds(t004 - t000).count() / 1e9;
        writePLYFile_VF(out_surface_path_, finalMesh_v_, finalMesh_fv_);
        std::cout << "------- tet search time "<< tetgenmesh::tet_search_time_st << std::endl;
        std::cout << "------- voxel pt ave nn num "<< LocalVipss::ave_voxel_nn_pt_num_ / LocalVipss::DistCallNum << std::endl;
        printf("------- nn search time: %f \n", local_vipss_.search_nn_time_sum_);
        printf("------- cal nn coordinate and hrbf time: %f \n", local_vipss_.pass_time_sum_);
        printf(" ------ voxel dist func val evaluated count : %d  \n", LocalVipss::DistCallNum);
        printf(" ------ voxel dist func val evaluated time : %f \n", LocalVipss::DistCallTime);
        printf("------- total surface time: %f \n", total_surface_time);
        G_VP_stats.neighbor_search_time_ += local_vipss_.search_nn_time_sum_;
        G_VP_stats.cal_nn_coordinate_and_hbrf_time_ += local_vipss_.pass_time_sum_;
        G_VP_stats.voxel_cal_num += LocalVipss::DistCallNum;
        G_VP_stats.nn_evaluate_count_ = LocalVipss::DistCallNum;
        G_VP_stats.average_neighbor_num_ = double(LocalVipss::ave_voxel_nn_pt_num_)/ double(LocalVipss::DistCallNum);
        // G_VP_stats.surface_total_time_ += total_surface_time;
        G_VP_stats.surface_total_time_ = total_surface_time;

    } else {
        rbf_api_.user_lambda_ = user_lambda_;
        // rbf_api_.user_lambda_ = 0.001;
        rbf_api_.outpath_ = data_dir_ + file_name_ + "/";
        rbf_api_.is_surfacing_ = true;
        rbf_api_.n_voxel_line_ = volume_dim_;
        rbf_api_.run_vipss(local_vipss_.out_pts_, newnormals_, s_func_vals_);
    }

    // if(use_hrbf_surface_)
    // {
    //     rbf_api_.user_lambda_ = user_lambda_;
    //     // rbf_api_.user_lambda_ = 0.001;
    //     rbf_api_.outpath_ = data_dir_ + file_name_ + "/";
    //     rbf_api_.is_surfacing_ = true;
    //     rbf_api_.n_voxel_line_ = volume_dim_;
    //     rbf_api_.run_vipss(local_vipss_.out_pts_, newnormals_, s_func_vals_);
    // } 
    // local_vipss_.
}


void VIPSSUnit::GenerateAdaptiveGrid()
{
    auto t000 = Clock::now();   
    // std::cout << " test val " << test_val << std::endl;
    std::array<size_t,3> resolution = {3, 3, 3};
    GenerateAdaptiveGridOut(resolution, local_vipss_.voro_gen_.bbox_min_, 
                            local_vipss_.voro_gen_.bbox_max_, out_dir_,  file_name_, adgrid_threshold_);
    auto t001 = Clock::now();
    double adgrid_gen_time = std::chrono::nanoseconds(t001 - t000).count() / 1e9;
    printf("adaptive grid generation time : %f ! \n", adgrid_gen_time);
}

void VIPSSUnit::SolveOptimizaiton()
{
    // auto tn0 = Clock::now();
    
    // auto tn1 = Clock::now();
    // double get_h_sub_block_time = std::chrono::nanoseconds(tn1 - tn0).count() / 1e9;
    // // printf("get H sub block time : %f ! \n", get_h_sub_block_time);
    // G_VP_stats.take_h_sub_block_time_ = get_h_sub_block_time;
    auto ts0 = Clock::now();
    Solver::open_log_ = true;
    if(user_lambda_ < 1e-12)
    {
        if (hard_constraints_)
        {
            OptUnitVipssNormalSimple();
        }
        else {
            OptUnitVipssNormalDirectSimple();
        }
    } else {
        if (hard_constraints_)
        {
            OptUnitVipssNormal();
        }
        else {
            OptUnitVipssNormalDirect();
        }
    }
    Final_H_.clear();
    local_vipss_.final_H_.clear();
    
    // printf("opt fun call time : %f \n", VIPSSUnit::opt_func_time_g);
#pragma omp parallel for
    for(int i = 0; i < newnormals_.size()/3; ++i)
    {
        double normal_len = sqrt(newnormals_[3*i] * newnormals_[3*i] 
        + newnormals_[3*i + 1] * newnormals_[3*i + 1] + newnormals_[3*i + 2] * newnormals_[3*i + 2]);

        newnormals_[3*i] /= normal_len;
        newnormals_[3*i + 1] /= normal_len;
        newnormals_[3*i + 2] /= normal_len;
    }

    auto ts1 = Clock::now();
    double solve_time = std::chrono::nanoseconds(ts1 - ts0).count() / 1e9;
    printf("opt solve time : %f ! \n", solve_time);
    printf("opt fun call count : %d \n", VIPSSUnit::opt_func_count_g);
    G_VP_stats.opt_solver_time_ = solve_time;
    G_VP_stats.opt_func_call_num_ = VIPSSUnit::opt_func_count_g;
}

void VIPSSUnit::Run()
{
    local_vipss_.out_dir_ = data_dir_ + "/" + file_name_ + "/";
    auto t00 = Clock::now();
    rbf_api_.Set_RBF_PARA();
    InitPtNormalWithLocalVipss();
    SolveOptimizaiton();
    BuildNNHRBFFunctions();
    auto new_pts = local_vipss_.octree_leaf_pts_;
    if(LocalVipss::use_octree_sample_)
    {
        auto ts00 = Clock::now();
        const auto&pts = local_vipss_.octree_split_leaf_pts_;
        // const auto&pts = local_vipss_.octree_leaf_pts_;
        std::cout << "split pts size : " << pts.size()/3 << std::endl;
        for(int i =0; i < pts.size()/3; ++i)
        {
            // printf("cur pt  : %f %f %f \n", pts[3*i], pts[3*i + 1], pts[3*i + 2] );
            R3Pt cur_pt(pts[3*i], pts[3*i + 1], pts[3*i + 2]);
            double cur_dist = LocalVipss::NNDistFunction(cur_pt);
            // printf("cur dist : %f \n", cur_dist );
            if(abs(cur_dist) > distfunc_threshold_)
            {
                new_pts.push_back(pts[3*i]);
                new_pts.push_back(pts[3*i + 1]);
                new_pts.push_back(pts[3*i + 2]);
            }
        }
        auto ts01 = Clock::now();
        double total_time = std::chrono::nanoseconds(ts01 - ts00).count()/1e9;
        printf("------- remaining pts dist function evaluation time : %f ! \n", total_time);
        // std::string octree_sample_path = out_dir_  + file_name_ +  "_octree_distSample.xyz";
        // writeXYZ(octree_sample_path, new_pts);
    }
    auto t01 = Clock::now();
    double total_time = std::chrono::nanoseconds(t01 - t00).count()/1e9;
    printf("total local vipss running time : %f ! \n", total_time);
    // std::string out_path  = local_vipss_.out_dir_ + local_vipss_.filename_  + "_opt";
    writePLYFile_VN(out_normal_path_, local_vipss_.out_pts_, newnormals_);

    // is_surfacing_ = false;
    if (is_surfacing_)
    {
        if(use_adgrid_)
        {
            GenerateAdaptiveGrid();
        } else {
            ReconSurface();
        }  
        local_vipss_.voro_gen_.SetInsertBoundaryPtsToUnused();
        // if(make_nn_const_neighbor_num_)
        // local_vipss_.voro_gen_.SetInsertBoundaryPtsToUnused();
    }
    
    // test_vipss_timing::test_local_vipss(input_data_path_);
    // test_vipss_timing::visual_distval_pt(input_data_path_, 200);
    
    std::string out_csv_file = out_dir_ + file_name_ + "_time_stats.txt";
    WriteStatsTimeCSV(out_csv_file, G_VP_stats);
}

void VIPSSUnit::CalEnergyWithGtNormal()
{
    std::string norm_path = "c:\\Users\\xiaji\\Documents\\projects\\sketches_results\\crab_out_normal_old.ply";
    std::vector<double> vertices;
    std::vector<double> normals;
    readPLYFile(norm_path, vertices, normals);

    if(user_lambda_ <= 1e-12)
    {
        size_t n = vertices.size()/3;
        Eigen::VectorXd arma_x(n*3);
        for(int i=0;i<n;++i){
            // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
            arma_x(i)     = normals[3*i];
            arma_x(i+n)   = normals[3*i + 1];
            arma_x(i+n*2) = normals[3*i + 2];
        }

        Eigen::VectorXd a2 = (arma_x.transpose() * local_vipss_.final_h_eigen_).transpose();
        double re = arma_x.dot(a2);
        std::cout << "final residual val : " << re << std::endl;
    } else {
        size_t n = vertices.size()/3;
        Eigen::VectorXd arma_x(n*4);
        for(int i=0;i<n;++i){
            // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
            arma_x(i)       = s_func_vals_[i];
            arma_x(i + n)   = normals[3*i];
            arma_x(i + n*2) = normals[3*i + 1];
            arma_x(i + n*3) = normals[3*i + 2];
        }

        Eigen::VectorXd a2 = (arma_x.transpose() * local_vipss_.final_h_eigen_).transpose();
        double re = arma_x.dot(a2);
        std::cout << "final residual val : " << re << std::endl;
    }
    

}