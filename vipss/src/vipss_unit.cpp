#include <chrono>
#include "vipss_unit.hpp"
#include "stats.h"

typedef std::chrono::high_resolution_clock Clock;

VP_STATS G_VP_stats;

void VIPSSUnit::InitPtNormalWithLocalVipss()
{
    std::string& data_dir =  data_dir_;
    
    local_vipss_.filename_ = file_name_;
    // l_vp.filename_ = "planck";
    local_vipss_.out_dir_ = data_dir + local_vipss_.filename_ + "/";
    std::string path = data_dir + local_vipss_.filename_ + "/" + local_vipss_.filename_ + ".ply";

    // local_vipss_.angle_threshold_ = 30;
    // local_vipss_.user_lambda_ = user_lambda_;
    local_vipss_.user_lambda_ = user_lambda_;
    // local_vipss_.max_iter_ = 30;
    local_vipss_.use_hrbf_surface_ = false;
    local_vipss_.angle_threshold_ = merge_angle_;
    // local_vipss_.volume_dim_ = 100;
    local_vipss_.Init(path);

    // return;
    if(init_with_cluster_merge_)
    {
        local_vipss_.InitNormalsWithMerge();
    } else {
        local_vipss_.InitNormals();
    }
    local_vipss_.ClearPartialMemory();
    
    local_vipss_.BuildMatrixH();
    npt_ = local_vipss_.points_.size();
    initnormals_ = local_vipss_.out_normals_;

    // printf("unit vipss J mat init time : %f \n", local_vipss_.vipss_api_.u_v_time);
}

// double acc_tim1;
// static int countopt = 0;

double optfunc_unit_vipss_simple(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    // printf("input point number : %llu \n", n);
    arma::vec arma_x(n*3);

    //(  sin(a)cos(b), sin(a)sin(b), cos(a)  )  a =>[0, pi], b => [-pi, pi];
    std::vector<double>sina_cosa_sinb_cosb(n * 4);
    for(size_t i=0;i<n;++i){
        size_t ind = i*4;
        sina_cosa_sinb_cosb[ind] = sin(x[i*2]);
        sina_cosa_sinb_cosb[ind+1] = cos(x[i*2]);
        sina_cosa_sinb_cosb[ind+2] = sin(x[i*2+1]);
        sina_cosa_sinb_cosb[ind+3] = cos(x[i*2+1]);
    }
    for(int i=0;i<n;++i){
        auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;

        if(drbf->axi_plane_ == AXI_PlANE::XYZ)
        {
            arma_x(i) = p_scsc[0] * p_scsc[3];
            arma_x(i+n) = p_scsc[0] * p_scsc[2];
            arma_x(i+n*2) = p_scsc[1];
        } else {
            arma_x(i) = p_scsc[0] * p_scsc[3];
            arma_x(i+n) = p_scsc[1];
            arma_x(i+n*2) = p_scsc[0] * p_scsc[2];
        }
        

        
    }
    arma::vec a2;
    a2 = drbf->Final_H_ * arma_x;
    if (!grad.empty()) {
        grad.resize(n*2);
        for(size_t i=0;i<n;++i){
            auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;

            if(drbf->axi_plane_ == AXI_PlANE::XYZ)
            {
                grad[i*2] = a2(i) * p_scsc[1] * p_scsc[3] + a2(i+n) * p_scsc[1] * p_scsc[2] - a2(i+n*2) * p_scsc[0];
                grad[i*2+1] = -a2(i) * p_scsc[0] * p_scsc[2] + a2(i+n) * p_scsc[0] * p_scsc[3];
            } else {
                grad[i*2] = a2(i) * p_scsc[1] * p_scsc[3] + a2(i+n*2) * p_scsc[1] * p_scsc[2] - a2(i+n) * p_scsc[0];
                grad[i*2+1] = -a2(i) * p_scsc[0] * p_scsc[2] + a2(i+n*2) * p_scsc[0] * p_scsc[3];
            }
        }
    }
    double re = arma::dot( arma_x, a2 );
    // printf("residual val : %f \n", re);
    // printf("Final_H_ non zero  : %d \n", drbf->Final_H_.n_nonzero);
    // countopt++;
    // acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);
    return re;
}


double optfunc_unit_vipss_simple_eigen(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    // printf("input point number : %llu \n", n);
    Eigen::VectorXd arma_x(n*3);

    //(  sin(a)cos(b), sin(a)sin(b), cos(a)  )  a =>[0, pi], b => [-pi, pi];
    std::vector<double>sina_cosa_sinb_cosb(n * 4);
    for(size_t i=0;i<n;++i){
        size_t ind = i*4;
        sina_cosa_sinb_cosb[ind] = sin(x[i*2]);
        sina_cosa_sinb_cosb[ind+1] = cos(x[i*2]);
        sina_cosa_sinb_cosb[ind+2] = sin(x[i*2+1]);
        sina_cosa_sinb_cosb[ind+3] = cos(x[i*2+1]);
    }
    for(int i=0;i<n;++i){
        auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;

        if(drbf->axi_plane_ == AXI_PlANE::XYZ)
        {
            arma_x(i) = p_scsc[0] * p_scsc[3];
            arma_x(i+n) = p_scsc[0] * p_scsc[2];
            arma_x(i+n*2) = p_scsc[1];
        } else {
            arma_x(i) = p_scsc[0] * p_scsc[3];
            arma_x(i+n) = p_scsc[1];
            arma_x(i+n*2) = p_scsc[0] * p_scsc[2];
        }
    }

    Eigen::VectorXd a2 = drbf->local_vipss_.final_h_eigen_ * arma_x;
    if (!grad.empty()) {
        grad.resize(n*2);
        for(size_t i=0;i<n;++i){
            auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;

            if(drbf->axi_plane_ == AXI_PlANE::XYZ)
            {
                grad[i*2] = a2(i) * p_scsc[1] * p_scsc[3] + a2(i+n) * p_scsc[1] * p_scsc[2] - a2(i+n*2) * p_scsc[0];
                grad[i*2+1] = -a2(i) * p_scsc[0] * p_scsc[2] + a2(i+n) * p_scsc[0] * p_scsc[3];
            } else {
                grad[i*2] = a2(i) * p_scsc[1] * p_scsc[3] + a2(i+n*2) * p_scsc[1] * p_scsc[2] - a2(i+n) * p_scsc[0];
                grad[i*2+1] = -a2(i) * p_scsc[0] * p_scsc[2] + a2(i+n*2) * p_scsc[0] * p_scsc[3];
            }
            

            
        }
    }
    double re =  arma_x.dot( a2 );
    // printf("residual val : %f \n", re);
    // printf("Final_H_ non zero  : %d \n", drbf->Final_H_.n_nonzero);
    // countopt++;
    // acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);
    return re;
}

double optfunc_unit_vipss_direct_simple(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    arma::vec arma_x(n*3);
    for(int i=0;i<n;++i){
        // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        arma_x(i)     = x[3*i];
        arma_x(i+n)   = x[3*i + 1];
        arma_x(i+n*2) = x[3*i + 2];
    }
    arma::vec a2;
    a2 = drbf->Final_H_ * arma_x;
    if (!grad.empty()) {
        // std::cout << " grad size " <<grad.size() << std::endl;
        grad.resize(n*3);
        for(size_t i=0;i<n;++i){
            grad[i*3]   = a2(i);
            grad[i*3+1] = a2(n+i);
            grad[i*3+2] = a2(2*n+i);
        }
    }
    double re = arma::dot(arma_x, a2);
    double alpha = 100.0;

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
    // arma::vec a2;
    // Eigen::VectorXd a2 = drbf->local_vipss_.final_h_eigen_ * arma_x;
    Eigen::VectorXd a2 = (arma_x.transpose() * drbf->local_vipss_.final_h_eigen_).transpose();
    if (!grad.empty()) {
        // std::cout << " grad size " <<grad.size() << std::endl;
        grad.resize(n*3);
        for(size_t i=0;i<n;++i){
            grad[i*3]   = a2(i);
            grad[i*3+1] = a2(n+i);
            grad[i*3+2] = a2(2*n+i);
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
    VIPSSUnit::opt_func_count ++;
    return re;
}


double optfunc_unit_vipss_direct(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    size_t u_size = 4;
    arma::vec arma_x(n * u_size);
    for(int i=0;i<n;++i){
        // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        arma_x(i)     = x[u_size*i ];
        arma_x(i+n)   = x[u_size*i + 1];
        arma_x(i+n*2) = x[u_size*i + 2];
        arma_x(i+n*3) = x[u_size*i + 3];
    }
    arma::vec a2;
    a2 = drbf->Final_H_ * arma_x;
    if (!grad.empty()) {
        // std::cout << " grad size " <<grad.size() << std::endl;
        grad.resize(n*u_size);
        for(size_t i=0;i<n;++i){
            grad[i*u_size]   = a2(i);
            grad[i*u_size+1] = a2(i + n);
            grad[i*u_size+2] = a2(i + 2*n);
            grad[i*u_size+3] = a2(i + 3*n);
        }
    }
    double re = arma::dot(arma_x, a2);
    double alpha = 100.0;

    for(size_t id =0; id < n; ++id)
    {
        double cur_re = x[u_size*id + 1] * x[u_size*id + 1] + x[u_size*id + 2] * x[u_size*id + 2] 
                        + x[u_size*id + 3] * x[u_size*id + 3] - 1;
        if(!grad.empty()) 
        {
            grad[u_size* id + 1] += alpha * 2 * x[u_size*id + 1] * cur_re; 
            grad[u_size* id + 2] += alpha * 2 * x[u_size*id + 2] * cur_re; 
            grad[u_size* id + 3] += alpha * 2 * x[u_size*id + 3] * cur_re; 
        }
        re += alpha* cur_re * cur_re;
    }
    // printf("res val : %f \n", re);
    return re;
}

int VIPSSUnit::opt_func_count = 0;


double optfunc_unit_vipss_direct_eigen(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    size_t u_size = 4;
    Eigen::VectorXd arma_x(n * u_size);
    for(int i=0;i<n;++i){
        // auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        arma_x(i)     = x[u_size*i ];
        arma_x(i+n)   = x[u_size*i + 1];
        arma_x(i+n*2) = x[u_size*i + 2];
        arma_x(i+n*3) = x[u_size*i + 3];
    }
   
    // Eigen::VectorXd a2 = drbf->local_vipss_.final_h_eigen_ * arma_x;
    Eigen::VectorXd a2 = (arma_x.transpose() * drbf->local_vipss_.final_h_eigen_).transpose();
    if (!grad.empty()) {
        // std::cout << " grad size " <<grad.size() << std::endl;
        grad.resize(n*u_size);
        for(size_t i=0;i<n;++i){
            grad[i*u_size]   = a2(i);
            grad[i*u_size+1] = a2(i + n);
            grad[i*u_size+2] = a2(i + 2*n);
            grad[i*u_size+3] = a2(i + 3*n);
        }
    }
    double re = arma_x.dot(a2);
    double alpha = 1000.0;
    arma::vec re_vec(n);
    size_t id;
// #pragma omp parallel for shared(x, alpha, grad, re_vec) private(id)
    for(id =0; id < n; ++id)
    {
        double cur_re = x[u_size*id + 1] * x[u_size*id + 1] + x[u_size*id + 2] * x[u_size*id + 2] 
                        + x[u_size*id + 3] * x[u_size*id + 3] - 1;
        if(!grad.empty()) 
        {
            grad[u_size* id + 1] += alpha * 2 * x[u_size*id + 1] * cur_re; 
            grad[u_size* id + 2] += alpha * 2 * x[u_size*id + 2] * cur_re; 
            grad[u_size* id + 3] += alpha * 2 * x[u_size*id + 3] * cur_re; 
        }
        // re += alpha* cur_re * cur_re;
        re_vec[id] = alpha* cur_re * cur_re;
    }
    re += arma::accu(re_vec);
    VIPSSUnit::opt_func_count ++;

    // printf("res val : %f \n", re);
    return re;
}
// double myconstraint(unsigned n, const double *x, double *grad, void *data)

double contraint_simple(const std::vector<double>&x, std::vector<double>&grad, void *fdata)
{
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t id = drbf->constraint_count_ % (x.size()/3);
    // std::cout << " grad size " <<grad.size() << std::endl;
    // drbf->constraint_count_ ++;
    size_t npt = x.size() /3;
    // if(grad.empty())  grad.resize(x.size());
    double re =0;
    for(size_t id =0; id < npt; ++id)
    {
        double cur_re = x[3*id] * x[3*id] + x[3*id + 1] * x[3*id + 1] + x[3*id + 2] * x[3*id + 2] - 1;
        if(!grad.empty()) 
        {
           
            grad[3* id] = 2 * x[3*id] * cur_re; 
            grad[3* id + 1] = 2 * x[3*id + 1] * cur_re; 
            grad[3* id + 2] = 2 * x[3*id + 2] * cur_re; 
        }
        re +=  cur_re * cur_re;
    }
    // double re = x[0] * x[0] -1;
    return re;
}


double optfunc_unit_vipss(const std::vector<double>&x, std::vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    VIPSSUnit *drbf = reinterpret_cast<VIPSSUnit*>(fdata);
    // size_t n = drbf->npt;
    size_t n = drbf->npt_;
    // printf("input point number : %llu \n", n);
    arma::vec arma_x(n*4);

    //(  sin(a)cos(b), sin(a)sin(b), cos(a)  )  a =>[0, pi], b => [-pi, pi];
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
    arma::vec a2;
    a2 = drbf->Final_H_ * arma_x;
    if (!grad.empty()) {
        grad.resize(n*3);
        for(size_t i=0;i<n;++i){
            auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
            grad[i*3] = 2 * a2[i];
            grad[i*3 + 1] = 2 * (a2(n + i) * p_scsc[1] * p_scsc[3] + a2(i+ n * 2) * p_scsc[1] * p_scsc[2] - a2(i+n*3) * p_scsc[0]);
            grad[i*3 + 2] = 2* (-a2(i + n) * p_scsc[0] * p_scsc[2] + a2(i+n* 2) * p_scsc[0] * p_scsc[3]);
        }
    }
    double re = arma::dot( arma_x, a2 );
    // printf("residual val : %f \n", re);
    // printf("Final_H_ non zero  : %d \n", drbf->Final_H_.n_nonzero);
    // countopt++;
    // acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);
    drbf->countopt_ ++;
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
    printf("finish init solver ! \n");
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
        printf("start the solver ! \n");
        Solver::nloptwrapper(lower,upper,optfunc_unit_vipss_simple,this,1e-7,3000,solver_);
        // callfunc_time = acc_time;
        // solve_time = sol.time;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;
    }
    newnormals_.resize(npt_*3);
    printf("-----------newnormal size : %lu \n", newnormals_.size());
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
                this,1e-7,3000,solver_);
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

    Solver::nloptwrapperDirect(lower,upper,optfunc_unit_vipss_direct_eigen,
                this, 1e-7, 3000, solver_);
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
    printf("finish init solver ! \n");
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
        printf("start the solver ! \n");
        Solver::nloptwrapper(lower,upper,optfunc_unit_vipss,this,1e-7,3000,solver_);
    }
    newnormals_.resize(npt_*3);
    arma::vec y(npt_ + 3 * npt_);
    for(size_t i=0;i<npt_;++i)y(i) = 0;
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

void VIPSSUnit::BuildLocalHRBFPerNode()
{
    // std::vector<double> cluster_pts;
    // for(auto pt : local_vipss_.points_)
    // {
    //     size_t id = 
    // }

    // for(size_t i =0; i < cluster_num; ++i)
    // {
    //     // std::vector<double> vts;
    //     // std::vector<size_t> p_ids;
    //     auto cluster_pt_ids = GetClusterPtIds(i);
    //     auto cluster_pt_vec = GetClusterVerticesFromIds(cluster_pt_ids);
    //     auto t3 = Clock::now();
    //     vipss_api_.build_unit_vipss(cluster_pt_vec);
    //     time_sum += vipss_api_.rbf_core_.bigM_inv_time;
    // }
}



void VIPSSUnit::ReconSurface()
{
    printf(" start ReconSurface \n");
    // local_vipss_.TestVoronoiPts();
    // local_vipss_.voro_gen_.BuildTetMeshTetCenterMap();
    auto t000 = Clock::now();
    local_vipss_.voro_gen_.GenerateVoroData();
    auto t001 = Clock::now();
    double generate_voroi_data_time = std::chrono::nanoseconds(t001 - t000).count() / 1e9;
    printf("generate_voroi_data_time  : %f ! \n", generate_voroi_data_time);
    G_VP_stats.generate_voroi_data_time_ += generate_voroi_data_time;

    // local_vipss_.PtPCA(local_vipss_.out_pts_);
    // local_vipss_.out_normals_ = newnormals_;
    // local_vipss_.OptimizeAdjacentMat();
    // TestSpectralClustering(local_vipss_.cluster_adjacent_mat_opt_);
    // local_vipss_.voro_gen_.BuildTetMeshTetCenterMap();
    // local_vipss_.voro_gen_.BuildPicoTree();
    // bool use_nn_interpolation = true;
    auto t002 = Clock::now();
    bool is_group_cluster = false;
    local_vipss_.is_group_cluster_ = is_group_cluster;
    if(LOCAL_HRBF_NN == hrbf_type_)
    {
        printf(" start use_nn_interpolation \n");
        local_vipss_.normals_ = newnormals_;
        local_vipss_.s_vals_ = s_func_vals_;
        local_vipss_.user_lambda_ = user_lambda_;
        if(is_group_cluster)
        {
            local_vipss_.GroupClustersWithDegree();
            std::string group_pt_path = data_dir_  + file_name_ + "/" + "group_pts.obj";
            std::cout << " save group pts to file : " << group_pt_path << std::endl;
            local_vipss_.SaveGroupPtsWithColor(group_pt_path);
            local_vipss_.BuildHRBFPerCluster();
        } else {
            local_vipss_.BuildHRBFPerNode();
        }    
        auto t003 = Clock::now();
        double build_HRBF_time = std::chrono::nanoseconds(t003 - t002).count() / 1e9;
        G_VP_stats.build_per_cluster_hrbf_total_time_ += build_HRBF_time;     
        // printf(" start set local vipss static ptr \n");
        local_vipss_.SetThis();
        size_t n_voxels_1d = volume_dim_;
        Surfacer sf;
        auto surf_time = sf.Surfacing_Implicit(local_vipss_.out_pts_, n_voxels_1d, false, LocalVipss::NNDistFunction);
        auto t004 = Clock::now();
        sf.WriteSurface(finalMesh_v_,finalMesh_fv_);
        std::string out_path = data_dir_ + "/" + file_name_  + "/nn_surface";
        writePLYFile_VF(out_path, finalMesh_v_, finalMesh_fv_);
        auto t005 = Clock::now();
        double surface_file_save_time = std::chrono::nanoseconds(t005 - t004).count() / 1e9;

        double total_surface_time = std::chrono::nanoseconds(t005 - t000).count() / 1e9;


        printf("------- nn search time: %f \n", local_vipss_.search_nn_time_sum_);
        printf("------- cal nn coordinate and hrbf time: %f \n", local_vipss_.pass_time_sum_);
        printf(" ------ voxel dist func val evaluated count : %d  \n", LocalVipss::DistCallNum);
        printf(" ------ voxel dist func val evaluated time : %f \n", LocalVipss::DistCallTime);
        printf("------- total surface time: %f \n", total_surface_time);
        G_VP_stats.neighbor_search_time_ += local_vipss_.search_nn_time_sum_;
        G_VP_stats.cal_nn_coordinate_and_hbrf_time_ += local_vipss_.pass_time_sum_;
        G_VP_stats.voxel_cal_num += LocalVipss::DistCallNum;
        G_VP_stats.surface_total_time_ += total_surface_time;

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

void VIPSSUnit::Run()
{
    local_vipss_.out_dir_ = data_dir_ + "/" + file_name_ + "/";
    auto t00 = Clock::now();
    rbf_api_.Set_RBF_PARA();
    InitPtNormalWithLocalVipss();
    auto tn0 = Clock::now();
if(1)
{
    if(user_lambda_ < 1e-10)
    {
         local_vipss_.final_h_eigen_ = local_vipss_.final_h_eigen_.block(npt_, npt_, 3 * npt_, 3 * npt_);
    } 
    auto tn1 = Clock::now();
    double build_H_time = std::chrono::nanoseconds(tn1 - tn0).count() / 1e9;
    printf("Build H matrix time : %f ! \n", build_H_time);

    auto ts0 = Clock::now();
    if(user_lambda_ < 1e-12)
    {
        // OptUnitVipssNormalSimple();
        OptUnitVipssNormalDirectSimple();
    } else {
        // OptUnitVipssNormal();
        OptUnitVipssNormalDirect();
    }

    Final_H_.clear();
    local_vipss_.final_H_.clear();

    auto ts1 = Clock::now();
    double solve_time = std::chrono::nanoseconds(ts1 - ts0).count() / 1e9;
    printf("opt solve time : %f ! \n", solve_time);
    printf("opt fun call count : %d \n", VIPSSUnit::opt_func_count);

    G_VP_stats.opt_solver_time_ += solve_time;
    G_VP_stats.opt_func_call_num_ += VIPSSUnit::opt_func_count;
}

    auto t01 = Clock::now();
    double total_time = std::chrono::nanoseconds(t01 - t00).count()/1e9;
    printf("total local vipss running time : %f ! \n", total_time);

    std::string out_path  = local_vipss_.out_dir_ + local_vipss_.filename_  + "_opt";
    printf("out size : %lu, %lu \n", local_vipss_.out_pts_.size(), newnormals_.size());
    writePLYFile_VN(out_path, local_vipss_.out_pts_, newnormals_);

    std::string color_out_path  = local_vipss_.out_dir_ + local_vipss_.filename_  + "_opt_color";
    output_opt_pts_with_color(local_vipss_.out_pts_, s_func_vals_,color_out_path);

    // printf("start to ReconSurface 000 \n");
    // std::string init_path  = local_vipss_.out_dir_ + local_vipss_.filename_  + "_init";
    // writePLYFile_VN(init_path, local_vipss_.out_pts_, initnormals_);
    // printf("start to ReconSurface 0001 \n");
    ReconSurface();
    std::string log_path = local_vipss_.out_dir_ + "stats.txt";
    WriteStatsLog(log_path, G_VP_stats);
}

