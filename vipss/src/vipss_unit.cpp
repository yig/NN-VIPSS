#include <chrono>
#include "vipss_unit.hpp"

typedef std::chrono::high_resolution_clock Clock;


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
    
    local_vipss_.BuildMatrixH();
    npt_ = local_vipss_.points_.size();
    initnormals_ = local_vipss_.out_normals_;

    printf("unit vipss J mat init time : %f \n", local_vipss_.vipss_api_.u_v_time);
}


void VIPSSUnit::BuildVipssUnitMatrixP()
{
    auto& cluster_pt_mat_vec = local_vipss_.cluster_pt_mat_vec_;
    auto& cluster_J_mat_vec  = local_vipss_.cluster_J_mat_vec_;
    // auto& cluster_pt_ids_vec = local_vipss_.cluster_pt_ids_vec_;

    auto& pts = local_vipss_.points_;
    npt_ = pts.size();
    Final_H_.resize(4 * npt_, 4 * npt_);
    temp_Hs_.resize(npt_);
// #pragma omp parallel for
    for(size_t i = 0; i < cluster_pt_mat_vec.size(); ++i)
    {
        const auto& unit_cluster_mat = cluster_pt_mat_vec[i];
        const auto& J_m = cluster_J_mat_vec[i];
        // const auto& current_ids = cluster_pt_ids_vec[i];
        if(user_lambda_ < 1e-10)
        {            
            arma::sp_mat sp_J(J_m); 
            arma::sp_mat temp_H = unit_cluster_mat.t() * sp_J * unit_cluster_mat;
            temp_Hs_[i] = temp_H;
        } else {
            size_t unit_npt = unit_cluster_mat.n_rows / 4; 
            arma::mat F(4 * unit_npt, 4 * unit_npt);
            arma::mat E(unit_npt, unit_npt);
            E.eye();
            F(0, 0, arma::size(unit_npt, unit_npt)) = E;
            arma::mat K = (F + J_m * (user_lambda_));
            arma::sp_mat sp_K(K);  
            arma::sp_mat temp_H = unit_cluster_mat.t() * sp_K * unit_cluster_mat;
            temp_Hs_[i] = temp_H;
        }
    }
    size_t cur_n = npt_;
    size_t cur_loop = 0;
    printf("npt_ %lu log 2 : %d \n", npt_, (int)log2(npt_));
    size_t loop_level = (size_t)log2(npt_);
    if(npt_ > pow(2, loop_level)) loop_level ++;
    for(size_t cur_loop = 1; cur_loop <= loop_level; ++cur_loop)
    {
        size_t cur_step = pow(2, cur_loop);
        // printf("-----current N : %llu \n", cur_n);
        // size_t remain =  cur_n % 2;
        if(cur_n / 2 == 0) break;
        cur_n = cur_n / 2 + cur_n % 2; 
        // printf("current N 0: %llu \n", cur_n);
    // #pragma omp parallel for
        for(size_t m_i = 0; m_i < cur_n; ++m_i)
        {
            if(m_i*cur_step + cur_step/2 < npt_)
            {
                temp_Hs_[m_i*cur_step] += temp_Hs_[m_i*cur_step + cur_step/2]; 
            }
        } 
    }
    Final_H_ = temp_Hs_[0];
    // auto& voro_gen = local_vipss_.voro_gen_;
    // auto& pt_ids = local_vipss_.cluster_core_pt_ids_;
    // auto& cluster_pt_id_map = voro_gen.point_id_map_;
     
    // Final_H_.resize(4 * npt_, 4 * npt_);
    // for(auto& pt : pts)
    // {
    //     const auto& cluster_pts = voro_gen.point_cluster_pts_map_[pt];
    //     std::vector<size_t> unit_pt_ids;
    //     std::vector<double> pts_data;
    //     size_t pt_id = cluster_pt_id_map[pt];
    //     local_vipss_.GetInitClusterPtIds(pt_id, pts_data, unit_pt_ids);
    //     size_t unit_npt = unit_pt_ids.size();
    //     arma::sp_mat unit_cluster_mat(unit_npt*4, npt_ * 4);
    //     for(size_t id = 0; id < unit_pt_ids.size(); ++id)
    //     {
    //         size_t pid = unit_pt_ids[id];
    //         unit_cluster_mat(id, pid) = 1.0;
    //         unit_cluster_mat(id + unit_npt,     pid + npt_ ) = 1.0;
    //         unit_cluster_mat(id + unit_npt * 2, pid + npt_ * 2) = 1.0;
    //         unit_cluster_mat(id + unit_npt * 3, pid + npt_ * 3) = 1.0;
    //     }
    //     rbf_api_.build_unit_vipss(pts_data);
    //     arma::sp_mat J_m(rbf_api_.rbf_core_.Minv);
    //     if(lambda_ < 1e-10)
    //     {
    //         Final_H_ += unit_cluster_mat.t() * J_m * unit_cluster_mat;
    //     } else {
    //         arma::sp_mat F(4 * unit_npt, 4 * unit_npt);
    //         arma::sp_mat E(unit_npt, unit_npt);
    //         E.eye();
    //         F(0, 0, arma::size(unit_npt, unit_npt)) = E;
    //         Final_H_ += unit_cluster_mat.t() * (F + J_m * lambda_) * unit_cluster_mat;
    //     }
    // }
    if(user_lambda_ < 1e-10)
    {
        // printf("Final_H_ rows : %llu, cols : %llu \n", Final_H_.n_rows, Final_H_.n_cols);
        // auto sub_H_ = Final_H_(0, 0, arma::size(npt_, npt_));
        // printf("sub_H_ non zero final : %d \n", sub_H_.n_nonzero);
        Final_H_ = Final_H_(npt_, npt_, arma::size(3 *npt_, 3 * npt_));
        // printf("Final_H_ rows : %llu, cols : %llu \n", Final_H_.n_rows, Final_H_.n_cols);
        // printf("Final_H_ non zero final : %d \n", Final_H_.n_nonzero);
    }    
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

    Solver::nloptwrapperDirect(lower,upper,optfunc_unit_vipss_direct_simple,
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

    Solver::nloptwrapperDirect(lower,upper,optfunc_unit_vipss_direct,
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
    local_vipss_.voro_gen_.GenerateVoroData();
    // local_vipss_.voro_gen_.BuildTetMeshTetCenterMap();
    // local_vipss_.voro_gen_.BuildPicoTree();
    // return;
    
    // bool use_nn_interpolation = true;
    if(LOCAL_HRBF_NN == hrbf_type_)
    {
        printf(" start use_nn_interpolation \n");
        local_vipss_.normals_ = newnormals_;
        local_vipss_.s_vals_ = s_func_vals_;
        local_vipss_.user_lambda_ = user_lambda_;
        local_vipss_.BuildHRBFPerNode(); 
        // printf(" start set local vipss static ptr \n");
        local_vipss_.SetThis();
        // printf(" finish set local vipss static ptr \n");
        // local_vipss_.testNNPtDist();

    if(0)
    {
        // double p0[3] = {0.266693, 0.369411, 0.0690456};
        // double p1[3] = {0.278074, 0.326238, 0.029064};
        // double p2[3] = {0.373304, 0.289721, 0.0855713};
        double p0[3] = {-0.142955, 0.147453, -0.273193};
        double p1[3] = {-0.145742, 0.241941, -0.243457};
        double p2[3] = {-0.146177, 0.201896, -0.202977};
        VoroPlane visual_plane(&p0[0], &p1[0], &p2[0]);
        // std::string plane_save_path = data_dir_ + file_name_ + "/visual_func_plane.obj";
        // visual_plane.SavePlane(plane_save_path);

        std::string visual_func_path = data_dir_ + file_name_ + "/visual_func_vals.obj";
        local_vipss_.VisualFuncValues(LocalVipss::NNDistFunction, visual_plane, visual_func_path);
    }
    
        size_t n_voxels_1d = 100;
        Surfacer sf;
        auto surf_time = sf.Surfacing_Implicit(local_vipss_.out_pts_, n_voxels_1d, false, LocalVipss::NNDistFunction);
        sf.WriteSurface(finalMesh_v_,finalMesh_fv_);
        printf(" ------ DistCallNum : %d  \n", LocalVipss::DistCallNum);
        printf(" ------ DistCallTime : %f \n", LocalVipss::DistCallTime);
        std::string out_path = data_dir_ + "/" + file_name_  + "/nn_surface";
        writePLYFile_VF(out_path, finalMesh_v_, finalMesh_fv_);
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
        Final_H_ = local_vipss_.final_H_(npt_, npt_, arma::size(3 *npt_, 3 * npt_));
    } else {
        Final_H_ = local_vipss_.final_H_;
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

    auto ts1 = Clock::now();
    double solve_time = std::chrono::nanoseconds(ts1 - ts0).count() / 1e9;
    printf("opt solve time : %f ! \n", solve_time);
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
}
