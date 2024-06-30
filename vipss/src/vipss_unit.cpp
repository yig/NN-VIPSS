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
    local_vipss_.user_lambda_ = lambda_;
    // local_vipss_.max_iter_ = 30;
    local_vipss_.use_hrbf_surface_ = false;
    // local_vipss_.volume_dim_ = 100;
    local_vipss_.Init(path);
    local_vipss_.InitNormals();
    initnormals_ = local_vipss_.out_normals_;
}


void VIPSSUnit::BuildVipssUnitMatrixP()
{
    auto& pts = local_vipss_.points_;
    auto& voro_gen = local_vipss_.voro_gen_;
    auto& pt_ids = local_vipss_.cluster_core_pt_ids_;
    auto& cluster_pt_id_map = voro_gen.point_id_map_;
    npt_ = pts.size(); 
    
    Final_H_.resize(4 * npt_, 4 * npt_);
    
    for(auto& pt : pts)
    {
        const auto& cluster_pts = voro_gen.point_cluster_pts_map_[pt];
        std::vector<size_t> unit_pt_ids;
        std::vector<double> pts_data;
        size_t pt_id = cluster_pt_id_map[pt];
        local_vipss_.GetInitClusterPtIds(pt_id,pts_data,  unit_pt_ids);
        
        size_t unit_npt = unit_pt_ids.size();
        arma::sp_mat unit_cluster_mat(unit_npt*4, npt_ * 4);
        // printf("----------unit_cluster_mat col size %llu \n", unit_cluster_mat.n_cols);

        for(size_t id = 0; id < unit_pt_ids.size(); ++id)
        {

            size_t pid = unit_pt_ids[id];
            // printf("cur p id : %d \n", pid);
            unit_cluster_mat(id, pid) = 1.0;
            unit_cluster_mat(id + unit_npt,     pid + npt_ ) = 1.0;
            unit_cluster_mat(id + unit_npt * 2, pid + npt_ * 2) = 1.0;
            unit_cluster_mat(id + unit_npt * 3, pid + npt_ * 3) = 1.0;
            // printf("cur pid + npt_ * 3 : %d \n", pid + npt_ * 3);
            // break;
        }

        // arma::sp_mat::const_row_col_iterator it     = unit_cluster_mat.begin_row_col();
        // arma::sp_mat::const_row_col_iterator it_end = unit_cluster_mat.end_row_col();

        // for (; it != it_end; ++it) {
        //     cout << "-------val: " << (*it)    << endl;
        //     cout << "row: " << it.row() << endl;
        //     cout << "col: " << it.col() << endl;
        // }
        // unit_matrix_vec_.push_back(unit_cluster_mat);   
        rbf_api_.build_unit_vipss(pts_data);
        // arma::mat minv = rbf_api_.rbf_core_.Minv;
        arma::sp_mat J_m(rbf_api_.rbf_core_.Minv);
        // printf("minv rows %d, cols %d", minv.n_rows, minv.n_cols);
        // printf("J_m rows %d, cols %d \n", J_m.n_rows, J_m.n_cols);
        if(lambda_ < 1e-10)
        {
            // printf("J_m non zero  : %d \n", J_m.n_nonzero);
            // printf("unit_npt : %d \n", unit_npt);
            Final_H_ += unit_cluster_mat.t() * J_m * unit_cluster_mat;
            // auto sub_H_ = Final_H_(0, 0, arma::size(npt_, npt_));
            // printf("sub_H_ non zero  : %d \n", sub_H_.n_nonzero);
            // arma::sp_mat temp_m = unit_cluster_mat.t()  * unit_cluster_mat;
            // printf("----temp_m non zero  : %d \n", temp_m.n_nonzero);
            // auto sub_m_ = temp_m(npt_, npt_, arma::size(3 * npt_, 3 * npt_));
            // printf("----sub_m_ non zero  : %d \n", sub_m_.n_nonzero);
        } else {
            arma::sp_mat F(4 * unit_npt, 4 * unit_npt);
            arma::sp_mat E(unit_npt, unit_npt);
            E.eye();
            F(0, 0, arma::size(unit_npt, unit_npt)) = E;
            Final_H_ += unit_cluster_mat.t() * (F + J_m * lambda_) * unit_cluster_mat;
        }
        
    }
    if(lambda_ < 1e-10)
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
        arma_x(i) = p_scsc[0] * p_scsc[3];
        arma_x(i+n) = p_scsc[0] * p_scsc[2];
        arma_x(i+n*2) = p_scsc[1];
    }
    arma::vec a2;
    a2 = drbf->Final_H_ * arma_x;
    if (!grad.empty()) {
        grad.resize(n*2);
        for(size_t i=0;i<n;++i){
            auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
            grad[i*2] = a2(i) * p_scsc[1] * p_scsc[3] + a2(i+n) * p_scsc[1] * p_scsc[2] - a2(i+n*2) * p_scsc[0];
            grad[i*2+1] = -a2(i) * p_scsc[0] * p_scsc[2] + a2(i+n) * p_scsc[0] * p_scsc[3];
        }
    }
    double re = arma::dot( arma_x, a2 );
    // printf("residual val : %f \n", re);
    // printf("Final_H_ non zero  : %d \n", drbf->Final_H_.n_nonzero);
    // countopt++;
    // acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);
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
            grad[i*3] = a2[i];
            grad[i*3 + 1] = a2(n + i) * p_scsc[1] * p_scsc[3] + a2(i+ n * 2) * p_scsc[1] * p_scsc[2] - a2(i+n*3) * p_scsc[0];
            grad[i*3 + 2] = -a2(i + n) * p_scsc[0] * p_scsc[2] + a2(i+n* 2) * p_scsc[0] * p_scsc[3];
        }
    }
    double re = arma::dot( arma_x, a2 );
    // printf("residual val : %f \n", re);
    // printf("Final_H_ non zero  : %d \n", drbf->Final_H_.n_nonzero);
    // countopt++;
    // acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);
    return re;
}

void VIPSSUnit::OptUnitVipssNormalSimple(){

    printf("start to call solver ! \n");
    solver_.solveval.resize(npt_ * 2);

    for(size_t i=0;i<npt_;++i){
        double *veccc = initnormals_.data()+i*3;
        {
            solver_.solveval[i*2] = atan2(sqrt(veccc[0]*veccc[0]+veccc[1]*veccc[1]),veccc[2] );
            solver_.solveval[i*2 + 1] = atan2( veccc[1], veccc[0]   );
        }
    }
    printf("finish init solver ! \n");
    if(1){
        std::vector<double>upper(npt_*2);
        std::vector<double>lower(npt_*2);
        for(int i=0;i<npt_;++i){
            upper[i*2] = 2 * M_PI_;
            upper[i*2 + 1] = 2 * M_PI_;

            lower[i*2] = -2 * M_PI_;
            lower[i*2 + 1] = -2 * M_PI_;
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
    arma::vec y(npt_ + 3 * npt_);
    for(size_t i=0;i<npt_;++i)y(i) = 0;
    for(size_t i=0;i<npt_;++i){
        double a = solver_.solveval[i*2], b = solver_.solveval[i*2+1];
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
            upper[i*3 + 1] = 2 * M_PI_;
            upper[i*3 + 2] = 2 * M_PI_;

            lower[i*3] = -1.0;
            lower[i*3 + 1] = -2 * M_PI_;
            lower[i*3 + 2] = -2 * M_PI_;
        }
        // countopt = 0;
        // acc_time = 0;
        //LocalIterativeSolver(sol,kk==0?normals:newnormals,300,1e-7);
        printf("start the solver ! \n");
        Solver::nloptwrapper(lower,upper,optfunc_unit_vipss,this,1e-7,3000,solver_);
        // callfunc_time = acc_time;
        // solve_time = sol.time;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;
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

void VIPSSUnit::ReconSurface()
{
    if(use_hrbf_surface_)
    {
        rbf_api_.user_lambda_ = lambda_;
        rbf_api_.outpath_ = data_dir_ + file_name_ + "/";
        rbf_api_.is_surfacing_ = true;
        rbf_api_.run_vipss(local_vipss_.out_pts_, newnormals_);
    }
}

void VIPSSUnit::Run()
{
    rbf_api_.Set_RBF_PARA();
    InitPtNormalWithLocalVipss();
    BuildVipssUnitMatrixP();
    if(lambda_ < 1e-10)
    {
        OptUnitVipssNormalSimple();
    } else {
        OptUnitVipssNormal();
    }
    
    std::string out_path  = local_vipss_.out_dir_ + local_vipss_.filename_  + "_opt";
    writePLYFile_VN(out_path, local_vipss_.out_pts_, newnormals_);

    ReconSurface();
}