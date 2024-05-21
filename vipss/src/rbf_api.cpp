#include "rbf_api.h"

void RBF_API::Set_RBF_PARA(){

    // std::cout << " Start to set rbf para ..... " << endl;
    RBF_Paras &para= para_;
    RBF_InitMethod initmethod = Lamnbda_Search;

    RBF_Kernal Kernal = XCube;
    int polyDeg = 1;
    double sigma = 0.9;
    double rangevalue = 0.001;

    para.Kernal = Kernal;para.polyDeg = polyDeg;
    para.sigma = sigma;para.rangevalue = rangevalue;
    para.Hermite_weight_smoothness = 0.0;
    para.Hermite_ls_weight = 0;
    para.Hermite_designcurve_weight = 00.0;
    para.Method = RBF_METHOD::Hermite_UnitNormal;

    para.InitMethod = initmethod;
    para.user_lamnbda = 0;
    para.isusesparse = false;
    // return para;
}


void RBF_API::run_vipss(std::vector<double> &Vs)
{
    // vector<double>Vs;
    // std::cout << "start to set rbf para "<< std::endl;
    // RBF_Paras para = this->Set_RBF_PARA();
    // std::cout << Vs.size() << endl;
    para_.user_lamnbda = user_lambda_;
    // std::cout << "start inject data "<< std::endl;
    RBF_Core rbf_core_;
    rbf_core_.InjectData(Vs, para_);
    // std::cout << "finish inject data "<< std::endl;
    rbf_core_.BuildK(para_);
    // std::cout << "finish build K "<< std::endl;
    rbf_core_.InitNormal(para_);
    rbf_core_.OptNormal(0);
    
    
    if(is_surfacing_){
        rbf_core_.Write_Hermite_NormalPrediction(outpath_ + "_normal", 1);
        rbf_core_.Surfacing(0,n_voxel_line_);
        rbf_core_.Write_Surface(outpath_ +"_surface");
    }
    if(is_outputtime_){
        rbf_core_.Print_TimerRecord_Single(outpath_ +"_time.txt");
    }
    normals_ = rbf_core_.newnormals;
    // rbf_core_.clear();
}