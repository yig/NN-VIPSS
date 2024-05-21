#pragma once
#include "rbfcore.h"

class RBF_API{
    public:
        RBF_API(){};
        ~RBF_API(){};
        void Set_RBF_PARA();
        void run_vipss(std::vector<double> &Vs);


    public:
        RBF_Paras para_;
        bool is_surfacing_ = true;
        bool is_outputtime_ = false;
        double user_lambda_ = 0.0001; 
        int n_voxel_line_ = 30;
        std::vector<double> normals_;
        
        std::string outpath_ = "./vipss";



};