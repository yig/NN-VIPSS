#include "incremental_vipss.h"
#include <chrono>

#include "readers.h"

typedef std::chrono::high_resolution_clock Clock;

void IncreVipss::Init()
{
    auto t0 = Clock::now();
    readXYZ(input_path_,in_pts_);
    pt_sampler_.init(in_pts_);
    pt_sampler_.FurthestSamplePointCloud(init_sample_num_);
    pt_sampler_.SplitSamplePts(key_pts_, auxi_pts_);
    if(save_init_sample_pts_)
    {
        std::string key_vs_out_path = out_dir_ + std::to_string(init_sample_num_) + "_key_vs.xyz";
        writeXYZ(key_vs_out_path, key_pts_);    
    }

    sample_tree_.InitTree(key_pts_);
    cur_sample_dist_ = sqrt(pt_sampler_.cur_min_dist);

    auto t1 = Clock::now();
    double pre_time = (std::chrono::nanoseconds(t1 - t0).count()/1e9);
    printf("incre vipss init time %f \n", pre_time);

    InitRBFPara();
}

void IncreVipss::InitRBFPara()
{
    
    vipss_api_.Set_RBF_PARA();
    vipss_api_.is_surfacing_ = is_surfacing_;
    vipss_api_.user_lambda_ = lambda_;
    vipss_api_.n_voxel_line_ = n_voxel_line_;
    vipss_api_.outpath_ = out_dir_;
    vipss_api_.incre_vipss_beta_ = residual_beta_;
}

void IncreVipss::SampleIncrePts()
{
    if(vipss_api_.auxi_dist_vec_.size() > 0) 
    {
        arma::uvec indices = arma::sort_index(vipss_api_.auxi_dist_vec_, "descend");
        printf("indices size : %llu \n", indices.size());
        std::vector<double> new_auxi_vs;
        size_t add_count = 0;
        double add_dist = cur_sample_dist_ * add_pt_dist_scale_; 
        for(size_t id = 0; id < indices.size(); ++id)
        {
            size_t cur_id = indices(id);
            if(add_count < incremental_num_)
            {
                MKdtree::Point_d cur_p(auxi_pts_[cur_id*3], auxi_pts_[3* cur_id + 1], auxi_pts_[3*cur_id + 2]);
                if(sample_tree_.InsertPt(cur_p, add_dist))
                {
                    add_count ++;
                    continue;
                }
            }
            new_auxi_vs.push_back(auxi_pts_[cur_id*3]);
            new_auxi_vs.push_back(auxi_pts_[cur_id*3 + 1]);
            new_auxi_vs.push_back(auxi_pts_[cur_id*3 + 2]);
        }
        // printf("finish adding new key pts \n");
        auxi_pts_ = new_auxi_vs;
        new_auxi_vs.clear();
        add_pt_dist_scale_ *= pt_dist_scale_decay_;
        key_pts_ = sample_tree_.GetTreePts();
        key_init_normals_ = vipss_api_.rbf_core_.EstimateNormals(key_pts_);
        if(save_normals_)
        {
            std::string key_vs_out_path = out_dir_ + std::to_string(key_pts_.size()/3) + "_key_pt_normal.xyz";
            // writePLYFile_VN(key_vs_out_path, key_pts_, key_init_normals_);
        }

    }

}

void IncreVipss::Run()
{
    for(size_t i = 0; i < iter_num_; ++i)
    {
        size_t key_v_num = key_pts_.size() / 3;
        auto all_pts = key_pts_;

        all_pts.insert(all_pts.end(),   std::make_move_iterator(auxi_pts_.begin()),
                                        std::make_move_iterator(auxi_pts_.end()) );
        if(!key_init_normals_.empty())
        {
            vipss_api_.use_input_normal_ = true;
            vipss_api_.rbf_core_.initnormals = key_init_normals_;
        }
        vipss_api_.incre_vipss_beta_ = residual_beta_;
        vipss_api_.run_vipss_for_incremental(all_pts, key_v_num);
        printf("finish vipss optimization \n");
        SampleIncrePts();
    }    
}