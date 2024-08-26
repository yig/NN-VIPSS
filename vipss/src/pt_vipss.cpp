#include "pt_vipss.h"
#include "readers.h"


void PtVipss::Init(const std::string& path)
{
    voro_gen_.loadData(path);
    voro_gen_.InitMesh();
    voro_gen_.BuildPtIdMap();
    pt_num_ = voro_gen_.pt_num_;
    pt_cluster_map_ = voro_gen_.BuildAllClusters();
    vipss_api_.Set_RBF_PARA();
    vipss_api_.is_surfacing_ = false;
    printf("finish init pt vipss and build pt cluster map ! \n");

}

void PtVipss::InitClusters()
{
    size_t c_count = 0;
    for(auto &ele: pt_cluster_map_)
    {
        pts_.push_back(ele.first);
        // std::string out_path = out_dir_ + out_name_ + "_cluster_" + std::to_string(c_count);
        // ele.second.OutputPtWithColor(out_path);
        ele.second.vipss_api_ = &vipss_api_;
        ele.second.CalcualteNormals();
        pt_clusters_.push_back(&ele.second);
        for(auto &pt : ele.second.group_pts_)
        {
            if(pt == ele.first) continue;
            if(pt_cluster_map_.find(pt) != pt_cluster_map_.end())
            {
                ele.second.adjacent_clusters_.insert(&pt_cluster_map_[pt]);
            }
        }
        c_count ++;
    }
    printf("finish init cluster normals and adjacent clusters ! \n");
    std::string out_path = out_dir_ + out_name_ + "_ptn0" ;
    UpdatePtsAndNormals();
    OutPutNormals(out_path);

    for(auto& cluster: pt_clusters_)
    {
        cluster->CalculateScoresWithAdjacentClusters();
    }

    // for(auto& cluster: pt_clusters_)
    // {
    //     printf("cluster score : %f \n", cluster->score_);
    // }

    printf("finish calculate adjacent clusters scores! \n");

    std::sort(pt_clusters_.begin(), pt_clusters_.end(), 
            [](PtCluster* a, PtCluster*b) {return a->score_ > b->score_;});
    
    MergeClusters();
    
}

void PtVipss::MergeClusters()
{
    std::set<PtCluster*> visited_clusters;
    merge_pairs_num_ = 0;
    for(auto& cluster: pt_clusters_)
    {
        if(cluster->score_ < angle_threshold_) continue;
        if(visited_clusters.find(cluster) != visited_clusters.end()) continue;
        visited_clusters.insert(cluster);
        // for(auto& ele : )
        std::vector<std::pair<PtCluster*, double>> adj_clusters_with_scores(cluster->adj_cluster_scores_.begin(), cluster->adj_cluster_scores_.end());

        std::sort(adj_clusters_with_scores.begin(), adj_clusters_with_scores.end(), 
            []( std::pair<PtCluster*, double>& a, std::pair<PtCluster*, double>& b)
        {
            return a.second > b.second;
        });
        PtCluster* to_merge_cluster = NULL; 
        for(auto &ele : adj_clusters_with_scores)
        {
            if(visited_clusters.find(ele.first) == visited_clusters.end())
            {
                to_merge_cluster = ele.first;
                break;
            }
        }
        if(to_merge_cluster)
        {
            cluster->MergeCluster(to_merge_cluster);
            visited_clusters.insert(to_merge_cluster);
            for(auto & pt : to_merge_cluster->key_pts_)
            {
                pt_cluster_map_[pt] = *cluster;
            }
        }
        merge_pairs_num_++;
    }
    printf("There are %d number of cluster pairs merged! \n", merge_pairs_num_);
}


void PtVipss::UpdateClusters()
{
    std::vector<PtCluster*> updated_clusters;
    for(auto& cluster: pt_clusters_)
    {
        if(cluster->merged_) continue;
        updated_clusters.push_back(cluster);
    }
    pt_clusters_ = updated_clusters;

    for(auto& cluster: pt_clusters_)
    {
        if(cluster->need_update_)
        {
            // cluster->InitNormals(vipss_api_);
            cluster->CalculateScoresWithAdjacentClusters();
            cluster->need_update_ = false;
        }
    }

    std::sort(pt_clusters_.begin(), pt_clusters_.end(), 
            [](PtCluster* a, PtCluster*b) {return a->score_ > b->score_;});

    MergeClusters();

}

void PtVipss::UpdatePtsAndNormals()
{
    
    normals_.clear();
    arma::vec3 init_n{0, 0, 0};
    for(auto& pt : pts_)
    {
        auto cluster = pt_cluster_map_[pt];
        if(cluster.pt_normal_map_.find(pt) == cluster.pt_normal_map_.end())
        {
            printf("current pt is not contained by a cluster \n");
            normals_.push_back(init_n);
        } else {
            auto cur_n = cluster.pt_normal_map_[pt];
            normals_.push_back(cur_n);
        }
    }
}

void PtVipss::OutPutNormals(const std::string& out_path)
{
    // std::string out_path = out_dir_ + out_name_ + "_ptn" ;
    std::vector<double> pts;
    std::vector<double> normals;
    for(auto &pt : pts_)
    {
        pts.push_back(pt[0]);
        pts.push_back(pt[1]);
        pts.push_back(pt[2]);
    }
    for(auto &nv : normals_ )
    {
        normals.push_back(nv[0]);
        normals.push_back(nv[1]);
        normals.push_back(nv[2]);
    }
    
    writePLYFile_VN(out_path, pts, normals);
}

void PtVipss::Run()
{
    InitClusters();
    printf("Finish init clusters! \n");
    std::string out_path = out_dir_ + out_name_ + "_ptn0" ;
    // UpdatePtsAndNormals();
    // OutPutNormals(out_path);

    printf("save path : %s", out_path.c_str());


    int iter = 0;
    while(merge_pairs_num_ != 0)
    {
        iter ++;
        MergeClusters();
        out_path = out_dir_ + out_name_ + "_ptn_" + std::to_string(iter);
        UpdatePtsAndNormals();
        OutPutNormals(out_path);
    }
    

}

