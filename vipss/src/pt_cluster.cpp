#include "pt_cluster.h"
#include <cmath>
#include "readers.h"

namespace cluster
{
    void PtCluster::UpdateAdjacentClusters()
    {
        std::set<PtCluster*> delete_clusters;
        for(auto &cluster : adjacent_clusters_)
        {
            if(cluster->merged_)
            {
                delete_clusters.insert(cluster);
            } 
        }
        for(auto &d_cluster : delete_clusters)
        {
            adjacent_clusters_.erase(d_cluster);
            if(d_cluster->merged_cluster_ != this)
            {
                adjacent_clusters_.insert(d_cluster->merged_cluster_);
            }
        }
    }

    void PtCluster::ReplaceAdjacentClustersWithMergedClusters(PtCluster* merged_cluster)
    {
        for(auto& adj_cluster : adjacent_clusters_)
        {
            if(adj_cluster->adjacent_clusters_.find(this) != adj_cluster->adjacent_clusters_.end())
            {
                adj_cluster->adjacent_clusters_.erase(this);
            }
            if(adj_cluster != merged_cluster)
            adj_cluster->adjacent_clusters_.insert(merged_cluster);
        }
    }

    void PtCluster::MergeCluster(PtCluster* cluster)
    {
        key_pts_.insert(cluster->key_pts_.begin(), cluster->key_pts_.end());
        group_pts_.insert(cluster->group_pts_.begin(), cluster->group_pts_.end());

        for(auto &pt : key_pts_)
        {
            if(group_pts_.find(pt) != group_pts_.end())
            {
                group_pts_.erase(pt);
            }
        }
        adjacent_clusters_.insert(cluster->adjacent_clusters_.begin(), cluster->adjacent_clusters_.end());
        if(adjacent_clusters_.find(this) != adjacent_clusters_.end())
        {
            adjacent_clusters_.erase(this);
        }

        if(adjacent_clusters_.find(cluster) != adjacent_clusters_.end())
        {
            adjacent_clusters_.erase(cluster);
        }

        for(auto& cluster : cluster->adjacent_clusters_)
        {
            cluster->ReplaceAdjacentClustersWithMergedClusters(this);
        }

        this->CalcualteNormals();
        
        cluster->merged_ = true;
        cluster->merged_cluster_ = this;
        need_update_ = true;

    }

    std::vector<P3d_Ptr> PtCluster::GetClusterPts()
    {
        std::set<P3d_Ptr> pts;
        for(auto &pt : key_pts_)
        {
            pts.insert(pt);
        }

        for(auto &pt : group_pts_)
        {
            pts.insert(pt);
        }
        std::vector<P3d_Ptr> pts_vec(pts.begin(), pts.end()); 
        return pts_vec;
    }
    
    std::vector<vec3d> PtCluster::EstimateNormals(RBF_API& rbf_api)
    {
        std::vector<P3d_Ptr> pts = GetClusterPts();
        std::vector<double> pts_vec(pts.size() *3);
        for(size_t i = 0; i < pts.size(); ++i)
        {
            pts_vec[3*i] = pts[i][0];
            pts_vec[3*i + 1] = pts[i][0 + 2];
            pts_vec[3*i + 2] = pts[i][0 + 2];
        }
        rbf_api.run_vipss(pts_vec);

        pt_normal_map_.clear();
        std::vector<vec3d> out_normals;
        for(size_t i = 0; i < pts.size(); ++i )
        {
            rbf_api.normals_[3*i];
            vec3d n_normal{rbf_api.normals_[3*i], rbf_api.normals_[3*i + 1], rbf_api.normals_[3*i + 2]};
            pt_normal_map_[pts[i]] = n_normal;
            out_normals.push_back(n_normal);
        }
        return out_normals;
    }

    double M_PI2  = 2*acos(0.0);
    double CalculateScores(std::vector<arma::vec3>& a_normals, std::vector<arma::vec3>& b_normals)
    {
        double min_project_p = 1.0;
        for(size_t i = 0; i < a_normals.size(); ++i)
        {
            double ab_proj = arma::dot(a_normals[i], b_normals[i]);
            min_project_p = min_project_p < ab_proj ? min_project_p : ab_proj;
        }
        // printf("ab_proj val min positive : %f\n", min_project_p);
        double min_project_n = 1.0;
        for(size_t i = 0; i < a_normals.size(); ++i)
        {
            double ab_proj = - 1.0 * arma::dot(a_normals[i],  b_normals[i]);
            min_project_n = min_project_n < ab_proj ? min_project_n : ab_proj;
        }
        // printf("ab_proj val min negative : %f\n", min_project_n);

        min_project_p =  std::max(min_project_n, min_project_p);
        // min_project_p = std::min(1.0, min_project_p);

        // printf("---------max  pro : %f\n", min_project_p);

        double angle = acos (min_project_p) * 180.0 / M_PI2 ;
        // printf("---------max angle : %f\n", angle);
        return angle;
    }


    double CalculateClusterPairScore(PtCluster& a, PtCluster &b)
    {
        std::vector<vec3d> a_normals;
        std::vector<vec3d> b_normals;

        for(auto& pt : a.key_pts_)
        {
            if(b.pt_normal_map_.find(pt) != b.pt_normal_map_.end())
            {
                a_normals.push_back(a.pt_normal_map_[pt]);
                b_normals.push_back(b.pt_normal_map_[pt]);
            }
        }

        for(auto& pt : b.key_pts_ )
        {
            if(a.pt_normal_map_.find(pt) != a.pt_normal_map_.end())
            {
                a_normals.push_back(a.pt_normal_map_[pt]);
                b_normals.push_back(b.pt_normal_map_[pt]);
            }
        }

        double score = CalculateScores(a_normals, b_normals);
        return score;
    }

    void PtCluster::CalcualteNormals()
    {
        if(!vipss_api_)
        {
            printf(" vipss_api_ has not been initiliazed! \n");
            return;
        }
        std::vector<P3d_Ptr> pts = GetClusterPts();
        std::vector<double> pts_vec(pts.size() *3);
        for(size_t i = 0; i < pts.size(); ++i)
        {
            pts_vec[3*i] = pts[i][0];
            pts_vec[3*i + 1] = pts[i][0 + 1];
            pts_vec[3*i + 2] = pts[i][0 + 2];
        }
        vipss_api_->run_vipss(pts_vec);
        pt_normal_map_.clear();
        for(size_t i = 0; i < pts.size(); ++i )
        {
            vec3d n_normal{vipss_api_->normals_[3*i], vipss_api_->normals_[3*i + 1], vipss_api_->normals_[3*i + 2]};
            pt_normal_map_[pts[i]] = n_normal;
        }
    }


    void PtCluster::CalculateScoresWithAdjacentClusters()
    {
        adj_cluster_scores_.clear();
        
        for(auto &cluster : adjacent_clusters_)
        {
            if(cluster == this) continue;
            if(cluster->scored_times_ > scored_times_) continue;
            double score = CalculateClusterPairScore(*this, *cluster);
            adj_cluster_scores_[cluster] = score;
            if(cluster->adj_cluster_scores_.find(this) != cluster->adj_cluster_scores_.end())
            {
                cluster->adj_cluster_scores_[this] = score;
            }
        }
        double score_sum = 0.0;
        for(auto& ele : adj_cluster_scores_)
        {
            score_sum += ele.second;
        }
        score_ = 0.0;
        if(!adj_cluster_scores_.empty())
        {
            score_ = score_sum / double(adj_cluster_scores_.size());
        }
        scored_times_ ++;
    }

    void PtCluster::OutputPtWithColor(const std::string& out_path)
    {
        std::vector<double> pts;
        std::vector<uint8_t> colors;
        for(auto &pt : key_pts_)
        {
            pts.push_back(pt[0]);
            pts.push_back(pt[1]);
            pts.push_back(pt[2]);

            colors.push_back(255);
            colors.push_back(0);
            colors.push_back(0);
        }

        for(auto &pt : group_pts_)
        {
            pts.push_back(pt[0]);
            pts.push_back(pt[1]);
            pts.push_back(pt[2]);

            colors.push_back(0);
            colors.push_back(0);
            colors.push_back(255);
        }
        writePLYFile_CO(out_path, pts, colors);
    }
}