#include "SimpleOctree.h"
#include <numeric>
#include <iostream>
#include <fstream>

namespace SimOctree{


inline double SquareDistance(const Point& p1, const Point& p2)
{
    double dx = p1[0] - p2[0];
    double dy = p1[1] - p2[1];
    double dz = p1[2] - p2[2];
    return dx * dx + dy * dy + dz * dz;
}

OctBBox::OctBBox(const Point& minP, const Point& maxP) 
:min_pt_(minP), max_pt_(maxP)  
{

}

inline bool IsInsideBox(const Point& pt, const OctBBox& box)
{
    if(pt[0] >= box.min_pt_[0] && pt[0] <= box.max_pt_[0])
    {
        if(pt[1] >= box.min_pt_[1] && pt[1] <= box.max_pt_[1])
        {
            if(pt[2] >= box.min_pt_[2] && pt[2] <= box.max_pt_[2])
            {
                return true;
            }
        }
    }
    return false;
}

void OctBBox::SaveBoxMesh(const std::string& path) const
{
    std::ofstream file(path);
    auto corners = GenerateCorners();
    for(const auto& pt : corners)
    {
        file <<"v " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
    }
    file.close();
}

std::vector<Point> OctBBox::GenerateCorners() const
{

    std::vector<Point> corner_pts;
    corner_pts.push_back(min_pt_);
    corner_pts.push_back({min_pt_[0], max_pt_[1], min_pt_[2]});
    corner_pts.push_back({max_pt_[0], max_pt_[1], min_pt_[2]});
    corner_pts.push_back({max_pt_[0], min_pt_[1], min_pt_[2]});
    corner_pts.push_back({min_pt_[0], min_pt_[1], max_pt_[2]});
    corner_pts.push_back({min_pt_[0], max_pt_[1], max_pt_[2]});
    corner_pts.push_back({max_pt_[0], max_pt_[1], max_pt_[2]});
    corner_pts.push_back({max_pt_[0], min_pt_[1], max_pt_[2]});
    return corner_pts;
}

std::vector<OctBBox> OctBBox::DividBBox() const
{
    std::array<double,3> center;
    center[0] = (min_pt_[0] + max_pt_[0]) * 0.5;
    center[1] = (min_pt_[1] + max_pt_[1]) * 0.5;
    center[2] = (min_pt_[2] + max_pt_[2]) * 0.5;

    std::vector<OctBBox> child_boxes;

    //     5-------6
    //    /|      /|
    //   4-|-----7 |
    //   | 1-----|-2
    //   |/      |/
    //   0-------3 
    Point corners[8];
    corners[0] = min_pt_;
    corners[1] = {min_pt_[0], max_pt_[1], min_pt_[2]};
    corners[2] = {max_pt_[0], max_pt_[1], min_pt_[2]};
    corners[3] = {max_pt_[0], min_pt_[1], min_pt_[2]};
    corners[4] = {min_pt_[0], min_pt_[1], max_pt_[2]};
    corners[5] = {min_pt_[0], max_pt_[1], max_pt_[2]};
    corners[6] = {max_pt_[0], max_pt_[1], max_pt_[2]};
    corners[7] = {max_pt_[0], min_pt_[1], max_pt_[2]};
    
    for(int i = 0; i < 8; ++i)
    {
        OctBBox newBox;
        newBox.GetCornersFromPts(center, corners[i]);
        child_boxes.push_back(newBox);
    }

    return child_boxes;
}

void OctBBox::GetCornersFromPts(const std::vector<double>& pts)
{
    size_t ptn = pts.size() / 3;
    min_pt_[0] = std::numeric_limits<double>::max();
    min_pt_[1] = std::numeric_limits<double>::max();
    min_pt_[2] = std::numeric_limits<double>::max();

    max_pt_[0] = std::numeric_limits<double>::min();
    max_pt_[1] = std::numeric_limits<double>::min();
    max_pt_[2] = std::numeric_limits<double>::min();

    for(size_t i = 0; i < ptn; ++i)
    {
        min_pt_[0] = min_pt_[0] < pts[3*i]     ? min_pt_[0] : pts[3*i];
        min_pt_[1] = min_pt_[1] < pts[3*i + 1] ? min_pt_[1] : pts[3*i + 1];
        min_pt_[2] = min_pt_[2] < pts[3*i + 2] ? min_pt_[2] : pts[3*i + 2];

        max_pt_[0] = max_pt_[0] > pts[3*i]     ? max_pt_[0] : pts[3*i];
        max_pt_[1] = max_pt_[1] > pts[3*i + 1] ? max_pt_[1] : pts[3*i + 1];
        max_pt_[2] = max_pt_[2] > pts[3*i + 2] ? max_pt_[2] : pts[3*i + 2];
    }
}

void OctBBox::GetCornersFromPts(const std::vector<Point>& pts)
{
    size_t ptn = pts.size();
    min_pt_[0] = std::numeric_limits<double>::max();
    min_pt_[1] = std::numeric_limits<double>::max();
    min_pt_[2] = std::numeric_limits<double>::max();

    max_pt_[0] = std::numeric_limits<double>::min();
    max_pt_[1] = std::numeric_limits<double>::min();
    max_pt_[2] = std::numeric_limits<double>::min();

    std::cout << " input pt size : " << ptn << std::endl;

    for(size_t i = 0; i < ptn; ++i)
    {
        min_pt_[0] = min_pt_[0] < pts[i][0] ? min_pt_[0] : pts[i][0];
        min_pt_[1] = min_pt_[1] < pts[i][1] ? min_pt_[1] : pts[i][1];
        min_pt_[2] = min_pt_[2] < pts[i][2] ? min_pt_[2] : pts[i][2];

        max_pt_[0] = max_pt_[0] > pts[i][0] ? max_pt_[0] : pts[i][0];
        max_pt_[1] = max_pt_[1] > pts[i][1] ? max_pt_[1] : pts[i][1];
        max_pt_[2] = max_pt_[2] > pts[i][2] ? max_pt_[2] : pts[i][2];
    }
}

void OctBBox::GetCornersFromPts(const Point& p1,const Point& p2)
{
    min_pt_[0] = p1[0] < p2[0] ? p1[0] : p2[0];
    min_pt_[1] = p1[1] < p2[1] ? p1[1] : p2[1];
    min_pt_[2] = p1[2] < p2[2] ? p1[2] : p2[2];

    max_pt_[0] = p1[0] > p2[0] ? p1[0] : p2[0];
    max_pt_[1] = p1[1] > p2[1] ? p1[1] : p2[1];
    max_pt_[2] = p1[2] > p2[2] ? p1[2] : p2[2];
}


std::array<double,3> TreeNode::GetCenter() const
{
    std::array<double,3> center;
    center[0] = (bbox_.min_pt_[0] + bbox_.max_pt_[0]) * 0.5;
    center[1] = (bbox_.min_pt_[1] + bbox_.max_pt_[1]) * 0.5;
    center[2] = (bbox_.min_pt_[2] + bbox_.max_pt_[2]) * 0.5;
    return center;
}

void SimpleOctree::InitOctTree(const std::vector<Point>& pts, int depth)
{
    size_t ptn = pts.size();
    pts_ = pts;
    std::vector<size_t> new_ids;
    for(size_t i = 0; i < ptn; ++i)
    {
        new_ids.push_back(i);
    }
    // std::cout << " new ids size " << ptn << std::endl;
    root_node_ = std::make_shared<TreeNode>();
    root_node_->bbox_.GetCornersFromPts(pts_);
    root_node_->depth_ = 0;
    max_depth_ = depth;

    // root_node_->bbox_.SaveBoxMesh("box_mesh.obj");

    DivideNode(root_node_, new_ids);

    // std::cout << "leaf_pids_ size : " << leaf_pids_.size() << std::endl;
} 

void SimpleOctree::DivideNode(std::shared_ptr<TreeNode> node, const std::vector<size_t>& pids)
{
    if(node->depth_ > max_depth_) return;
    std::array<double,3> center = node->GetCenter();
    auto child_boxes = node->bbox_.DividBBox();


    for(int i = 0; i < 8; ++i)
    {
        std::vector<size_t> child_pids;
        node->childNodes[i] = std::make_shared<TreeNode>();
        node->childNodes[i]->parentNode = node;
        node->childNodes[i]->depth_ = node->depth_ + 1;
        node->childNodes[i]->bbox_ = child_boxes[i]; 
        auto& cur_box = child_boxes[i]; 
        for(const auto pid : pids)
        {
            if(IsInsideBox(pts_[pid], cur_box))
            {
                child_pids.push_back(pid);
            } 
        }
        if(node->childNodes[i]->depth_ == max_depth_ || child_pids.empty())
        {
            node->childNodes[i]->is_leaf_ = true;
            leaf_pids_[node->childNodes[i]] = child_pids;
            continue;
        }
        DivideNode(node->childNodes[i], child_pids);
    }

}



void SimpleOctree::GetLeafPts()
{
    leaf_pts_.clear();
    for(auto &ele : leaf_pids_)
    {
        Point closet_pt;
        auto cur_node = ele.first;
        auto pids = ele.second;
        auto center = cur_node->GetCenter();
        double min_dist = std::numeric_limits<double>::max();
        if(pids.empty()) continue;
        for(const auto& pid : pids)
        {
            double cur_dist = SquareDistance(center, pts_[pid]);
            min_dist = min_dist < cur_dist? min_dist : cur_dist;
            closet_pt = pts_[pid];
        }
        leaf_pts_.push_back(closet_pt);
    }
}

}