#pragma once
#include <array>
#include <vector>
#include <string>
#include <unordered_map> 
#include <memory>


namespace SimOctree{

typedef std::array<double,3> Point;

// double Distance(const Point& p1, const Point& p2);

struct OctBBox{
    Point min_pt_;
    Point max_pt_;

    OctBBox(){};

    OctBBox(const Point& minP, const Point& maxP);
    void GetCornersFromPts(const std::vector<double>& pts);
    void GetCornersFromPts(const std::vector<Point>& pts);
    void GetCornersFromPts(const Point& p1,const Point& p2);
    std::vector<OctBBox> DividBBox() const; 
    void SaveBoxMesh(const std::string& path) const;
    std::vector<Point> GenerateCorners() const;  
};

inline bool IsInsideBox(const Point& pt, const OctBBox& box);

struct TreeNode{
    OctBBox bbox_;
    int depth_;
    Point GetCenter() const;
    std::shared_ptr<TreeNode> parentNode; 
    std::shared_ptr<TreeNode> childNodes[8];
    bool is_leaf_ = false;
};

typedef std::unordered_map<std::shared_ptr<TreeNode>, std::vector<size_t>> LeafNodeMap;
class SimpleOctree{
    public:
        SimpleOctree() {};
        ~SimpleOctree() {};

        // void InitOctTree(const std::vector<Point>& pts, int depth = 5);
        void InitOctTree(const std::vector<double>& pts, int depth = 5);
        void DivideNode(std::shared_ptr<TreeNode> node, const std::vector<size_t>& pids);
        std::vector<double> GetLeafMapPts();
        std::vector<double> GetSplitLeafMapPts();
        void SplitLeafNode(int iter_num);
    private: 
        std::vector<double> GetLeafPts(LeafNodeMap& leaf_pids);
        LeafNodeMap SplitLeafNode(LeafNodeMap& leaf_pids);

        // void SplitLeafNode();
        double *pt_data_;
    public:
        // std::unordered_map<SimpleOctTree*, int> depth;
        std::unordered_map<std::shared_ptr<TreeNode>, std::vector<size_t> > leaf_pids_;
        std::shared_ptr<TreeNode> root_node_;
        // std::vector<Point> pts_;
        int max_depth_ = 6; 
        // std::vector<Point> leaf_pts_;
        // std::vector<Point> split_leaf_pts_;
        std::unordered_map<std::shared_ptr<TreeNode>, std::vector<size_t> > split_leaf_pids_;
};

}