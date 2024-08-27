// #pragma once

// #include <CGAL/Simple_cartesian.h>
// #include <CGAL/point_generators_3.h>
// #include <CGAL/Orthogonal_k_neighbor_search.h>
// #include <CGAL/Search_traits_3.h>
// #include <CGAL/Euclidean_distance.h>
// #include <list>
// #include <cmath>
// #include <CGAL/pca_estimate_normals.h>
// #include <CGAL/jet_estimate_normals.h>
// #include <CGAL/vcm_estimate_normals.h>
// #include <CGAL/mst_orient_normals.h>
// #include <CGAL/property_map.h>
// #include <CGAL/IO/read_points.h>
// #include <CGAL/IO/write_ply_points.h>
// #include <CGAL/IO/write_xyz_points.h>
// #include <CGAL/compute_average_spacing.h>
// #include <CGAL/jet_smooth_point_set.h>



// namespace MKdtree{

// typedef CGAL::Simple_cartesian<double> K;
// typedef K::Point_3 Point_d;
// typedef K::Vector_3 Vector_d;
// typedef CGAL::Search_traits_3<K> TreeTraits;
// typedef CGAL::Euclidean_distance<TreeTraits> Distance;

// typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
// typedef Neighbor_search::Tree Tree;

// // Point with normal vector stored in a std::pair.
// typedef std::pair<Point_d, Vector_d> PointVectorPair;
// // Concurrency
// typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// void EstimatePtNormals(const std::string& fname, size_t nb_neighbors,
//                         std::vector<double>& pts, std::vector<double>& normals,
//                         double& point_avg_space);


// void EstimatePtAvgDist(const std::vector<double>& pts, double& point_avg_space);

// // const int D = 3;
// // typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > K;
// // typedef K::Point_d Point_d;
// // typedef CGAL::Search_traits_d<K,CGAL::Dimension_tag<D> >  Traits;
// // typedef CGAL::Euclidean_distance<Traits> Distance;
// // typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
// // // typedef CGAL::Kd_tree<Traits> Tree;
// // typedef K_neighbor_search::Tree                                         Tree;



// class MyKDTree
// {
//     public:
//         MyKDTree() {};
//         ~ MyKDTree() {};

//         void InitTree(const std::vector<double> &pts);
//         bool InsertPt(double x, double y, double z, double dist_threshold);
//         bool InsertPt(Point_d& pt, double dist_threshold);
//         bool FindNearestPt(double x, double y, double z, Point_d& pt, double& dist);
//         std::vector<Point_d> NeighborSearchRadius(double x, double y, double z, double radius);

//         double CalPtSampleDist(double x, double y, double z);

//         void RemovePt(double x, double y, double z);
//         std::vector<double> GetTreePts() const;
    

//     private:
//         void InsertPt(const Point_d& pt);
//         void InsertPt(double x, double y, double z);

//     public: 
//         std::vector<double> origin_pts_;
//         std::shared_ptr<Tree> tree_ptr_;


// };

// }

