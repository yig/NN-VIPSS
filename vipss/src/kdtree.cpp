#include "kdtree.hpp"

namespace MKdtree{

    void EstimatePtAvgDist(const std::vector<double>& pts, double& point_avg_space)
    {
        std::list<Point_d> points;
        for(size_t i = 0; i < pts.size()/3; ++i)
        {
            Point_d cur_p(pts[3*i], pts[3*i + 1], pts[3*i + 2]);
            points.push_back(cur_p);
        }

        
        int nb_neighbors = 4;
        point_avg_space = CGAL::compute_average_spacing<Concurrency_tag>
          (points, nb_neighbors);
           
    }

    void OrientPtNormals( size_t nb_neighbors, std::vector<double>& pts, std::vector<double>& normals)
    {
        std::list<PointVectorPair> points;
        for(size_t i = 0; i < pts.size()/3; ++i)
        {
            Point_d cur_p(pts[3*i], pts[3*i + 1], pts[3*i + 2]);
            Vector_d cur_n(normals[3*i], normals[3*i + 1], normals[3*i + 2]);

            points.push_back(PointVectorPair(cur_p, cur_n));
        }
        std::list<PointVectorPair>::iterator unoriented_points_begin =
            CGAL::mst_orient_normals(points, nb_neighbors,
                                    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                                        .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
        // Optional: delete points with an unoriented normal
        // if you plan to call a reconstruction algorithm that expects oriented normals.
        points.erase(unoriented_points_begin, points.end());

        

    }

    void EstimatePtNormals(const std::string& fname, size_t nb_neighbors, 
                            std::vector<double>& pts, std::vector<double>& normals, 
                            double& point_avg_space)
    {
        std::list<PointVectorPair> points;
        if(!CGAL::IO::read_points(fname,
                                    std::back_inserter(points),
                                    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())))
        {
            std::cerr << "Error: cannot read file " << fname<< std::endl;
            // return EXIT_FAILURE;
        }

        double spacing
        = CGAL::compute_average_spacing<Concurrency_tag>
          (points, nb_neighbors,
           CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()));
        point_avg_space = spacing;
        CGAL::jet_estimate_normals<Concurrency_tag>
            (points, nb_neighbors,
            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                            .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
        // CGAL::IO::write_PLY("oni_copy.ply", points,
        //                     CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
        //                                     .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())
        //                                     .stream_precision(17));
        
        // Orients normals.
        // Note: mst_orient_normals() requires a range of points
        // as well as property maps to access each point's position and normal.
        std::list<PointVectorPair>::iterator unoriented_points_begin =
            CGAL::mst_orient_normals(points, nb_neighbors,
                                    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                                        .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
        // Optional: delete points with an unoriented normal
        // if you plan to call a reconstruction algorithm that expects oriented normals.
        points.erase(unoriented_points_begin, points.end());

        unoriented_points_begin =
            CGAL::mst_orient_normals(points, nb_neighbors,
                                    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                                        .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
        points.erase(unoriented_points_begin, points.end());
        

        for(auto& p : points)
        {
            pts.push_back(p.first.x());
            pts.push_back(p.first.y());
            pts.push_back(p.first.z());

            normals.push_back(p.second.x());
            normals.push_back(p.second.y());
            normals.push_back(p.second.z());

            // std::cout << p.first << std::endl;
            // std::cout << p.second << std::endl;
            // break;
        }

        // CGAL::IO::write_XYZ("oni_copy.xyz", points,
                            //   CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                            //                    .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())
                            //                    .stream_precision(17));
    }


void MyKDTree::InitTree(const std::vector<double> &pts)
{
    origin_pts_ = pts;
    std::vector<Point_d> points;
    for(size_t i = 0; i < pts.size() / 3; ++i )
    {
        Point_d p(pts[3*i], pts[3*i + 1], pts[3*i + 2]);
        points.push_back(p);
    }
    tree_ptr_ = std::make_shared<Tree>(points.begin(), points.end());
}


void MyKDTree::InsertPt(const Point_d& pt)
{
    tree_ptr_->insert(pt);
    origin_pts_.push_back(pt.x());
    origin_pts_.push_back(pt.y());
    origin_pts_.push_back(pt.z());
}


void MyKDTree::InsertPt(double x, double y, double z)
{
    Point_d p(x, y, z);
    InsertPt(p);
    // tree_ptr_->insert(p);
}

bool MyKDTree::InsertPt(double x, double y, double z, double dist_threshold)
{
    Point_d target_p;
    double dist;
    if(FindNearestPt(x, y, z, target_p, dist))
    {
        if(dist >= dist_threshold)
        {
            InsertPt(x, y, z);
            return true;
        }
    }
    return false;
}

bool MyKDTree::InsertPt(Point_d& pt, double dist_threshold)
{
    return InsertPt(pt.x(), pt.y(), pt.z(), dist_threshold);
}

bool MyKDTree::FindNearestPt(double x, double y, double z, Point_d& target_p, double& dist)
{
    Point_d query(x, y, z);
    Neighbor_search search(*tree_ptr_, query, 1);
    Neighbor_search::iterator it = search.begin();

    if(it != search.end())
    {
        target_p = it->first;
        dist = std::sqrt(it->second);
        
        return true;    
    }

    return false;
}

std::vector<Point_d> MyKDTree::NeighborSearchRadius(double x, double y, double z, double radius)
{
    Point_d query(x, y, z);
    size_t K = origin_pts_.size();
    // std::cout << "K size " << K << std::endl;
    Neighbor_search N1(*tree_ptr_, query, K);
    Neighbor_search::iterator it = N1.begin();

    std::vector<Point_d> nearestPts; 
    double squred_dist = radius * radius;
    // std::cout << "squred_dist  " << squred_dist << std::endl;
    while(it != N1.end())
    {

        // std::cout << "dist : " << it->second << std::endl;
        if(it->second <= squred_dist)
        {
            nearestPts.push_back(it->first);
        } else {
            // break;
        }
        it++;
    }
    return nearestPts;
}

double MyKDTree::CalPtSampleDist(double x, double y, double z)
{
    Point_d query(x, y, z);
    // std::cout << "-------------query pt " << query << std::endl;
    // std::cout << "pt num " << origin_pts_.size() / 3 << std::endl;
    Neighbor_search search(*tree_ptr_, query, 3);
    Neighbor_search::iterator it = search.begin();
    double dist_sum = 0;
    size_t count = 0;
    while(it != search.end())
    {
        // std::cout << "dist sum " << it->second << std::endl;
        if(sqrt(it->second) > 1e-8)
        {
            dist_sum += sqrt(it->second);
            count ++;
        }
        ++it;
    }
    if(count > 0)
    {
        dist_sum = dist_sum / double(count); 
    }
    return dist_sum;
}

void MyKDTree::RemovePt(double x, double y, double z)
{
    Point_d p(x, y, z);
    tree_ptr_->remove(p);
}

std::vector<double> MyKDTree::GetTreePts() const
{
    auto pt_iter = tree_ptr_->begin();
    std::vector<double> kdPts;
    while(pt_iter != tree_ptr_->end())
    {
        auto& cur_pt = *pt_iter;
        kdPts.push_back(cur_pt.x());
        kdPts.push_back(cur_pt.y());
        kdPts.push_back(cur_pt.z());
        pt_iter ++;
    }
    return kdPts;
}


}