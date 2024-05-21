#include "orient_normal.h"

namespace ORIENT {

    void OrientPointNormals(std::vector<double>& pts, std::vector<double>& normals)
    {
        std::list<PointVectorPair> points;
        for(size_t i = 0; i < pts.size()/3; ++i)
        {
            Point cur_p(pts[3*i], pts[3*i + 1], pts[3*i + 2]);
            Vector cur_n(normals[3*i], normals[3*i + 1], normals[3*i + 2]);

            PointVectorPair pn(cur_p, cur_n);
            points.push_back(pn);
        }

        // Orients normals.
        // Note: mst_orient_normals() requires an iterator over points
        // as well as property maps to access each point's position and normal.
        int nb_neighbors = 6;
        std::list<PointVectorPair>::iterator unoriented_points_begin =
            CGAL::mst_orient_normals(points, nb_neighbors,
                                    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                                        .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
        // Optional: delete points with an unoriented normal
        // if you plan to call a reconstruction algorithm that expects oriented normals.
        points.erase(unoriented_points_begin, points.end());

        pts.clear();
        normals.clear();
        for(const auto& ele : points)
        {
            pts.push_back(ele.first.x());
            pts.push_back(ele.first.y());
            pts.push_back(ele.first.z());

            normals.push_back(ele.second.x());
            normals.push_back(ele.second.y());
            normals.push_back(ele.second.z());
        }
        

    }
} 