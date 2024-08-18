#include "local_vipss.hpp"

struct volume_pt
{
    double data[3]; 
};


class Volume {
    public: 
        Volume(){};
        ~Volume(){delete[] data_;};
        void SetDim(int dim);
        void SetDims(int nx, int ny, int nz);
        void InitData();
        double GetVal(int x_id, int y_id, int z_id) const;
        volume_pt GetPt(int x_id, int y_id, int z_id) const;
        void AssignVal(int x_id, int y_id, int z_id, double val) const;
        void AssignOrigin(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z);
        void AssignValsWithLocalVipss(LocalVipss& local_vipss);
    
    public:
        int dim_x_ = 100;
        int dim_y_ = 100;
        int dim_z_ = 100;
        double* data_;
        double origin_x_ = 0;
        double origin_y_ = 0;
        double origin_z_ = 0;
        double voxel_len_x_ = 0.01;
        double voxel_len_y_ = 0.01;
        double voxel_len_z_ = 0.01;
        int voxel_num_;

}; 