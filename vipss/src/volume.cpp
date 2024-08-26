#include "volume.h"
#include <algorithm>

void Volume::SetDim(int dim)
{
    dim_x_ = dim;
    dim_y_ = dim;
    dim_z_ = dim;
    InitData();
}

void Volume::SetDims(int nx, int ny, int nz)
{
    dim_x_ = nx;
    dim_y_ = ny;
    dim_z_ = nz;
    InitData();
}

void Volume::InitData()
{
    voxel_num_ = dim_x_ * dim_y_ * dim_z_;
    data_ = new double(voxel_num_);
}

double Volume::GetVal(int x_id, int y_id, int z_id) const
{
    int id = z_id * (dim_x_ * dim_y_) + y_id * dim_x_ + x_id;
    return data_[id];
}

volume_pt Volume::GetPt(int x_id, int y_id, int z_id) const
{
    volume_pt vx_pt;
    vx_pt.data[0] = origin_x_ + (x_id + 0.5) * voxel_len_x_;
    vx_pt.data[1] = origin_y_ + (y_id + 0.5) * voxel_len_y_;
    vx_pt.data[2] = origin_z_ + (z_id + 0.5) * voxel_len_z_;
    return vx_pt;
}

void Volume::AssignVal(int x_id, int y_id, int z_id, double val) const
{
    int id = z_id * (dim_x_ * dim_y_) + y_id * dim_x_ + x_id;
    data_[id] = val;
}

void Volume::AssignOrigin(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z)
{
    double cx = (min_x + max_x) / 2.0;
    double cy = (min_y + max_y) / 2.0;
    double cz = (min_z + max_z) / 2.0;

    double len_x = (max_x - min_x);
    double len_y = (max_y - min_y);
    double len_z = (max_z - min_z);

    voxel_len_x_ = len_x / dim_x_;
    voxel_len_y_ = len_y / dim_y_;
    voxel_len_z_ = len_z / dim_z_;
    double scale = 1.5 * 0.5;
    origin_x_ = cx - len_x * scale;
    origin_y_ = cy - len_y * scale;
    origin_z_ = cz - len_z * scale;
}


void Volume::AssignValsWithLocalVipss(LocalVipss& local_vipss)
{
    for(int ix = 0; ix < dim_x_; ++ix)
    {
        for(int iy = 0; iy < dim_y_; ++iy)
        {
            for(int iz = 0; iz < dim_z_; ++iz)
            {
                auto vpt = GetPt(ix, iy, iz);
                local_vipss.NatureNeighborDistanceFunction(&vpt.data[0]);
            }
        }
    }
}