import numpy as np
import trimesh
import os
import math
import kdtree


from plyfile import PlyData, PlyElement
# mesh = trimesh.load_mesh(mesh_dir)

def write_ply(pts, ptns, out_dir):
    p_num = len(pts)
    new_vertices = []
    for i in range(p_num):
        pt = (pts[i].x(),pts[i].y(), pts[i].z(), ptns[i].x(),ptns[i].y(), ptns[i].z())
        new_vertices.append(pt)
    new_vertices = np.array(new_vertices, dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
                                                 ('nx', 'f4'), ('ny', 'f4'), ('nz', 'f4')])
    pt_data = PlyElement.describe(new_vertices, 'vertex')
    PlyData([pt_data], text=True).write(out_dir)


def read_ptns(ply_path):
    plydata = PlyData.read(ply_path)
    pt_num = plydata.elements[0].count
    pt_data = plydata.elements[0]
    ptns = []
    for i in range(pt_num):
        cur_pt = list(pt_data[i])
        ptns.append(cur_pt)

    ptns = np.array(ptns)
    
    # print("ptns size ", ptns.shape)
    # print(ptns)
    return ptns

        
def data_convert():
    # pt_dir1 = r'../data/dragon/dragon_sample1k.ply'
    # pt_dir1 = r'../data/doghead_100/doghead_100.ply'
    # pt_dir1 = r'../data/arma_100k/arma_100k.ply'
    pt_dir1 = r'../data/dragon/dragon.ply'
    # pt_dir1 = r'/home/jjxia/Documents/projects/VIPSS_SP/VIPSS_M/data/test/dragon_stand/dragon32k.ply'
    ptns1 = read_ptns(pt_dir1)
    npt = ptns1.shape[0]
    # noise = np.random.normal(0,npt,3) 
    noise = np.random.normal(0, 0.01, (npt,3))
    ptns1[:,:3] = ptns1[:,:3] + noise
    new_vertices = []
    for i in range(npt):
        ptn = ptns1[i, :]
        pt = (ptn[0], ptn[1], ptn[2], ptn[3], ptn[4], ptn[5])
        # print(pt)
        new_vertices.append(pt)
    
    new_vertices = np.array(new_vertices, dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
                                                 ('nx', 'f4'), ('ny', 'f4'), ('nz', 'f4')])
    pt_data = PlyElement.describe(new_vertices, 'vertex')
    out_dir = r'../data/dragon/dragon_sample.ply'
    PlyData([pt_data], text=True).write(out_dir)


if __name__ == "__main__":
    data_convert()
    # parse_contours(' ')