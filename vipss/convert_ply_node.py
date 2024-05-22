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

def write_ply_color(pts, ptns, colors, out_dir):
    p_num = len(pts)
    new_vertices = []
    for i in range(p_num):
        pt = (pts[i].x(),pts[i].y(), pts[i].z(), ptns[i].x(),ptns[i].y(), ptns[i].z(), colors[i].x(), colors[i].y(), colors[i].z())
        new_vertices.append(pt)
    new_vertices = np.array(new_vertices, dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
                                                 ('nx', 'f4'), ('ny', 'f4'), ('nz', 'f4'), 
                                                 ('red', 'u1'), ('green', 'u1'),('blue', 'u1')])
    
    pt_data = PlyElement.describe(new_vertices, 'vertex')
    PlyData([pt_data], text=True).write(out_dir)

def write_ply_color( ptns, colors, out_dir):
    p_num = ptns.shape[0]
    new_vertices = []
    for i in range(p_num):
        ptn = ptns[i, :]
        
        pt = (ptn[0], ptn[1], ptn[2], ptn[3], ptn[4], ptn[5], colors[i], 0, 0)
        new_vertices.append(pt)
    new_vertices = np.array(new_vertices, dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
                                                 ('nx', 'f4'), ('ny', 'f4'), ('nz', 'f4'), 
                                                 ('red', 'u1'), ('green', 'u1'),('blue', 'u1')])
    
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
    pt_dir1 = r'../data/doghead200.ply'
    # pt_dir1 = r'/home/jjxia/Documents/projects/VIPSS_SP/VIPSS_M/data/test/dragon_stand/dragon32k.ply'
    ptns1 = read_ptns(pt_dir1)
    pts = ptns1[:,:3] 
    
    f = open(r'../data/doghead200.node','w')
    pt_num = pts.shape[0]
    header = str(pt_num) + ' 3 0 0 \n'
    f.write(str(pt_num) + ' 3 0 0 \n')
    for i in range(pt_num):
        data_line = str(i) + ' ' + str(pts[i][0]) + ' ' + str(pts[i][1]) + ' ' + str(pts[i][2]) + ' \n'
        f.write(data_line)
        # f.write('\n')
    f.close
    
    


if __name__ == "__main__":
    data_convert()
    # parse_contours(' ')