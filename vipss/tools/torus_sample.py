import numpy as np
import os

# Function to save points and normals into a .xyz file
def save_xyz(file_path, points):
    # Concatenate points and normals (each row will have [x, y, z, nx, ny, nz])

    
    # Save as a space-separated .xyz file
    np.savetxt(file_path, points, fmt='%.6f')  # Saving with 6 decimal precision

def load_xyz(file_path):
    data = np.loadtxt(file_path)
    points = data[:, :3]  # Extract first three columns as point coordinates (x, y, z)
    return points



torus_vcircle_path = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus_sample\blender\torus_mont256.xyz'
torus_hcircle_path = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus_sample\blender\torus_mont2560.xyz'
torus_disk_path  = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus_sample\blender\torus_mont25600.xyz'
torus_mont_path = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus_sample\blender\torus_mont256000.xyz'

save_dir = r'C:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus_sample\sample'

torus_vcircle = load_xyz(torus_vcircle_path)
torus_hcircle = load_xyz(torus_hcircle_path)
torus_disk = load_xyz(torus_disk_path)
torus_mont = load_xyz(torus_mont_path)

name_list = ['mont256',  'mont2k', 'mont25k', 'mont256k']

torus_list = [torus_vcircle,  torus_hcircle,  torus_disk, torus_mont]


for i in range(4):
    for j in range(i+1,4):
        # print(i, " ", j)
        sample_pts = []
        sample_name = name_list[i] + '_' + name_list[j] + '.xyz'
        for pt in torus_list[i]:
            if pt[0] > 0:
                sample_pts.append(pt)
        for pt in torus_list[j]:
            if pt[0] < 0:
                sample_pts.append(pt)
        out_path = os.path.join(save_dir, sample_name)
        save_xyz(out_path, sample_pts)



sample_pts = []
sample_name =  'combine_all.xyz'
for pt in torus_list[0]:
    if pt[0] > 0 and pt[1] > 0:
        sample_pts.append(pt)
for pt in torus_list[1]:
    if pt[0] < 0 and pt[1] > 0:
        sample_pts.append(pt)
for pt in torus_list[2]:
    if pt[0] < 0 and pt[1] < 0:
        sample_pts.append(pt)
for pt in torus_list[3]:
    if pt[0] > 0 and pt[1] < 0:
        sample_pts.append(pt)
out_path = os.path.join(save_dir, sample_name)
save_xyz(out_path, sample_pts)
