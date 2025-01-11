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

plate_disk1k_path = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\plate_sample\blend\plate_disk2k.xyz'
plate_disk100k_path = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\plate_sample\blend\plate_disk40k.xyz'
plate_mont1k_path  = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\plate_sample\blend\plate_mont3k.xyz'
plate_mont100k_path = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\plate_sample\blend\plate_mont60k.xyz'

save_dir = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\plate_sample\blend'

plate_disk1k     = load_xyz(plate_disk1k_path)
plate_disk100k   = load_xyz(plate_disk100k_path)
plate_mont1k     = load_xyz(plate_mont1k_path)
plate_mont100k   = load_xyz(plate_mont100k_path)

name_list = ['disk2k',  'disk40k', 'mont3k', 'mont60k']
torus_list = [plate_disk1k,  plate_disk100k,  plate_mont1k, plate_mont100k]

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
