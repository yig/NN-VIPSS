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

kitten_sample1_path = r'./data/kitten_sample/kitten_2000_uniform.xyz'
kitten_sample2_path = r'./data/kitten_sample/kitten_250_uniform.xyz'

save_dir = r'./data/kitten_sample/kitten_250u_2000u_combine.xyz'
kitten_sample1  = load_xyz(kitten_sample1_path)
kitten_sample2  = load_xyz(kitten_sample2_path)

name_list = ['disk2k',  'disk40k', 'mont3k', 'mont60k']
# torus_list = [plate_disk1k,  plate_disk100k]



sample_pts = []
# sample_name =  'combine_all.xyz'
for pt in kitten_sample1:
    if pt[0] > 0:
        sample_pts.append(pt)

for pt in kitten_sample2:
    if pt[0] < 0:
        sample_pts.append(pt)

# for pt in torus_list[2]:
#     if pt[0] < 0 and pt[1] < 0:
#         sample_pts.append(pt)
# for pt in torus_list[3]:
#     if pt[0] > 0 and pt[1] < 0:
        # sample_pts.append(pt)
# out_path = os.path.join(save_dir, sample_name)

save_xyz(save_dir, sample_pts)
