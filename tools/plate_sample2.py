import numpy as np
import os

# Function to save points and normals into a .xyz file
def save_xyz(file_path, points):
    # Concatenate points and normals (each row will have [x, y, z, nx, ny, nz])
    # Save as a space-separated .xyz file
    np.savetxt(file_path, points, fmt='%.6f')  # Saving with 6 decimal precision

def load_xyz(file_path):
    data = np.loadtxt(file_path)
    points = data[:, :6]  # Extract first three columns as point coordinates (x, y, z)
    return points

plate_disk10k_path = r'./data/thin_plate/plate_2000.xyz'
plate_disk500_path = r'./data/thin_plate/plate_500n.xyz'

save_dir = r'./data/thin_plate'

plate_disk10k   = load_xyz(plate_disk10k_path)
plate_disk500   = load_xyz(plate_disk500_path)

# for i in range(4):
#     for j in range(i+1,4):
#         # print(i, " ", j)
#         sample_pts = []
#         sample_name = name_list[i] + '_' + name_list[j] + '.xyz'
#         for pt in torus_list[i]:
#             if pt[0] > 0:
#                 sample_pts.append(pt)
#         for pt in torus_list[j]:
#             if pt[0] < 0:
#                 sample_pts.append(pt)
#         out_path = os.path.join(save_dir, sample_name)
#         save_xyz(out_path, sample_pts)

sample_pts = []
sample_name =  'plate_comb_1000.xyz'
for pt in plate_disk10k:
    if pt[1] < 0 :
        sample_pts.append(pt)
for pt in plate_disk500:
    if pt[1] > 0:
        sample_pts.append(pt)

out_path = os.path.join(save_dir, sample_name)
save_xyz(out_path, sample_pts)
