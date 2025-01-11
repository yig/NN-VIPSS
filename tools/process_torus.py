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

torus_path = r'/home/jjxia/Downloads/torus_500.xyz'
save_dir = r'./data/torus'

torus_pts  = load_xyz(torus_path)


sample_pts = []
sample_name =  'torus_two_parts.xyz'
for pt in torus_pts:
    if pt[1] < -0.33 or pt[1] > 0.33:
        sample_pts.append(pt)

out_path = os.path.join(save_dir, sample_name)
save_xyz(out_path, sample_pts)

