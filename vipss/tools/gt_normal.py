import numpy as np
import os

# Read the .xyz file
def read_xyz(file_path):
    # Load the point cloud data from the file
    points = np.loadtxt(file_path, )
    return points

# Generate random normals and normalize them
def generate_random_normals(num_points, points):
    # Generate random values from a normal distribution
    normals = np.random.normal(size=(num_points, 3))
    normals[:, 2] *= 0
    normals[:, :2] = points[:, :2]
    
    # Normalize each normal vector to unit length
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / norms  * 1.0
    normals = points - normals

    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / norms 
    return normals


# Save points and normals to a new file (e.g., .xyz format with normals)
def save_xyz_with_normals(points, normals, output_file_path):
    # Concatenate points and normals horizontally
    data_with_normals = np.hstack((points, normals))
    # Save the result with six decimal places
    np.savetxt(output_file_path, data_with_normals, fmt="%.6f", )

# Define file paths
input_file_path = r"c:\Users\xiaji\Documents\projects\sketches_results\ipsr_torus_out_normal.xyz"  # Replace with your input file path
output_file_path = r"c:\Users\xiaji\Documents\projects\sketches_results\ipsr_torus_gt.xyz"  # Replace with your desired output file path

# # Read the input point cloud
# points = read_xyz(input_file_path)

# # Generate random normals for each point
# normals = generate_random_normals(len(points), points)

# # Save the points and generated normals to the output file
# save_xyz_with_normals(points, normals, output_file_path)

print(f"Point cloud with normals saved to {output_file_path}")

file_dir = r'C:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus_sample\sampling'
out_dir = r'C:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\results\gt\torus_sample'
namelist = os.listdir(file_dir)

for name in namelist:
    file_path = os.path.join(file_dir, name)
    out_path = os.path.join(out_dir, name)
    points = read_xyz(file_path)
    # Generate random normals for each point
    normals = generate_random_normals(len(points), points)

    # Save the points and generated normals to the output file
    save_xyz_with_normals(points, normals, out_path)