import numpy as np
import os

# Read the .xyz file
def read_xyz(file_path):
    # Load the point cloud data from the file
    points = np.loadtxt(file_path, )
    return points[:,:3]

# Generate random normals and normalize them
def generate_random_normals(num_points, points, axis):
    # Generate random values from a normal distribution
    normals = np.random.normal(size=(num_points, 3))
    
    normals[:,:3] = points[:, :3]
    normals[:, 1] *= 0

    print(normals)


    
    # Normalize each normal vector to unit length
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / norms * 0.75 * 2
    normals = points[:, :3] - normals

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
# input_file_path = r"c:\Users\xiaji\Documents\projects\sketches_results\ipsr_torus_out_normal.xyz"  # Replace with your input file path
input_file_path = r'/home/jjxia/Documents/projects/3D_pointcloud_dataset/results/gt/torus/circle_torus_points_n.xyz'
output_file_path = r'/home/jjxia/Documents/projects/3D_pointcloud_dataset/results/gt/torus/circle_torus_points_n2.xyz'  # Replace with your desired output file path

# # Read the input point cloud
points = read_xyz(input_file_path)

axis = 1
# # Generate random normals for each point
normals = generate_random_normals(len(points), points, axis )

# # Save the points and generated normals to the output file
save_xyz_with_normals(points, normals, output_file_path)

print(f"Point cloud with normals saved to {output_file_path}")

# file_dir = r'C:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus_sample\sampling'
# out_dir = r'C:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\results\gt\torus_sample'
# namelist = os.listdir(file_dir)

# for name in namelist:
    # file_path = os.path.join(file_dir, name)
    # out_path = os.path.join(out_dir, name)

# file_path = './data/torus/torus_two_parts.xyz'
# out_path = './data/torus/torus_two_parts_normal.xyz'
# points = read_xyz(file_path)
# # Generate random normals for each point
# normals = generate_random_normals(len(points), points)

# # Save the points and generated normals to the output file
# save_xyz_with_normals(points, normals, out_path)