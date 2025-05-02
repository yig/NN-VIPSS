import numpy as np
import os

def read_xyz_with_normals(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                values = list(map(float, line.strip().split()))
                if len(values) == 6:
                    data.append(values)
                else:
                    raise ValueError("Each line must have 6 values: x y z nx ny nz")
    return np.array(data)

def scale_points_uniform(points, target_min=-1.0, target_max=1.0):
    center = (points.max(axis=0) + points.min(axis=0)) / 2.0
    centered = points - center
    scale = np.max(np.abs(centered))
    scaled = centered / (scale + 1e-8)  # avoid division by zero
    return scaled

def write_xyz_with_normals(filename, points, normals):
    with open(filename, 'w') as file:
        for p, n in zip(points, normals):
            file.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f} {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n")

# # Main usage
# input_file = 'input.xyz'
# output_file = 'output_scaled.xyz'
# # Example usage
# input_file = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\elbow_scan\6007_9_9.xyz'
# output_file = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\elbow_scan\6007_9_9_normalized.xyz'

input_dir = '/home/jjxia/Documents/projects/VIPSS_LOCAL/data_old/new_gt'
out_dir = '/home/jjxia/Documents/projects/VIPSS_LOCAL/data/gt_normal'

file_names = os.listdir(input_dir)
for file_name in file_names : 

    input_file = os.path.join(input_dir, file_name)
    output_file = os.path.join(out_dir, file_name)
    
    data = read_xyz_with_normals(input_file)
    points = data[:, :3]
    normals = data[:, 3:]

    scaled_points = scale_points_uniform(points)
    write_xyz_with_normals(output_file, scaled_points, normals)

    # file_path = os.path.join(input_dir, file_name)
    # out_file = os.path.join(out_dir, file_name)
    # # Read the original XYZ file
    # atoms, coords = read_xyz(file_path)
    # # Normalize the coordinates
    # normalized_coords = normalize_coordinates(coords)
    # # Write the normalized data back to a new XYZ file
    # write_xyz(out_file,  normalized_coords)
    # print(f"Normalized XYZ file written to {out_file}")