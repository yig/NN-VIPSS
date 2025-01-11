import numpy as np
import os
import trimesh

from scipy.spatial import cKDTree
# Function to save points and normals into a .xyz file
def save_xyz_with_normals(file_path, points, normals):
    # Concatenate points and normals (each row will have [x, y, z, nx, ny, nz])
    data = np.hstack((points, normals))
    
    # Save as a space-separated .xyz file
    np.savetxt(file_path, data, fmt='%.6f')  # Saving with 6 decimal precision

def load_xyz(file_path):
    data = np.loadtxt(file_path)
    points = data[:, :3]  # Extract first three columns as point coordinates (x, y, z)
    return points

# Load the original point cloud with normals from the .xyz file
def load_xyz_with_normals(file_path):
    data = np.loadtxt(file_path)
    points = data[:, :3]  # First three columns are the point coordinates (x, y, z)
    normals = data[:, 3:6]  # Next three columns are the normals (nx, ny, nz)
    return points, normals

# Load the original point cloud with normals
infile = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\thin_plate_400k.xyz'


# Build a KDTree from the original point cloud



# def ipsr_pt_gt_normal():

out_dir = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\results\gt\high_genus'
in_dir = r'C:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\dense\high_genus\withgt'

filenames = os.listdir(in_dir)

in_dir = r'C:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\results\ipsr\high_genus\withgt'

for name in filenames : 

    infile = os.path.join(out_dir, name)
    original_points, original_normals = load_xyz_with_normals(infile)
    kdtree = cKDTree(original_points)

    newname = name.split('.')[0] + ".ply_n.ply"
    infile = os.path.join(in_dir, newname)
    # sampled_points = load_xyz(infile)
    mesh = trimesh.load(infile)
    # Extract the vertex positions (which are the point cloud data)
    sampled_points = mesh.vertices
    # Find the closest original point for each sampled point
    _, closest_indices = kdtree.query(sampled_points)

    # Get the normals of the nearest points in the original point cloud
    sampled_normals = original_normals[closest_indices]
    outfile = os.path.join(out_dir, name.split('.')[0] + "ipsr.xyz")
    # Save the sampled points and normals to a new .xyz file
    save_xyz_with_normals(outfile, sampled_points, sampled_normals)

    print("Sampled point cloud with normals saved to .xyz file!")


# psr_pt_gt_normal()