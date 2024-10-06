import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the .xyz file with points and normals
def load_xyz_with_normals(file_path):
    # Assume file format is x y z nx ny nz
    data = np.loadtxt(file_path)
    points = data[:, :3]   # Extract points (x, y, z)
    normals = data[:, 3:6] # Extract normals (nx, ny, nz)
    return points, normals

# Calculate normal differences
def calculate_normal_differences(normals1, normals2):
    differences = normals1 - normals2
    magnitudes = np.linalg.norm(differences, axis=1)  # Calculate the magnitude of differences
    print(np.max(magnitudes))
    print(np.min(magnitudes))
    return differences, magnitudes

# Visualize points and normals differences using quiver
def visualize_normal_differences(points, normals1, normals2, magnitudes):
    # fig = plt.figure()
    fig = plt.figure(figsize=(8, 6)) 

    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.grid(False)
    ax.set_axis_off()

    diff = np.abs(np.sum(normals1 * normals2, axis=1, keepdims=True))
    new_normals1 =   normals1 * (1 - diff)
    new_normals2 =   normals2 * (1 - diff)

    # Plot points
    autoPlot = ax.scatter(points[:, 0], points[:, 1], points[:, 2], c=magnitudes, cmap='seismic', s=10, alpha=0.9)
    
    # Plot normal vectors from the first set
    ax.quiver(points[:, 0], points[:, 1], points[:, 2],
              new_normals1[:, 0], new_normals1[:, 1], new_normals1[:, 2], length=0.1, color='red', label='Normals 1', alpha=0.7,  linewidth=0.5)
    
    # Plot normal vectors from the second set
    ax.quiver(points[:, 0], points[:, 1], points[:, 2],
              new_normals2[:, 0], new_normals2[:, 1], new_normals2[:, 2], length=0.1, color='blue', label='Normals 2',  alpha=0.3, linewidth=0.5)

    # Add labels and color bar for magnitudes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.colorbar(autoPlot, label="Normal Difference Magnitude")
    ax.view_init(elev=80, azim=0)
    
    # plt.title('Normal Differences Visualization')
    plt.show()
    fig.savefig("custom_size_plot.png", dpi=300)

# File paths
file1 = r'c:\Users\xiaji\Documents\projects\IsoConstraints\results\torus_parallelcut_oriented.xyz'  # Replace with your first file path
file2 = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\torus\torus_parellelcut_n.xyz'  # Replace with your second file path

# Load points and normals from both files
points1, normals1 = load_xyz_with_normals(file1)
points2, normals2 = load_xyz_with_normals(file2)

# Ensure both point clouds have the same number of points
assert points1.shape == points2.shape, "Point clouds must have the same number of points."

# Calculate normal differences and magnitudes
normal_differences, magnitudes = calculate_normal_differences(normals1, normals2)

# Visualize the normal differences
visualize_normal_differences(points1, normals1, normals2, magnitudes)
