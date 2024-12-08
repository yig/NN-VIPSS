import trimesh
import numpy as np
import os

# Load your 3D mesh file (replace 'your_mesh.obj' with your file)
mesh_path = r'./data/kitten.ply'
mesh = trimesh.load(mesh_path)
out_dir = r'./data/kitten'

# Number of points to sample
nums = list(map(lambda x: x * 10000,range(1, 11)))
names = list(map(lambda x: 'kitten_'+str(x)+'0k.xyz', range(1, 11)))

for num, name in zip(nums, names):
# Sample points evenly from the surface of the mesh
    points, _ = trimesh.sample.sample_surface_even(mesh, num)
    out_path = os.path.join(out_dir, name)
    # Save the sampled points to a file (optional)
    np.savetxt(out_path, points)

# # Visualize the points (optional)
# cloud = trimesh.points.PointCloud(points)
# cloud.show()
