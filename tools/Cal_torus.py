import trimesh
import numpy as np

# Load the point cloud from an .obj file
# point_cloud = trimesh.load(r"c:\Users\xiaji\Documents\projects\sketches_results\ipsr_torus.xyz", file_type='xyz')
# path = "./data/torus/gt.obj"

path = "./data/torus/gt.obj"

point_cloud = trimesh.load(path)

# Get the axis-aligned bounding box (AABB)
bounding_box = point_cloud.bounds

# The bounding_box will return two corners of the AABB: [min(x, y, z), max(x, y, z)]
min_corner = bounding_box[0]
max_corner = bounding_box[1]

# Print the bounding box
print(f"Min corner: {min_corner}")
print(f"Max corner: {max_corner}")

# You can also compute the dimensions of the bounding box
bbox_dimensions = max_corner - min_corner
print(f"Bounding box dimensions: {bbox_dimensions}")