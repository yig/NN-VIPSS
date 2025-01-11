import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters for the 3D circle
radius = 1.0  # Radius of the circle
center = np.array([0.0, 0.0, 0.0])  # Center of the circle
normal = np.array([0.0, 0.0, 1.0])  # Normal vector to the circle's plane

# Function to generate 100 points uniformly on the circle
def sample_circle_3d(center, normal, radius, num_points=100):
    # Normalize the normal vector
    normal = normal / np.linalg.norm(normal)
    
    # Find a vector perpendicular to the normal
    perp_vector = np.array([1, 0, 0]) if abs(normal[0]) < 1.0 else np.array([0, 1, 0])
    u = np.cross(normal, perp_vector)
    u = u / np.linalg.norm(u)
    
    # Find another perpendicular vector
    v = np.cross(normal, u)
    
    # Generate angles for uniform sampling
    angles = np.linspace(0, 2 * np.pi, int(num_points), endpoint=False)
    
    # Generate points on the circle
    points = [
        center + radius * (np.cos(angle) * u + np.sin(angle) * v)
        for angle in angles
    ]
    
    return np.array(points)



# Parameters for the 3D circle
# radius = 1.0  # Radius of the circle
centerlist = [  [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.75, 0.0, 0.0],
                [0.0, 0.75, 0.0],
                [-0.75, 0.0, 0.0],
                [0.0, -0.75, 0.0]]
radius_list = [0.5, 1.0, 0.25, 0.25, 0.25, 0.25]

normal_list = [[0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0]]
  
# normal = np.array([0.0, 0.0, 1.0])  # Normal vector to the circle's plane
# Generate points
num_list = np.array([2,2, 1, 1, 1,1])



for pn in range(1,11):
    point_list = []
    target_sample_num = 1000 * pn
    base_num = target_sample_num / np.sum(num_list)
    pnum_list = num_list * base_num - np.array([2, 2, -1, -1, -1, -1])
    # pnum_list = num_list * base_num - np.array([2, 2, -1, -1, -1, -1])
    for i in range(6):
        points = sample_circle_3d(centerlist[i], normal_list[i], radius_list[i], num_points=pnum_list[i])
        point_list.append(points)

    # Save the sampled points to an .xyz file
    output_file = f"./data/torus_wire/torus_wire_{pn}k.xyz"

    # Write points to the file
    with open(output_file, "w") as f:
        for points in point_list:
            for point in points:
                f.write(f"{point[0]} {point[1]} {point[2]}\n")