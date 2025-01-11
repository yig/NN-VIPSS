import numpy as np
from scipy.spatial import distance_matrix
import os
import kdtree

def load_obj(file_path):
    vertices = []
    faces = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('v '):
                parts = line.strip().split()
                vertex = [float(parts[1]), float(parts[2]), float(parts[3])]
                vertices.append(vertex)
            elif line.startswith('f '):
                faces.append(line.strip())
    
    return np.array(vertices), faces

def save_obj(file_path, vertices, faces):
    with open(file_path, 'w') as f:
        for v in vertices:
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for face in faces:
            f.write(face + '\n')

def save_xyz(file_path, vertices):
    with open(file_path, 'w') as xyz_file:
        for i in range(len(vertices)):
            vertex = vertices[i]
            # Write vertex and normal data to XYZ file
            xyz_file.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")


def remove_close_vertices(vertices, threshold):
    # dist_matrix = distance_matrix(vertices, vertices)
    # to_remove = set()
    
    # for i in range(len(vertices)):
    #     for j in range(i + 1, len(vertices)):
    #         if dist_matrix[i, j] < threshold:
    #             to_remove.add(j)

    # filtered_vertices = np.array([v for i, v in enumerate(vertices) if i not in to_remove])
    # np.random.shuffle(filtered_vertices)
    emptyTree = kdtree.create(dimensions=3)
    results = []
    for vertex in vertices:
        dist = emptyTree.search_knn(vertex, 1)
        if len(dist) > 0:
            if dist[0][1] < threshold:
                continue
        emptyTree.add(vertex)
        results.append(vertex)
    # print(results)
    np.random.shuffle(results)
    return results

# Paths to input and output files
# input_obj_file = r'c:\Users\xiaji\Documents\projects\3d-sketches-curated-dataset\scaffolds3d_dog.obj'
# output_obj_file = 'scaffolds3d_dog.obj'
distance_threshold = 0.000001  # Set your threshold for "close" points

input_dir = r'c:\Users\xiaji\Documents\projects\sketches'
filenames = os.listdir(input_dir)
output_dir = r'c:\Users\xiaji\Documents\projects\sketches_refine'

for filename in filenames:
    input_obj_file = input_dir + "/" + filename
    # Load vertices and faces from the obj file
    vertices, faces = load_obj(input_obj_file)
    # Remove vertices that are closer than the threshold
    # filtered_vertices = remove_close_vertices(vertices, distance_threshold)
    # noise = np.random.normal(0, 0.001, vertices.shape)
    # filtered_vertices = vertices + noise
    filtered_vertices = remove_close_vertices(vertices, distance_threshold)
    output_xyz_file = output_dir + "/" + filename.split('.')[0] + '.xyz'
    # Save the filtered obj file
    save_xyz(output_xyz_file, filtered_vertices)

    print(f"Filtered {len(vertices) - len(filtered_vertices)} close vertices.")
