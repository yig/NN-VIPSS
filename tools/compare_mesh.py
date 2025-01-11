import numpy as np
import trimesh
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

def compute_vertex_chamfer_distance(mesh1, mesh2):
    """
    Compute Chamfer Distance between the vertices of two meshes.

    Args:
        mesh1 (trimesh.Trimesh): The first mesh.
        mesh2 (trimesh.Trimesh): The second mesh.

    Returns:
        tuple: Chamfer distance (float), per-vertex distances for mesh1.
    """
    # Extract vertices
    vertices1 = mesh1.vertices
    vertices2 = mesh2.vertices

    # Build KDTree for mesh2 vertices
    tree2 = cKDTree(vertices2)

    # Compute distances from vertices of mesh1 to mesh2
    distances1, _ = tree2.query(vertices1, k=1)

    # Chamfer distance (one-directional)
    chamfer = np.mean(distances1**2)
    return chamfer, distances1

def colorize_mesh(mesh, vertex_distances, colormap='viridis'):
    """
    Apply colors to the mesh vertices based on distances using a colormap.

    Args:
        mesh (trimesh.Trimesh): The input mesh.
        vertex_distances (np.ndarray): Per-vertex distances.
        colormap (str): Name of the matplotlib colormap.

    Returns:
        trimesh.Trimesh: Mesh with colored vertices.
    """
    # Normalize distances for coloring
    normalized = (vertex_distances - np.min(vertex_distances)) / (np.max(vertex_distances) - np.min(vertex_distances))
    cmap = plt.get_cmap(colormap)
    colors = cmap(normalized)[:, :3]  # RGB colors (ignore alpha)
    mesh.visual.vertex_colors = (colors * 255).astype(np.uint8)
    return mesh


# Example usage
if __name__ == "__main__":
    # Load two meshes
    path2 = './out/adgrid/kitten_sample_ghrbf/kitten_100u_2000u_combine_mesh.obj'
    path1 = './out/adgrid/kitten_sample_lv/kitten_100u_2000u_combine_mesh.obj'

    mesh1 = trimesh.load(path1)  # Replace with your mesh file path
    mesh2 = trimesh.load(path2)  # Replace with your mesh file path
    
    # Compute Chamfer Distance and per-vertex distances
    chamfer, vertex_distances1 = compute_vertex_chamfer_distance(mesh1, mesh2)
    print(f"Chamfer Distance: {chamfer}")

    # Colorize the first mesh based on distances
    colored_mesh1 = colorize_mesh(mesh1, vertex_distances1)

    # Save and visualize the result
    colored_mesh1.export("mesh1_colored.ply")
    print("Colored mesh saved as 'mesh1_colored.ply'")
    # colored_mesh1.show()
