import numpy as np

def load_xyz(file_path):
    """
    Load a point cloud from an .xyz file.
    
    Args:
        file_path (str): Path to the .xyz file.
    
    Returns:
        np.ndarray: A numpy array of shape (N, 3) representing point coordinates.
    """
    return np.loadtxt(file_path)

def resize_coordinates(points, scale_factor):
    """
    Resize (scale) the coordinates of the point cloud.
    
    Args:
        points (np.ndarray): Array of shape (N, 3) with point coordinates.
        scale_factor (float): Scaling factor to resize coordinates.
    
    Returns:
        np.ndarray: Scaled point cloud.
    """
    return points * scale_factor

def save_xyz(file_path, points):
    """
    Save a point cloud to an .xyz file.
    
    Args:
        file_path (str): Path to save the .xyz file.
        points (np.ndarray): Array of shape (N, 3) with point coordinates.
    """
    np.savetxt(file_path, points, fmt='%.6f')

# Example usage
input_file = "./data/kitten_sample/kitten_1000_nu.xyz"   # Path to the input file
output_file = "./data/kitten_sample_resize/kitten_1000_nu.xyz" # Path to the output file
scale_factor = 2.0         # Scale factor for resizing coordinates

# Load, resize, and save the point cloud
point_cloud = load_xyz(input_file)
scaled_point_cloud = resize_coordinates(point_cloud, scale_factor)
save_xyz(output_file, scaled_point_cloud)

print(f"Point cloud resized and saved to {output_file}")
