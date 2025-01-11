import numpy as np

# Function to load an .xyz file
def load_xyz_file(file_path):
    """
    Reads a point cloud from an .xyz file.
    Assumes the file has rows of `x y z` coordinates.
    """
    return np.loadtxt(file_path, usecols=(0, 1, 2))  # Read first 3 columns (x, y, z)

# Function to save a noisy point cloud back to an .xyz file
def save_xyz_file(file_path, point_cloud):
    """
    Saves a 3D point cloud to an .xyz file.
    """
    np.savetxt(file_path, point_cloud, fmt="%.6f")  # Save with 6 decimal precision


noise_percentage = 0.03
np_str = str(int(noise_percentage * 100))
# Example usage
file_path = "./data/kitten_small_vipss/kitten_4000.xyz"  # Replace with your input file path
output_path = f"data/noise_kitten/kitten_4k/kitten_4k_{np_str}.xyz"  # Path to save noisy point cloud



# Load the point cloud
point_cloud = load_xyz_file(file_path)
# Compute the bounding box size (range) along each axis
min_vals = np.min(point_cloud, axis=0)
max_vals = np.max(point_cloud, axis=0)
ranges = max_vals - min_vals  # Range for each dimension
max_len = np.max(ranges) * noise_percentage

noise_range = max_len * np.ones(ranges.shape)

noise = np.random.uniform(-noise_range, noise_range, point_cloud.shape)

# Define the percentage of noise (10% of the range)
# noise_percentage = 0.0025
# std_dev = noise_percentage * ranges
# max_val = max_len * noise_percentage


# print(f'std ranges : {ranges}')
# print(f'std dev : {std_dev}')

# Generate Gaussian noise
# noise = np.random.normal(0, std_dev, point_cloud.shape)
# noise = np.clip(noise,None, max_val)

# Add noise to the point cloud
noisy_point_cloud = point_cloud + noise

# Save the noisy point cloud to a file
save_xyz_file(output_path, noisy_point_cloud)

print(f"Noisy point cloud saved to: {output_path}")
