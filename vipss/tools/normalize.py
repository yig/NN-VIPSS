import numpy as np

# Function to read an XYZ file
def read_xyz(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        # num_atoms = int(lines[0].strip())  # First line contains the number of atoms
        # comment = lines[1].strip()         # Second line is a comment
        atoms = []
        coordinates = []

        for line in lines[2:]:
            parts = line.split()
            # atoms.append(parts[0])
            coordinates.append([float(parts[0]), float(parts[1]), float(parts[2])])

        return atoms, np.array(coordinates)

# Function to normalize the coordinates
def normalize_coordinates(coords):
    # Translate coordinates to have mean at origin
    mean = np.mean(coords, axis=0)
    centered_coords = coords - mean
    
    # Scale the coordinates so that the maximum absolute value in any dimension is 1
    max_val = np.max(np.abs(centered_coords))
    normalized_coords = centered_coords / max_val

    return normalized_coords

# Function to write the normalized coordinates back to an XYZ file
def write_xyz(filename,  coordinates):
    with open(filename, 'w') as file:
        for coord in  coordinates:
            file.write(f"{coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

# Example usage
input_file = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\elbow_scan\6007_9_9.xyz'
output_file = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\elbow_scan\6007_9_9_normalized.xyz'

# Read the original XYZ file
atoms, coords = read_xyz(input_file)

# Normalize the coordinates
normalized_coords = normalize_coordinates(coords)

# Write the normalized data back to a new XYZ file
write_xyz(output_file,  normalized_coords)

print(f"Normalized XYZ file written to {output_file}")