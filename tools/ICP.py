import open3d as o3d
import numpy as np

# Load the source and target point clouds from .xyz files
def load_point_cloud(file_path):
    # Load point cloud from .xyz file
    point_cloud = np.loadtxt(file_path)
    point_cloud_o3d = o3d.geometry.PointCloud()
    point_cloud_o3d.points = o3d.utility.Vector3dVector(point_cloud[:, :3])  # Assuming .xyz file contains x, y, z columns
    return point_cloud_o3d

# Load the source and target point clouds
sampledPtfile = r'c:\Users\xiaji\Documents\projects\3D_pointcloud_dataset\contours\sparse\curve.xyz'
infile = r'c:\Users\xiaji\Downloads\Livewire_demo\Livewire_demo\test data\hand_normal.xyz'
source = load_point_cloud(sampledPtfile)
target = load_point_cloud(infile)

# Preprocessing (optional): downsample the point clouds for faster processing
source_down = source.voxel_down_sample(voxel_size=0.01)
target_down = target.voxel_down_sample(voxel_size=0.01)

# Estimate normals for both point clouds
source_down.estimate_normals(o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))
target_down.estimate_normals(o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))



# Perform ICP registration
# Use point-to-point ICP (alternative: point-to-plane ICP is also available)
threshold = 0.05  # Distance threshold for ICP
initial_transform = np.identity(4)  # Initial guess is identity


# Coarse alignment (large voxel size)
voxel_size_coarse = 0.05
source_coarse = source.voxel_down_sample(voxel_size=voxel_size_coarse)
target_coarse = target.voxel_down_sample(voxel_size=voxel_size_coarse)

# Fine alignment (small voxel size)
voxel_size_fine = 0.01
source_fine = source.voxel_down_sample(voxel_size=voxel_size_fine)
target_fine = target.voxel_down_sample(voxel_size=voxel_size_fine)

# Coarse alignment
icp_coarse_result = o3d.pipelines.registration.registration_icp(
    source_coarse, target_coarse, threshold, initial_transform,
    o3d.pipelines.registration.TransformationEstimationPointToPoint()
)

# # Use the transformation from the coarse alignment for the fine alignment
# icp_fine_result = o3d.pipelines.registration.registration_icp(
#     source_fine, target_fine, threshold, icp_coarse_result.transformation,
#     o3d.pipelines.registration.TransformationEstimationPointToPlane(),
#     o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=2000)
# )

threshold = 0.01

icp_fine_result = o3d.pipelines.registration.registration_icp(
    source_down, target_down, threshold, icp_coarse_result.transformation,
    o3d.pipelines.registration.TransformationEstimationPointToPlane(),
    o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=2000)
)
# # Apply the ICP algorithm
# icp_result = o3d.pipelines.registration.registration_icp(
#     source_down, target_down, threshold, initial_transform,
#     o3d.pipelines.registration.TransformationEstimationPointToPoint()
# )

# Print the transformation matrix
print("Transformation Matrix:")
print(icp_fine_result.transformation)

# Transform the source point cloud using the estimated transformation
source.transform(icp_fine_result.transformation)

# Visualize the aligned point clouds
o3d.visualization.draw_geometries([source, target])

# Save the transformed source point cloud
o3d.io.write_point_cloud("aligned_source_point_cloud.xyz", source)

print("Source point cloud aligned and saved.")
