import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import copy
import trimesh

lv_folder = 'kitten_k_lv'
stats_dir = "./out/adgrid"
out_root = './out'

file_ids = range(1, 11) 
file_ids = list(map(lambda x: (x * 10), file_ids))
file_names = list(map(lambda x: str(x) + 'k', file_ids))
# vipss_dir = os.path.join(stats_dir, vipss_folder)
lv_dir = os.path.join(stats_dir, lv_folder)
thresholds = [0.0012, 0.00145, 0.00163, 0.00177, 0.0019, 0.00204, 0.00212, 0.00222, 0.00229, 0.00236]

lv_v_nums = []
lv_f_nums = []

for i in range(len(file_names)):
    cur_name = file_names[i]
    t_str = str(thresholds[i])
    mesh_path = os.path.join(lv_dir, 'kitten_' + cur_name + '_' + t_str +  '_mesh.obj') 
    mesh = trimesh.load(mesh_path)
    lv_v_nums.append(mesh.vertices.shape[0])
    lv_f_nums.append(mesh.faces.shape[0])

p_v_nums = []
p_f_nums = []

poisson_dir = "./out/adgrid/kitten_poisson"
for i in range(len(file_names)):
    cur_name = file_names[i]
    t_str = str(thresholds[i])
    mesh_path = os.path.join(poisson_dir, 'kitten_' + cur_name + '.ply') 
    mesh = trimesh.load(mesh_path)
    p_v_nums.append(mesh.vertices.shape[0])
    p_f_nums.append(mesh.faces.shape[0])




plt.figure(figsize=(16, 12))
# plt.plot(file_ids, lv_v_nums, marker='o', label='adgrid vertice number')
# plt.plot(file_ids, lv_f_nums, marker='o', label='adgrid face number')

# plt.plot(file_ids, p_v_nums, marker='o', label='poisson vertice number')
# plt.plot(file_ids, p_f_nums, marker='o', label='poisson face number')

# evaluate_times = all_level_data[octree_depth[0]]['evaluate count'] 
# print('evalte times: ', evaluate_times)
# plt.plot(file_ids, vipss_time_list, marker='o', label='vipss time')

# vipss_surface_time = [a + b for a, b in zip(vipss_vals_dict[vipss_keys[0]], vipss_vals_dict[vipss_keys[1]])]
# plt.plot(file_ids, vipss_surface_time, marker='o', label='Global HRBF surface time')
# for key in vipss_keys:
#     plt.plot(file_ids, vipss_vals_dict[key], marker='o', linestyle='--', label=vipss_keys_map[key])


thresholds = [0.0012, 0.00145, 0.00163, 0.00177, 0.0019, 0.00204, 0.00212, 0.00222, 0.00229, 0.00236]
plt.plot(file_ids, thresholds, marker='o', label='adgrid threshold')


# Add labels, title, and legend
plt.xlabel("Input Point Number(k)")
# plt.ylabel("number")
plt.ylabel("threshold")
# plt.title("Adaptive grid mesh and poisson mesh vertex and face number")
plt.title("Adaptive grid threshold")
plt.legend()

# Display grid for better readability
plt.grid(True)

# Show the plot
plt.show()

        

# # Data
# x = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 9000]
# y1 = [float(value) for value in [' 0.0052255', ' 0.0031252', ' 0.0088357', ' 0.0068443', ' 0.0182718', ' 0.0108671', ' 0.0120147', ' 0.027068', ' 0.0182355', ' 0.0315393', ' 0.0302044', ' 0.0298475', ' 0.0326081']]
# y2 = [float(value) for value in [' 0.0742291', ' 0.0745208', ' 0.104073', ' 0.101811', ' 0.153238', ' 0.149737', ' 0.17805', ' 0.264459', ' 0.234297', ' 0.299448', ' 0.363024', ' 0.420034', ' 0.455142']]
# y3 = [float(value) for value in [' 0.0618171', ' 0.106564', ' 0.165052', ' 0.193553', ' 0.24084', ' 0.277491', ' 0.298805', ' 0.361136', ' 0.432766', ' 0.502843', ' 0.626014', ' 0.710662', ' 0.851845']]

# # Calculate global Y-axis range
# all_y_values = y1 + y2 + y3  # Combine all Y data
# y_min = int(min(all_y_values) * 10) / 10  # Round down for consistency
# y_max = int(max(all_y_values) * 10 + 1) / 10  # Round up for consistency
# y_step = 0.1  # Define step size for Y-ticks

# # Create the plot
# plt.figure(figsize=(10, 6))
# plt.plot(x, y1, marker='o', label="y1 (Series 1)")
# plt.plot(x, y2, marker='s', label="y2 (Series 2)")
# plt.plot(x, y3, marker='^', label="y3 (Series 3)")

# # Set consistent Y-axis for all plots
# plt.ylim(y_min, y_max)
# plt.yticks(np.arange(y_min, y_max + y_step, y_step))

# # Add labels, title, and legend
# plt.xlabel("X-axis (Input X values)")
# plt.ylabel("Y-axis")
# plt.title("Line Graph with Consistent Y-Axis")

# # Place the legend outside the graph
# plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')

# # Display grid for better readability
# plt.grid(True)

# # Adjust layout to prevent clipping of the legend
# plt.tight_layout()

# # Save the plot to a file (e.g., PNG format)
# plt.savefig('line_graph.png', format='png')  # Specify the filename and format

# # Show the plot
# plt.show()
