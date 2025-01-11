import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import copy



def get_stats_all_folder(folder_dir, file_names):
    vipss_vals_dict = {}
    for name in file_names:
        stats_path = os.path.join(folder_dir, 'kitten_' + name + '_time_stats.txt') 
        with open(stats_path, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header if the file has one
            for row in reader:
                # Take the first column as the key
                key = row[0]
                # keys.append(key)
                vipss_vals_dict[key] = []      
            # print(keys)
            break

    for name in file_names:
        stats_path = os.path.join(folder_dir, 'kitten_' + name + '_time_stats.txt') 
        with open(stats_path, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header if the file has one
            for row in reader:
                # Take the first column as the key
                key = row[0]
                vipss_vals_dict[key].append(float(row[1]))
    return vipss_vals_dict


# octree_depth = ['kitten_origin_mp', 'kitten_octree4_mp', 'kitten_octree5_mp', 'kitten_octree6_mp']
# octree_depth = ['kitten_origin', 'kitten_octree4', 'kitten_octree5', 'kitten_octree6']
# vipss_folder = 'kitten_small_ghrbf'
lv_folder1 = 'kitten_k_lv_001'
lv_folder2 = 'kitten_k_lv_002'
lv_folder3 = 'kitten_k_lv_003'
# lv_folder1 = 'kitten_k_lv_old01'
# lv_folder2 = 'kitten_k_lv_old02'
# lv_folder3 = 'kitten_k_lv_old03'
stats_dir = "./out/adgrid"
out_root = './out'

file_ids = range(1, 11) 
file_ids = list(map(lambda x: (x * 10), file_ids))
file_names = list(map(lambda x: str(x) + 'k', file_ids))
lv_dir1 = os.path.join(stats_dir, lv_folder1)
lv_dir2 = os.path.join(stats_dir, lv_folder2)
lv_dir3 = os.path.join(stats_dir, lv_folder3)

lv_vals_dict1 =  get_stats_all_folder(lv_dir1, file_names)
lv_vals_dict2 =  get_stats_all_folder(lv_dir2, file_names)
lv_vals_dict3 =  get_stats_all_folder(lv_dir3, file_names)

lv_surface_tracker_dir1 = "./out/surface_track"
lv_surface_tracker_dict1 =  get_stats_all_folder(lv_surface_tracker_dir1, file_names)


vipss_keys = ['global HRBF coefficient', 'adgrid generation']
vipss_keys_map = {'global HRBF coefficient': 'Global HRBF coefficient calculation', 'adgrid generation':'Global HRBF adgrid generation'}
lv_keys = ['build nn HRBF', 'adgrid generation']
lv_keys_map = {'build nn HRBF': 'Local Vipss Coefficient Calculation', 'adgrid generation':'Local Vipss Evaluation'}

lv_st_keys = ['build nn HRBF', 'NN surface']
lv_st_keys_map = {'build nn HRBF': 'Local Vipss Coefficient Calculation', 'NN surface':'Local Vipss Surface Tracker Evaluation'}

lv_st_surface_time = [a + b  for a, b in zip(lv_surface_tracker_dict1[lv_st_keys[0]], lv_surface_tracker_dict1[lv_st_keys[1]])]




plt.figure(figsize=(20, 12))

nn_time = [(a + b + c) / 3.0 for a, b, c in zip(lv_vals_dict1[lv_keys[0]], lv_vals_dict2[lv_keys[0]], lv_vals_dict3[lv_keys[0]])]

adgrid_time = [(a + b + c) / 3.0 for a, b, c in zip(lv_vals_dict1[lv_keys[1]], lv_vals_dict2[lv_keys[1]], lv_vals_dict3[lv_keys[1]])]


lv_surface_time = [a + b  for a, b in zip(nn_time, adgrid_time)]


plt.plot(file_ids, lv_surface_time, marker='o', label='LV Adgrid Surface time')

# plt.plot(file_ids, nn_time, marker='o', linestyle='--', label='build nn HRBF')
# plt.plot(file_ids, adgrid_time, marker='o', linestyle='--', label='Local Vipss Evaluation')

plt.plot(file_ids, lv_st_surface_time, marker='o', label='LV Surface Tracker time')
# plt.plot(file_ids, lv_surface_tracker_dict1[lv_st_keys[0]], marker='o', linestyle='--',label='LV Surface Tracker Build NN time')
# plt.plot(file_ids, lv_surface_tracker_dict1[lv_st_keys[1]], marker='o', linestyle='--', label='LV Surface Tracker Evaluation time')

# for key in lv_keys:
#     plt.plot(file_ids, lv_vals_dict[key], marker='o', linestyle='--', label=lv_keys_map[key])
poisson_path =  './out/adgrid/kitten_poisson/time.txt'
with open(poisson_path, 'r') as file:
    loaded_list = [float(line.strip()) for line in file]
    plt.plot(file_ids, loaded_list, marker='o', label='Screen Poisson surface time')
    print(loaded_list)

# Add labels, title, and legend
plt.xlabel("Input Point Number(k)")
plt.ylabel("surface time(s)")
plt.title("Surface time comparison between Local Vipss and Screened Poisson")
plt.legend()

# Display grid for better readability
plt.grid(True)

# Show the plot
# plt.show()
fig_save_path = './out/figure/surface_time_lv_and_screend_poisson_only.png'
plt.savefig(fig_save_path)
        

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
