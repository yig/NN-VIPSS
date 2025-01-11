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

data_dir = "./data/kitten_small_vipss"
# vipss_folder = 'kitten_small_vipss'
# lv_folder = 'kitten_small_lv'
lv_stats_dir = "./out/kitten_timing/lv2"
lv_stats_dir0 = "./out/kitten_timing/lv0"
lv_stats_dir1 = "./out/kitten_timing/lv1"

file_ids = range(1, 9) 
file_ids = list(map(lambda x: (x * 500), file_ids))
file_names = list(map(lambda x: str(x), file_ids))

print(file_names)

lv_vals_dict =  get_stats_all_folder(lv_stats_dir, file_names)
lv_vals_dict0 =  get_stats_all_folder(lv_stats_dir0, file_names)
lv_vals_dict1 =  get_stats_all_folder(lv_stats_dir1, file_names)


vipss_keys = ['global HRBF coefficient', 'adgrid generation']
vipss_keys_map = {'global HRBF coefficient': 'Global HRBF coefficient calculation', 'adgrid generation':'Global HRBF adgrid generation'}
# lv_keys = ['init normal', 'construct Hmat', 'optimization', 'generate voro', 'opt num', 'opt per iter time']

# lv_keys = ['init normal', 'construct Hmat', 'optimization', 'generate voro']
lv_keys = ['init normal', 'construct Hmat', 'optimization', 'tetgen triangulation']
lv_keys_map = {'init normal' : 'LV init normal', 'construct Hmat' : 'LV build H', 'optimization' : 'LV optimization', 'tetgen triangulation' : 'LV tetgen triangulation'}


plt.figure(figsize=(16, 12))

all_time = [(a + b + c) for a, b, c in zip( lv_vals_dict[lv_keys[0]],  lv_vals_dict[lv_keys[1]],  lv_vals_dict[lv_keys[2]])]
all_time0 = [(a + b + c) for a, b, c in zip( lv_vals_dict0[lv_keys[0]],  lv_vals_dict0[lv_keys[1]],  lv_vals_dict0[lv_keys[2]])]
all_time1 = [(a + b + c) for a, b, c in zip( lv_vals_dict1[lv_keys[0]],  lv_vals_dict1[lv_keys[1]],  lv_vals_dict1[lv_keys[2]])]
all_time_ave =  [(a + b + c)/3.0 for a, b, c in zip(all_time, all_time0, all_time1)]

plt.plot(file_ids, all_time_ave, marker='o', label='LOCAL VIPSS Normal Estimation')

k1 = 'tetgen triangulation'
k2 = 'init normal'
lv_vals_dict0[k2] = [(a - b) for a, b in zip(lv_vals_dict0[k2],  lv_vals_dict0[k1])]
lv_vals_dict1[k2] = [(a - b) for a, b in zip(lv_vals_dict1[k2],  lv_vals_dict1[k1])]
lv_vals_dict[k2] = [(a - b) for a, b in zip(lv_vals_dict[k2],  lv_vals_dict[k1])]

for key in lv_keys:
    means = [(a + b + c) / 3 for a, b, c in zip( lv_vals_dict[key],  lv_vals_dict0[key],  lv_vals_dict1[key])]
    plt.plot(file_ids, means, marker='o', linestyle='--', label=lv_keys_map[key])

# Add labels, title, and legend
plt.xlabel("Input Point Number")
plt.ylabel("Local VIPSS time(s)")
plt.title("Local Vipss Running Time - Hard Constraints")
plt.legend()
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
