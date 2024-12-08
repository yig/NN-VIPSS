import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import copy

data_dir = "../../data/planck_sample"

octree_depth = ['octree_depth_0', 'octree_depth_3', 'octree_depth_4', 'octree_depth_5']
stats_dir = "../../out/octree_depth_0"
out_root = '../../out'

file_ids = range(1, 11) 
file_ids = list(map(lambda x: (x * 10), file_ids))
file_names = list(map(lambda x: str(x) + 'k', file_ids))

print(file_names)

keys = []
all_vals_dict = {}

for name in file_names:
    stats_path = os.path.join(stats_dir, 'planck_' + name + '_time_stats.txt') 
    with open(stats_path, mode='r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header if the file has one
        for row in reader:
            # Take the first column as the key
            key = row[0]
            keys.append(key)
            all_vals_dict[key] = []        
        print(keys)
        break

all_level_data = {}
for octree_level in octree_depth:
    all_level_data[octree_level] =  copy.deepcopy(all_vals_dict)
    octree_out_dir = os.path.join(out_root, octree_level)
    # file_names = os.listdir(octree_out_dir)
    for name in file_names:

        stats_path = os.path.join(octree_out_dir, 'planck_' + name + '_time_stats.txt') 
        with open(stats_path, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header if the file has one
            for row in reader:
                # Take the first column as the key
                key = row[0]
                all_level_data[octree_level][key].append(float(row[1]))
        # print(all_level_data[octree_level])

non_display_keys = ['ave evaluate nn num', 'ave cluster size','octree dummy pt num','tetgen triangulation', 'init normal', 'construct Hmat', 'optimization', 'opt num', 'opt per iter time', 'generate voro', 'build nn HRBF', 'evaluate count']


plt.figure(figsize=(16, 12))
    
for key in keys:
    if key in non_display_keys:
        continue
    for octree_level in octree_depth:
        # print(all_level_data[octree_level][key])
        plt.plot(file_ids, all_level_data[octree_level][key], marker='o', label=key+' '+octree_level)


# Add labels, title, and legend
plt.xlabel("Input Point Number(k)")
plt.ylabel("nn surface time(s)")
plt.title("Local vipss NN surface time with different level of octree")
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
