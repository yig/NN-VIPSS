import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import copy



def get_stats_all_folder_vipss(folder_dir, file_names):
    vipss_vals_dict = {}
    for name in file_names:
        stats_path = os.path.join(folder_dir, 'kitten_' + name + '_time.txt') 
        with open(stats_path, mode='r') as file:
            reader = csv.reader(file, delimiter=':')
            next(reader)  # Skip the header if the file has one
            for row in reader:
                # Take the first column as the key
                key = row[0]
                # keys.append(key)
                vipss_vals_dict[key] = []      
            # print(keys)
            break
    for name in file_names:
        stats_path = os.path.join(folder_dir, 'kitten_' + name + '_time.txt') 
        with open(stats_path, mode='r') as file:
            reader = csv.reader(file, delimiter=':')
            next(reader)  # Skip the header if the file has one
            for row in reader:
                # Take the first column as the key
                key = row[0]
                val_s = row[1].replace('s', '')
                vipss_vals_dict[key].append(float(val_s))
    return vipss_vals_dict

data_dir = "./data/kitten_small_vipss"
# vipss_folder = 'kitten_small_vipss'
# lv_folder = 'kitten_small_lv'
lv_stats_dir = "./out/kitten_timing/vipss0"

file_ids = range(1, 9) 
file_ids = list(map(lambda x: (x * 500), file_ids))
file_names = list(map(lambda x: str(x), file_ids))

print(file_names)

vipss_vals_dict = get_stats_all_folder_vipss(lv_stats_dir, file_names)
# lv_vals_dict =  get_stats_all_folder(lv_stats_dir, file_names)


vipss_keys = ['setup_time (Compute H)', 'init_time (Optimize g/Eigen)', 'solve_time (Optimize g/LBFGS)']
vipss_keys_map = {'setup_time (Compute H)': 'VIPSS Compute H', 'init_time (Optimize g/Eigen)': 'VIPSS Init Normals', 'solve_time (Optimize g/LBFGS)' : 'VIPSS Optimization'}
# lv_keys = ['init normal', 'construct Hmat', 'optimization', 'generate voro', 'opt num', 'opt per iter time']

plt.figure(figsize=(16, 12))

all_time = [(a + b + c) for a, b, c in zip( vipss_vals_dict[vipss_keys[0]],  vipss_vals_dict[vipss_keys[1]],  vipss_vals_dict[vipss_keys[2]])]
plt.plot(file_ids, all_time, marker='o', label='VIPSS Normal Estimation')
for key in vipss_keys:
    plt.plot(file_ids, vipss_vals_dict[key], marker='o', linestyle='--', label=vipss_keys_map[key])


# Add labels, title, and legend
plt.xlabel("Input Point Number")
plt.ylabel("VIPSS time(s)")
plt.title("Vipss Running Time")
plt.legend()

# Display grid for better readability
plt.grid(True)

# Show the plot
plt.show()
