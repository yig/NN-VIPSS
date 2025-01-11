import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import copy



def get_stats_all_folder(folder_dir, file_names):
    vipss_vals_dict = {}
    for name in file_names:
        stats_path = os.path.join(folder_dir,  name + '_time_stats.txt') 
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
        stats_path = os.path.join(folder_dir, name + '_time_stats.txt') 
        with open(stats_path, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header if the file has one
            for row in reader:
                # Take the first column as the key
                key = row[0]
                vipss_vals_dict[key].append(float(row[1]))
    return vipss_vals_dict

data_dir = "./data/kitten_large_n"
# vipss_folder = 'kitten_small_vipss'
# lv_folder = 'kitten_small_lv'
# lv_stats_dir0 = "./out/kitten_timing_l"
lv_stats_dir0 = "./out/torus_wire_soft"
# lv_stats_dir1 = "./out/kitten_timing/lv_large1"
# lv_stats_dir2 = "./out/kitten_timing/lv_large2"

# lv_stats_dir0 = "./out/kitten_timing/lv_large_soft0"
# lv_stats_dir1 = "./out/kitten_timing/lv_large_soft1"
# lv_stats_dir2 = "./out/kitten_timing/lv_large_soft2"
# lv_stats_dir0 = "./out/kitten_timing/lv0"
# lv_stats_dir1 = "./out/kitten_timing/lv1"
file_name_root = 'torus_wire_' 

file_ids = range(1, 11) 
file_ids = list(map(lambda x: (x), file_ids))
file_names = list(map(lambda x: file_name_root + str(x) + 'k', file_ids))

print(file_names)

lv_vals_dict0 =  get_stats_all_folder(lv_stats_dir0, file_names)
# lv_vals_dict1 =  get_stats_all_folder(lv_stats_dir1, file_names)
# lv_vals_dict2 =  get_stats_all_folder(lv_stats_dir2, file_names)


# vipss_keys = ['global HRBF coefficient', 'adgrid generation']
# vipss_keys_map = {'global HRBF coefficient': 'Global HRBF coefficient calculation', 'adgrid generation':'Global HRBF adgrid generation'}
# lv_keys = ['init normal', 'construct Hmat', 'optimization', 'generate voro', 'opt num', 'opt per iter time']

# lv_keys = ['init normal', 'construct Hmat', 'optimization', 'generate voro']
# lv_keys = ['init normal', 'construct Hmat', 'optimization', 'tetgen triangulation']
# lv_keys = ['optimization']
# lv_keys_map = {'init normal' : 'LV init normal', 'construct Hmat' : 'LV build H', 'optimization' : 'LV optimization', 'tetgen triangulation' : 'LV tetgen triangulation'}

lv_keys = ['optimization', 'opt num']
# 
# lv_keys = ['ave cluster size', 'cluster size std deviation']
lv_keys_map = {'init normal' : 'init normal partial vipss', 'construct Hmat' : 'build H mat', 
                'optimization' : 'optimization with soft constraints', 
                'tetgen triangulation' : 'tetgen triangulation', 
                'opt num': 'optimization Iteration Number',
                'ave cluster size' : 'ave cluster size',
                'cluster size std deviation' : 'cluster size std deviation'}

plt.figure(figsize=(16, 10))

# all_time = [(a + b + c) for a, b, c in zip( lv_vals_dict2[lv_keys[0]],  lv_vals_dict2[lv_keys[1]],  lv_vals_dict2[lv_keys[2]])]
# all_time0 = [(a + b + c) for a, b, c in zip( lv_vals_dict0[lv_keys[0]],  lv_vals_dict0[lv_keys[1]],  lv_vals_dict0[lv_keys[2]])]
# all_time1 = [(a + b + c) for a, b, c in zip( lv_vals_dict1[lv_keys[0]],  lv_vals_dict1[lv_keys[1]],  lv_vals_dict1[lv_keys[2]])]
# all_time_ave =  [(a + b + c)/3.0 for a, b, c in zip(all_time, all_time0, all_time1)]

# plt.plot(file_ids, all_time0, marker='o', label='Normal Estimation Total Time')
k1 = 'tetgen triangulation'
k2 = 'init normal'
# lv_vals_dict0[k2] = [(a - b) for a, b in zip(lv_vals_dict0[k2],  lv_vals_dict0[k1])]
# lv_vals_dict1[k2] = [(a - b) for a, b in zip(lv_vals_dict1[k2],  lv_vals_dict1[k1])]
# lv_vals_dict2[k2] = [(a - b) for a, b in zip(lv_vals_dict2[k2],  lv_vals_dict2[k1])]

ave_iter =  [(a /b) for a, b in zip(lv_vals_dict0[lv_keys[0]],  lv_vals_dict0[lv_keys[1]])]

plt.plot(file_ids, ave_iter, marker='o', linestyle='--', label='time per iteration')
# ave_list = []
# for val in lv_vals_dict0[lv_keys[0]]:
#     if val > 0:
#         ave_list.append(val)

# plt.plot(file_ids, ave_list, marker='o', linestyle='--', label='ave cluster size')
# plt.plot(file_ids, lv_vals_dict0[lv_keys[1]], marker='o', linestyle='--', label=lv_keys_map[lv_keys[1]])
for key in lv_keys:

    # means = [(a + b + c) / 3 for a, b, c in zip( lv_vals_dict2[key],  lv_vals_dict0[key],  lv_vals_dict1[key])]
    plt.plot(file_ids, lv_vals_dict0[key], marker='o', linestyle='--', label=lv_keys_map[key])
    print(lv_vals_dict0[key])


plt.tick_params(axis='both', which='major', labelsize=14)  # Font size for major ticks
plt.tick_params(axis='both', which='minor', labelsize=10)  # Font size for minor ticks (if any)
# Add labels, title, and legend
plt.xlabel("Input Point Number(k)", fontsize=14)
plt.ylabel(" Time (s) / Count ", fontsize=14)
plt.title("Torus Wire Scalable VIPSS Optimization Timing -- soft constraints", fontsize=18)
# plt.ylabel(" Number ", fontsize=14)
# plt.title("Torus Wire Scalable VIPSS Stats -- soft constraints", fontsize=18)

plt.legend()

# Display grid for better readability
plt.grid(True)

save_path = './out/figure/torus_wire_opt_num_timing.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')

