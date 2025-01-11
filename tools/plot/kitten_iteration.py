import os
import csv
import matplotlib.pyplot as plt
import numpy as np

def get_stats_all_folder(folder_dir, file_names):
    vipss_vals_dict = {}
    for name in file_names:
        stats_path = os.path.join(folder_dir, name) 
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
        stats_path = os.path.join(folder_dir, name) 
        with open(stats_path, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header if the file has one
            for row in reader:
                # Take the first column as the key
                key = row[0]
                vipss_vals_dict[key].append(float(row[1]))
    return vipss_vals_dict



indata_dir = './data/kitten_noise'
data_dir = '../../VIPSS_LOCAL/out/kitten_timing/lv_noise_soft'
out_dir = './out/kitten_timing/lv_noise_soft'

file_ids = range(1, 11) 
file_ids = list(map(lambda x: (x * 10), file_ids))
file_names = list(map(lambda x: str(x) + 'k', file_ids))
time_res = []

file_names = ['kitten_10k_0', 'kitten_10k_0.5', 'kitten_10k_1', 'kitten_10k_2', 'kitten_10k_3']
lambda_list = [0, 0.0001, 0.001, 0.01, 0.1]

lambda_list_str = map(lambda x: str(x), lambda_list) 

key = 'opt num'
plt.figure(figsize=(16, 12))

for file_name in file_names:
    iter_nums = []
    cur_names = []
    for l in lambda_list:
        cur_name = file_name + '_' + str(l) +'_time_stats.txt'
        cur_names.append(cur_name)
        # cur_path = os.path.join(data_dir, cur_name)
    time_dict = get_stats_all_folder(out_dir, cur_names) 
    print((time_dict[key]))
    value_list = [1, 2, 3, 4, 5]
    plt.plot(value_list, time_dict[key], marker='o', label=file_name)
    
# plt.xscale('log')
plt.xticks([ 1, 2, 3, 4, 5], labels=["0", "0.0001", "0.001", "0.01", "0.1"])
# Add labels, title, and legend
plt.xlabel("Lambda Value")
plt.ylabel("Iteration Number")
plt.title("Local Vipss Optimization Iteration Number - Soft Constraints")
plt.legend()
plt.grid(True)
# Show the plot
plt.show()

