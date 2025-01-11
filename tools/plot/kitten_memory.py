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



# indata_dir = './data/kitten_large_n'
data_dir = './out/memory'

file_ids = range(1, 11) 
file_ids = list(map(lambda x: (x * 100), file_ids))
file_names = list(map(lambda x: 'kitten_' + str(x) + 'k_time_stats.txt', file_ids))
key = 'max memory'
time_dict = get_stats_all_folder(data_dir, file_names) 
print((time_dict[key]))
plt.figure(figsize=(16, 10))

plt.plot(file_ids, time_dict[key], marker='o', label='Scalable Vipss Memory Usage')
    
# plt.xscale('log')
# plt.xticks([ 1, 2, 3, 4, 5], labels=["0", "0.0001", "0.001", "0.01", "0.1"])
# Add labels, title, and legend
plt.tick_params(axis='both', which='major', labelsize=14)  # Font size for major ticks
plt.tick_params(axis='both', which='minor', labelsize=10)  # Font size for minor ticks (if any)

plt.xlabel("Input Point Number(k)", fontsize=18)
plt.ylabel("Memory Size(G)",fontsize=18 )
plt.title("Scalable VIPSS Memory Usage", fontsize=20)
plt.legend()
plt.grid(True)
# Show the plot
# plt.show()
save_path = './out/figure/kitten_memory.png'
plt.savefig(save_path)

