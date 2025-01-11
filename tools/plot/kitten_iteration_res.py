import os
import csv
import matplotlib.pyplot as plt
import numpy as np


res_vect = []
# stats_path = 'out/kitten_timing/noise_lamda/kitten_10k_2_res_0.01.txt'
stats_path = 'out/test4/kitten_10k_3_res_80000.txt'
count = 0
with open(stats_path, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if the file has one
    for row in reader:
        # Take the first column as the key
        # print(row[0])
        count += 1
        if float(row[0]) > 50000: 
            continue
            print(float(row[0]))
        if count > 20000 : 
            if float(row[0]) > 1000: 
                continue

        res_vect.append(float(row[0]))    


res_vect = res_vect[:len(res_vect)]

file_ids = range(0, len(res_vect)) 

key = 'opt num'
plt.figure(figsize=(16, 9))
plt.plot(file_ids, res_vect, label='10k 3% noise lambda 0.01')
    
# plt.xscale('log')
# plt.xticks([ 1, 2, 3, 4, 5], labels=["0", "0.0001", "0.001", "0.01", "0.1"])
# Add labels, title, and legend

plt.tick_params(axis='both', which='major', labelsize=14)  # Font size for major ticks
plt.tick_params(axis='both', which='minor', labelsize=10)  # Font size for minor ticks (if any)

plt.xlabel("Iter Count", fontsize=14)
plt.ylabel("LOCAL VIPSS Energy", fontsize=14)
plt.title("Local Vipss Optimization Energy S init as Initial Cluster VIPSS Dist Vals - Hard Constraints", fontsize=16)
plt.legend()
plt.grid(True)
# Show the plot
# plt.show()

save_path = './out/figure/sv_res_kitten_10k_3_0.01_s_init_v.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')


