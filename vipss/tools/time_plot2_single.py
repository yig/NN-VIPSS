import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import copy


#  --------- vipss ave test pt time : 0.000288049
# stats_path = r'../build/dist_time_no_octree_sample_MP.csv'
stats_path = r'../build/dist_time.csv'
origin_path = r'../build/dist_time_no_octree_sample.csv'
octree3_path = r'../build/dist_time_octree3.csv'
octree4_path = r'../build/dist_time_octree_nomp.csv'

origin_time_vals = []
origin_dist_vals = []
count = 0
with open(origin_path, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if the file has one
    for row in reader:
        # Take the first column as the key
        if float(row[1]) > 0.005: continue
        origin_time_vals.append(float(row[1]))
        origin_dist_vals.append(float(row[0]))


o3_time_vals = []
o3_dist_vals = []
count = 0
with open(octree3_path, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if the file has one
    for row in reader:
        # Take the first column as the key
        if float(row[1]) > 0.005: continue
        o3_time_vals.append(float(row[0]))
        o3_dist_vals.append(float(row[1]))

o4_time_vals = []
o4_dist_vals = []
count = 0
with open(octree4_path, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if the file has one
    for row in reader:
        # Take the first column as the key
        if float(row[1]) > 0.005: continue
        o4_time_vals.append(float(row[0]))
        o4_dist_vals.append(float(row[1]))

        # count += 1
        # if count > 1000 : 
        #     break

        # all_level_data[octree_level][key].append(float(row[1]))




plt.figure(figsize=(12, 8))

plt.scatter(origin_time_vals, origin_time_vals, s=1, alpha=0.3)

# plt.scatter(origin_time_vals, origin_time_vals, s=1, alpha=0.3, c='blue' , label='local vipss')

# plt.scatter(o3_time_vals, o3_time_vals, s=1, alpha=0.3, c='red' , label='local vipss')
# plt.scatter(o4_time_vals, o4_time_vals, s=1, alpha=0.3, c='green' , label='local vipss')

# Add labels, title, and legend
plt.xlabel("distance to the surface")
plt.ylabel("distance function evaluate time(s)")
plt.title("Local vipss octree sample and NO OPENMP)")
plt.legend()

# Show the plot
plt.show()