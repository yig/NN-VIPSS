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
with open(origin_path, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if the file has one
    for row in reader:
        # Take the first column as the key
        if float(row[1]) > 0.005: continue
        origin_dist_vals.append(float(row[0]))
        origin_time_vals.append(float(row[1]))


o3_time_vals = []
o3_dist_vals = []
with open(octree3_path, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if the file has one
    for row in reader:
        # Take the first column as the key
        if float(row[1]) > 0.003: continue
        o3_dist_vals.append(float(row[0]))
        o3_time_vals.append(float(row[1]))

o4_time_vals = []
o4_dist_vals = []
with open(octree4_path, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if the file has one
    for row in reader:
        # Take the first column as the key
        if float(row[1]) > 0.003: continue
        o4_dist_vals.append(float(row[0]))
        o4_time_vals.append(float(row[1]))




plt.figure(figsize=(12, 8))
plt.scatter(origin_dist_vals, origin_time_vals, s=1, alpha=0.1, c='blue', label='local vipss no octree')

plt.scatter(o3_dist_vals, o3_time_vals, s=1, alpha=0.1, c='red', label='octree depth 3')
plt.scatter(o4_dist_vals, o4_time_vals, s=1, alpha=0.1, c='green', label='octree depth 4')


# Add labels, title, and legend
plt.xlabel("distance to the surface")
plt.ylabel("distance function evaluate time(s)")
plt.title("Local vipss octree sample and NO OPENMP)")
plt.legend()

# Show the plot
plt.show()