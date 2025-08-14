#!/usr/bin/env python3
"""Minimal test for matplotlib 3D visualization"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Simple test data
hkl_list = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1]]
sf_values = [4.32, 4.32, 4.32, 3.05, 3.05]

print(f"HKL list: {hkl_list}")
print(f"SF values: {sf_values}")

# Convert to numpy arrays
hkl_array = np.array(hkl_list)
sf_array = np.array(sf_values)

# Extract coordinates
h_coords = hkl_array[:, 0]
k_coords = hkl_array[:, 1]
l_coords = hkl_array[:, 2]

print(f"H coords: {h_coords}")
print(f"K coords: {k_coords}")
print(f"L coords: {l_coords}")

# Create figure
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Calculate dot sizes
min_size = 20
max_size = 200
normalized_sf = sf_array / sf_array.max()
dot_sizes = min_size + (max_size - min_size) * normalized_sf

print(f"Dot sizes: {dot_sizes}")

# Create scatter plot
scatter = ax.scatter(h_coords, k_coords, l_coords, 
                    s=dot_sizes, c=sf_array, cmap='viridis', 
                    alpha=0.7, edgecolors='black', linewidths=0.5)

# Add labels
ax.set_xlabel('H (r.l.u.)')
ax.set_ylabel('K (r.l.u.)')
ax.set_zlabel('L (r.l.u.)')
ax.set_title('Structure Factors Test')

# Add colorbar
fig.colorbar(scatter, ax=ax, label='|Structure Factor|', shrink=0.6)

plt.tight_layout()
plt.show()

print("Plot should be displayed now")
