import matplotlib.pyplot as plt
import numpy as np

# Create some example data
np.random.seed(0)
data = np.random.rand(10, 10)

# Create a scatter plot or imshow plot
plt.imshow(data, cmap='seismic')  # 'seismic' ranges from red to blue

# Add a colorbar with the same colormap
plt.colorbar(label='Color scale', cmap='seismic')

# Show the plot
plt.title("Colorbar from Red to Blue")
plt.show()