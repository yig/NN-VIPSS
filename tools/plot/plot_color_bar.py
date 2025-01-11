import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Define the RGBColor class
class RGBColor:
    def __init__(self, r, g, b):
        self.r = r
        self.g = g
        self.b = b

    @staticmethod
    def interpolate(c1, c2, t):
        """Linear interpolation between two RGB colors."""
        return RGBColor(
            c1.r + t * (c2.r - c1.r),
            c1.g + t * (c2.g - c1.g),
            c1.b + t * (c2.b - c1.b),
        )

def blend(u):
    """
    Custom blend function based on the Mathematica code.
    
    Args:
        u (float): Input value in the range [0, 1].
        
    Returns:
        RGBColor: Interpolated color for the given u.
    """
    if u < 0.0 or u > 1.0:
        raise ValueError("u must be in the range [0, 1]")
    
    # Define the keypoints and corresponding colors
    points = [
        (0.0, RGBColor(0.0, 0.0, 9.0 / 16.0)),  # RGBColor[0, 0, 9/16]
        (1.0 / 9.0, RGBColor(0.0, 0.0, 1.0)),   # Blue
        (23.0 / 63.0, RGBColor(0.0, 1.0, 1.0)), # Cyan
        (13.0 / 21.0, RGBColor(1.0, 1.0, 0.0)), # Yellow
        (47.0 / 63.0, RGBColor(1.0, 0.647, 0.0)), # Orange
        (55.0 / 63.0, RGBColor(1.0, 0.0, 0.0)),   # Red
        (1.0, RGBColor(0.5, 0.0, 0.0)),           # RGBColor[1/2, 0, 0]
    ]

    # Find the interval for u
    for i in range(len(points) - 1):
        t1, c1 = points[i]
        t2, c2 = points[i + 1]
        if t1 <= u <= t2:
            t = (u - t1) / (t2 - t1)  # Normalize u in the interval [t1, t2]
            return RGBColor.interpolate(c1, c2, t)

    # Default (should never be reached)
    return points[-1][1]

def create_colormap(resolution=512):
    """
    Generate a colormap by evaluating the Blend function.
    
    Args:
        resolution (int): Number of colors in the colormap.
        
    Returns:
        ListedColormap: Matplotlib colormap object.
    """
    colors = [
        blend(u / (resolution - 1)) for u in range(resolution)
    ]
    rgb_values = [(c.r, c.g, c.b) for c in colors]
    return ListedColormap(rgb_values)

def plot_colorbar(colormap):
    """
    Plot a color bar using the provided colormap.
    
    Args:
        colormap (ListedColormap): The custom colormap.
    """
    fig, ax = plt.subplots(figsize=(8, 1))
    fig.subplots_adjust(bottom=0.5)

    # Create a color bar
    norm = plt.Normalize(vmin=0, vmax=0.01)
    colorbar = plt.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=colormap),
        cax=ax,
        orientation='horizontal',
    )
    colorbar.set_label('HRBR Surface Distance')
    plt.show()

if __name__ == "__main__":
    # Create and plot the colormap
    cmap = create_colormap()
    plot_colorbar(cmap)
