import scipy as spy 
import pywavefront
import numpy as np
from scipy.spatial import ConvexHull

obj_path = r'out/split_pts.obj'
scene = pywavefront.Wavefront(obj_path)
print(scene.vertices) 

pts = np.array(scene.vertices)
print(pts.shape)

hull = ConvexHull(pts)
print(hull.volume)





