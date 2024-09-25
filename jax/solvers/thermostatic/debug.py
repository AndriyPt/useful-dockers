import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation

def tensor_map_func(x, y):
    radius = 0.4
    result = np.where((np.power(x - 1.0, 2) + np.power(y - 1.0, 2)) < radius, 
                        (np.cos(0.5 * np.pi / radius * (np.power(x - 1.0, 2) + np.power(y - 1.0, 2))) + 1.0),
                        1.0)
    return result

# def tensor_map_func(x, y):
#     radius = 0.4
#     return  (np.cos(0.5 * np.pi / radius * (np.power(x - 1.0, 2) + np.power(y - 1.0, 2))) + 1.0)

def f(x, y):
	return np.sin(np.sqrt(x ** 2 + y ** 2))

x = np.linspace(0, 2, 30)
y = np.linspace(0, 2, 30)
X, Y = np.meshgrid(x, y)
Z = tensor_map_func(X, Y)

tri = Triangulation(X.ravel(), Y.ravel())

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

ax.plot_trisurf(tri, Z.ravel(), cmap='cool', edgecolor='none', alpha=0.8)

ax.set_title('Surface Triangulation Plot of Lambda', fontsize=14)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.set_zlabel('z', fontsize=12)

plt.show()
