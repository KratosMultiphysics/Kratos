import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Define the parameterization
def r(s, t):
    sigma_x = 0.8
    sigma_y = 2
    mu_x = 2*np.cos(8*s)
    mu_y = 2*np.sin(8*s)
    x = mu_x + s * sigma_x * np.cos(t)
    y = mu_y + s * sigma_y * np.sin(t)
    z = 1 / (2 * np.pi * sigma_x * sigma_y) * np.exp(-s**2 / 2)
    return x, y, z, mu_x, mu_y

# Create grid for s and t
s = np.linspace(0.001, 4, 200) 
t = np.linspace(0, 2 * np.pi, 80)  # t from 0 to 2pi
S, T = np.meshgrid(s, t)

# Compute x, y, z coordinates
X, Y, Z, _, _ = r(S, T)

# Create the 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, Z, c=Z, cmap='viridis', s=2)

# Set labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Density')
ax.set_title('Bivariate Normal Distribution with Shifting Ellipse Centers')

x_range = X.max() - X.min()
y_range = Y.max() - Y.min()
max_range = max(x_range, y_range) / 2
x_center = (X.max() + X.min()) / 2
y_center = (Y.max() + Y.min()) / 2
ax.set_xlim(x_center - max_range, x_center + max_range)
ax.set_ylim(y_center - max_range, y_center + max_range)

# Show plot
plt.show()

# Plot sections of the surface to verify that they are ellipses.
# Also, plot their centers.
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot()
# for s_val in s[::15]:
#     X, Y, Z, mu_x, mu_y = r(s=s_val, t=np.linspace(0, 2*np.pi, 300))
#     ax.plot(X,Y)
#     # ax.plot(mu_x, mu_y,"x", markersize=5)
# 
# ax.set_aspect("equal")
# 
# plt.show()