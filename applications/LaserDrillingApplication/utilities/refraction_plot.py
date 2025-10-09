import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Set up the plotting style
plt.figure(figsize=(12, 8))

# Create x-axis data points
x = np.linspace(-4, 4, 1000)

# Generate Gaussian curve
mean = 0
std = 1
gaussian = norm.pdf(x, mean, std)
gaussian = -gaussian * 3 + 2  # Adjust curve position and amplitude

plt.rcParams.update({
    "text.usetex": False,
    "savefig.format": "png",
    "font.size": 20.0,
    "axes.labelsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20,
    "figure.titlesize": 20,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]
})
# Plot material interface (Gaussian curve)
plt.plot(x, gaussian, 'k-', linewidth=2, label='Material interface')

# Fill regions above and below the curve to represent different materials
plt.fill_between(x, gaussian, 4, color='lightyellow', alpha=0.3, label='Air (n=1)')
plt.fill_between(x, gaussian, -1, color='lightblue', alpha=0.3, label='Resin (n=1.5)')

# Define starting points for incident light
light_sources = np.array([-2, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2])  # x-coordinates of light sources
y_start = 3.5  # y-coordinate of light sources

# Calculate and plot incident and refracted rays
for x_source in light_sources:
    # Calculate the intersection point on the Gaussian curve
    y_intersect = norm.pdf(x_source, mean, std) * -3 + 2

    # Plot incident light ray (vertical line)
    plt.plot([x_source, x_source], [y_start, y_intersect],
             'y-', linewidth=2, alpha=0.8)

    # Calculate the slope of the Gaussian curve at the intersection point (derivative)
    dx = 0.0001
    dy = (norm.pdf(x_source + dx, mean, std) - norm.pdf(x_source, mean, std)) * -3
    gaussian_slope = dy / dx

    # Calculate the angle of the normal to the Gaussian curve (with respect to the vertical direction)
    normal_angle = np.arctan(gaussian_slope)
    normal_angle_in_deg = normal_angle * 180 / np.pi
    print(normal_angle_in_deg)

    # Calculate the incident angle (angle between the incident ray and the normal)
    incident_angle = np.abs(normal_angle)  # Light ray is vertical

    # Define refractive indices for air and resin
    n1, n2 = 1.0, 1.5  # Refractive indices for air and resin

    # Calculate the refracted angle using Snell's law
    refraction_angle = np.arcsin(n1 / n2 * np.sin(incident_angle))

    # Calculate endpoint for refracted ray
    length = 2  # Length of the refracted ray
    x_end = x_source + length * np.sin(refraction_angle) * np.sign(gaussian_slope)
    y_end = y_intersect - length * np.cos(refraction_angle)

    # Plot refracted ray
    plt.plot([x_source, x_end], [y_intersect, y_end],
             'r-', linewidth=2, alpha=0.8)

# Set up chart properties
plt.grid(True, alpha=0.3)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Refraction of light through a symmetric Gaussian interface')
plt.legend()
plt.axis('equal')
plt.ylim(-1, 4)
plt.tight_layout()
# Display the plot
plt.show()