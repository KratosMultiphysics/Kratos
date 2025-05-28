import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as GammaFunction


def ComputePeakFluenceSuperGaussian(parameters):
    omega_0 = parameters["omega_0"]
    b = parameters["gaussian_order"]
    average_laser_power = parameters["average_laser_power"]
    pulse_frequency = parameters["pulse_frequency"]

    Q = average_laser_power / pulse_frequency
    print(f"{Q=}")
    F_p = Q * b * 2 ** (2 / b) / (2 * np.pi * omega_0**2 * GammaFunction(2 / b))
    return F_p  # J/mm2


def CauchyPulse(r, parameters):
    gamma = parameters["gamma"]
    F_p = parameters["F_p"]
    denominator = 1 + np.power(r / gamma, 2)
    return F_p / denominator


def FluenceSuperGaussianNotNormalized(r, parameters):
    F_p = parameters["F_p"]
    omega_0 = parameters["omega_0"]
    n = parameters["gaussian_order"]

    fluence = np.exp(-2.0 * (r / omega_0) ** n)  # *F_p

    return fluence


n = 100
x_arr = np.linspace(0, 0.03, n)

v = 4

parameters = {
    "gamma": 0.005,
    "omega_0": 0.018,
    "gaussian_order": 2,
    "average_laser_power": 5,
    "pulse_frequency": 2e5,
}

F_p = ComputePeakFluenceSuperGaussian(parameters)
parameters["F_p"] = F_p
print(f"{F_p=}")

# y_arr = CauchyPulse(x_arr, parameters)
y_arr = FluenceSuperGaussianNotNormalized(x_arr, parameters)

print(x_arr)
print(y_arr)
plt.plot(x_arr, y_arr)
plt.show()

# Combine x and y into a 2D array
data = np.column_stack((x_arr, y_arr))

# Save to CSV
np.savetxt("fluence_table.csv", data, delimiter=",", header="", fmt="%s", comments="")
