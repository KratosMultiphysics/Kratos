import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as GammaFunction
from KratosMultiphysics.LaserDrillingApplication.laserdrilling_transient_solver_ablation_plus_thermal import (
    LaserDrillingTransientSolverAblationPlusThermal,
)


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
    omega_0 = parameters["omega_0"]
    n = parameters["gaussian_order"]

    fluence = np.exp(-2.0 * (r / omega_0) ** n)  # *F_p

    return fluence


n = 100
x_arr = np.linspace(0, 0.03, n)

v = 4

parameters = {
    "gamma": 0.005,
    "omega_0": 0.02,
    "gaussian_order": 2,
    "average_laser_power": 5,
    "pulse_frequency": 2e5,
}

Q = parameters["average_laser_power"] / parameters["pulse_frequency"]

F_p = ComputePeakFluenceSuperGaussian(parameters)
parameters["F_p"] = F_p
print(f"{F_p=}")


# Ablation threshold
delta_pen = 10.0  # J/cm2 ?
q_ast = 5e-4  # m ?
F_thresh = delta_pen * q_ast

LDTSAPT = LaserDrillingTransientSolverAblationPlusThermal()
CauchyPulseNormalized = LDTSAPT.NormalizeAxisymmetricFunction(
    lambda r: CauchyPulse(r, parameters), rmin=0, rmax=x_arr.max()
)

y_arr = CauchyPulse(x_arr, parameters)
# y_arr = FluenceSuperGaussianNotNormalized(x_arr, parameters)
y_arr_normalized = CauchyPulseNormalized(x_arr)

y_pulse = y_arr_normalized * Q

print("x_arr")
for x in x_arr:
    print(f"{x},")

print("Cauchy")
for y in y_arr:
    print(f"{y},")

print("Cauchy normalized")
for y in y_arr_normalized:
    print(f"{y},")

plt.plot(x_arr, y_arr, ".", label="Raw curve")
plt.plot(x_arr, y_pulse, ".", label=f"Cauchy pulse {Q=}")
plt.hlines(y=F_thresh, xmin=x_arr.min(), xmax=x_arr.max(), label="delta_pen*q_ast")
plt.legend()

plt.show()

# Combine x and y into a 2D array
data = np.column_stack((x_arr, y_arr))

# Save to CSV
np.savetxt("fluence_table.csv", data, delimiter=",", header="", fmt="%s", comments="")
