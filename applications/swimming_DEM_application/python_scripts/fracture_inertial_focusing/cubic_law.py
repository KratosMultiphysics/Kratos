import math
import numpy as np
from matplotlib import pyplot as plt
import root_finder
import random
import pylab

L0 = 2
H0 = 1
U0 = 0.05
epsilon = H0 / L0
epsilon_inv = 1.0 / epsilon
A1 = 0.2
A2 = 0.4
DeltaX = 0.
n_points = 1000
n_streamlines = 12
n_periods = 1
Omega = 2 * math.pi / L0 * n_periods
n_particles = 1000

delta_0 = (A1 + A2) / (2 * H0)
#alpha = 2 * math.pi * DeltaX / L0
alpha = DeltaX
gamma = (A2 - A1) / (A2 + A1)
cos_semi_alpha = math.cos(0.5 * alpha)
sin_semi_alpha = math.sin(0.5 * alpha)

   
def Phi1(X):
    return - 0.5 * H0 + A1 * math.sin(Omega * X - 0.5 * DeltaX)

def Phi2(X):
    return   0.5 * H0 + A2 * math.sin(Omega * X + 0.5 * DeltaX) 


class HorizontalDistributionMinusObjective:
    def __init__(self, L0, H0, A1, A2, DeltaX):
        self.C_phase = (A2 - A1) * math.cos(0.5 * DeltaX)
        self.C = 1.0 / (H0 * L0 + 1.0 / Omega * (A1 * math.cos(Omega * L0 - 0.5 * DeltaX) - A2 * math.cos(Omega * L0 + 0.5 * DeltaX) + self.C_phase))    
    def SetObjective(self, x):
        self.x = x
    def f(self, y):
        return self.C * (H0 * y + 1.0 / Omega * (A1 * math.cos(Omega * y - 0.5 * DeltaX) - A2 * math.cos(Omega * y + 0.5 * DeltaX) + self.C_phase)) - self.x
    def df(self, y):
        return self.C * (H0 - A1 * math.sin(Omega * y - 0.5 * DeltaX) + A2 * math.sin(Omega * y + 0.5 * DeltaX))

class phi_function:
    def __init__(self):
        pass
    def f(self, x):
        arg = Omega * L0 * x
        return delta_0 * (math.sin(arg) * cos_semi_alpha + gamma * math.cos(arg) * sin_semi_alpha)
    def df(self, x):
        arg = Omega * L0 * x
        return Omega * L0 * delta_0 * (math.cos(arg) * cos_semi_alpha - gamma * math.sin(arg) * sin_semi_alpha)
    def df2(self, x):
        arg = Omega * L0 * x       
        return (Omega * L0) ** 2 * delta_0 * (- math.sin(arg) * cos_semi_alpha - gamma * math.cos(arg) * sin_semi_alpha)
    
class h_function:
    def __init__(self):
        pass
    def f(self, x):
        arg = Omega * L0 * x
        return 0.5 + delta_0 * (math.cos(arg) * sin_semi_alpha + gamma * math.sin(arg) * cos_semi_alpha)
    def df(self, x):
        arg = Omega * L0 * x
        return Omega * L0 * delta_0 * (- math.sin(arg) * sin_semi_alpha + gamma * math.cos(arg) * cos_semi_alpha)
    def df2(self, x):
        arg = Omega * L0 * x
        return (Omega * L0) ** 2 * delta_0 * (- math.cos(arg) * sin_semi_alpha - gamma * math.sin(arg) * cos_semi_alpha)

def z_to_eta(x, z):
    phi = phi_function()
    h = h_function()
    eta = (z - phi.f(x)) / h.f(x)
    return x, eta

def eta_to_z(x, eta):
    phi = phi_function()
    h = h_function()
    z = eta * h.f(x) + phi.f(x)
    return x, z

my_h_func = h_function()
my_phi_func = phi_function()

def velocity_order_1(x, eta):    
    h = my_h_func.f(x)
    h_inv = 1.0 / h
    hdx = my_h_func.df(x)
    hdx2 = my_h_func.df2(x)
    phi = my_phi_func.f(x)
    phidx = my_phi_func.df(x)
    phidx2 = my_phi_func.df2(x)
    etadx = - (phidx + eta * hdx) * h_inv
    etadx2 = h_inv ** 2 * (hdx * (phidx + eta * hdx) - h * (phidx2 + etadx * hdx + eta * hdx2))
    one_minus_eta_2 = (1.0 - eta ** 2)
    
    UX =   0.75 * U0 * one_minus_eta_2 * h_inv
    UZ = - 0.75 * U0 * epsilon * etadx * one_minus_eta_2
    
    uxdx = - 0.75 * h_inv ** 2 * (2 * eta * etadx * h + one_minus_eta_2 * hdx)
    uxdz = - 1.5 * eta * h_inv ** 2
    uzdx = 0.75 * (2 * eta * etadx ** 2 - etadx2 * one_minus_eta_2)
    
    UXdX = U0 / L0 * uxdx
    UZdX = U0 * epsilon / L0 * uzdx
    UXdZ = U0 / H0 * uxdz
    UZdZ = - UXdX
    
    DUX = UX * UXdX + UZ * UXdZ
    DUZ = UX * UZdX + UZ * UZdZ
       
    return UX, UZ, DUX, DUZ

def GetFlowVariables(X, Z):
    x, eta = z_to_eta(X / L0, Z / H0)
    UX, UZ, DUX, DUZ = velocity_order_1(x, eta)
    D = min(abs(Phi1(X) - Z), abs(Phi2(X) - Z))
    return UX, UZ, DUX, DUZ, D
    
x_points = np.linspace(0, L0, n_points)
x_points = [x for x in x_points]
phi_1 = [Phi1(x) for x in x_points]
phi_2 = [Phi2(x) for x in x_points]
plt.axis('equal')
plt.plot(x_points, phi_1, color='k', linewidth=2)
plt.plot(x_points, phi_2, color='k', linewidth=2)

# Generate random positions
randoms_horizontal_guesses = np.random.uniform(0, 1.0 , n_particles)
objective_function = HorizontalDistributionMinusObjective(L0, H0, A1, A2, DeltaX)
i_value = 0
randoms_horizontal = [0.0 for value in randoms_horizontal_guesses]
randoms_vertical = [0.0 for value in randoms_horizontal]

for value in randoms_horizontal_guesses:
    objective_function.SetObjective(value)
    corrected_value = root_finder.FindRootBisection(objective_function, value)
    randoms_horizontal[i_value] = corrected_value
    randoms_vertical[i_value] = np.random.uniform(Phi1(corrected_value), Phi2(corrected_value))
    i_value += 1

# Streamlines
eta_values = np.linspace(-1, 1, n_streamlines)
eta_values = [value for value in eta_values]
eta_values = eta_values[1:-1]

for value in eta_values:
    streamline_Z_values = [H0 * eta_to_z(x / L0, value)[1] for x in x_points]
    plt.plot(x_points, streamline_Z_values, color='b', linestyle='dashed')
plt.scatter(randoms_horizontal, randoms_vertical)
print('delta_0 = ', delta_0)
print('alpha = ', alpha)
print('gamma = ', gamma)

for i in range(n_particles):
    X = randoms_horizontal[i]
    Z = randoms_vertical[i]
    UX, UZ, DUX, DUZ, D = GetFlowVariables(X, Z)
    pylab.arrow(X, Z, DUX, DUZ, fc = "k", ec = "k", head_width = 0.05, head_length = 0.1)
plt.show()
