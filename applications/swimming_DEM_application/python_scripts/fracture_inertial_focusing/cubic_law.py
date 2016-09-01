import math
import numpy as np
from matplotlib import pyplot as plt
import root_finder
import random
import pylab

class ProblemParameters:
    def __init__(self):
        self.L0 = 0.01
        self.H0 = 0.005
        self.U0 = 0.01
        self.epsilon = self.H0 / self.L0
        self.A1 = 0.001
        self.A2 = 0.0005
        self.DeltaX = 0.2
        self.n_points = 1000
        self.n_streamlines = 12
        self.n_periods = 4
        self.Omega = 2 * math.pi / self.L0 * self.n_periods
        self.n_particles = 1000

        self.delta_0 = (self.A1 + self.A2) / (2 * self.H0)
        #alpha = 2 * math.pi * self.DeltaX / self.L0
        self.alpha = self.DeltaX
        self.gamma = (self.A2 - self.A1) / (self.A2 + self.A1)
        self.cos_semi_alpha = math.cos(0.5 * self.alpha)
        self.sin_semi_alpha = math.sin(0.5 * self.alpha)

pp = ProblemParameters()

def GetProblemParameters():
    return pp

def Phi1(X):
    return - 0.5 * pp.H0 + pp.A1 * math.sin(pp.Omega * X - 0.5 * pp.DeltaX)

def Phi2(X):
    return   0.5 * pp.H0 + pp.A2 * math.sin(pp.Omega * X + 0.5 * pp.DeltaX) 


class HorizontalDistributionMinusObjective:
    def __init__(self, L0, H0, A1, A2, DeltaX):
        self.C_phase = (A2 - A1) * math.cos(0.5 * DeltaX)
        self.C = 1.0 / (H0 * L0 + 1.0 / pp.Omega * (A1 * math.cos(pp.Omega * L0 - 0.5 * DeltaX) - A2 * math.cos(pp.Omega * L0 + 0.5 * DeltaX) + self.C_phase))    
    def SetObjective(self, x):
        self.x = x
    def f(self, y):
        return self.C * (pp.H0 * y + 1.0 / pp.Omega * (pp.A1 * math.cos(pp.Omega * y - 0.5 * pp.DeltaX) - pp.A2 * math.cos(pp.Omega * y + 0.5 * pp.DeltaX) + self.C_phase)) - self.x
    def df(self, y):
        return self.C * (pp.H0 - pp.A1 * math.sin(pp.Omega * y - 0.5 * pp.DeltaX) + pp.A2 * math.sin(pp.Omega * y + 0.5 * pp.DeltaX))

class phi_function:
    def __init__(self):
        pass
    def f(self, x):
        arg = pp.Omega * pp.L0 * x
        return pp.delta_0 * (math.sin(arg) * pp.cos_semi_alpha + pp.gamma * math.cos(arg) * pp.sin_semi_alpha)
    def df(self, x):
        arg = pp.Omega * pp.L0 * x
        return pp.Omega * pp.L0 * pp.delta_0 * (math.cos(arg) * pp.cos_semi_alpha - pp.gamma * math.sin(arg) * pp.sin_semi_alpha)
    def df2(self, x):
        arg = pp.Omega * pp.L0 * x       
        return (pp.Omega * pp.L0) ** 2 * pp.delta_0 * (- math.sin(arg) * pp.cos_semi_alpha - pp.gamma * math.cos(arg) * pp.sin_semi_alpha)
    
class h_function:
    def __init__(self):
        pass
    def f(self, x):
        arg = pp.Omega * pp.L0 * x
        return 0.5 + pp.delta_0 * (math.cos(arg) * pp.sin_semi_alpha + pp.gamma * math.sin(arg) * pp.cos_semi_alpha)
    def df(self, x):
        arg = pp.Omega * pp.L0 * x
        return pp.Omega * pp.L0 * pp.delta_0 * (- math.sin(arg) * pp.sin_semi_alpha + pp.gamma * math.cos(arg) * pp.cos_semi_alpha)
    def df2(self, x):
        arg = pp.Omega * pp.L0 * x
        return (pp.Omega * pp.L0) ** 2 * pp.delta_0 * (- math.cos(arg) * pp.sin_semi_alpha - pp.gamma * math.sin(arg) * pp.cos_semi_alpha)

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
    
    UX =   0.75 * pp.U0 * one_minus_eta_2 * h_inv
    UZ = - 0.75 * pp.U0 * pp.epsilon * etadx * one_minus_eta_2
    
    uxdx = - 0.75 * h_inv ** 2 * (2 * eta * etadx * h + one_minus_eta_2 * hdx)
    uxdz = - 1.5 * eta * h_inv ** 2
    uzdx = 0.75 * (2 * eta * etadx ** 2 - etadx2 * one_minus_eta_2)
    
    UXdX = pp.U0 / pp.L0 * uxdx
    UZdX = pp.U0 * pp.epsilon / pp.L0 * uzdx
    UXdZ = pp.U0 / pp.H0 * uxdz
    UZdZ = - UXdX
    
    DUX = UX * UXdX + UZ * UXdZ
    DUZ = UX * UZdX + UZ * UZdZ
       
    return UX, UZ, DUX, DUZ

def GetFlowVariables(X, Z):
    x, eta = z_to_eta(X / pp.L0, Z / pp.H0)
    UX, UZ, DUX, DUZ = velocity_order_1(x, eta)
    D = min(abs(Phi1(X) - Z), abs(Phi2(X) - Z))
    return UX, UZ, DUX, DUZ, D
    
x_points = np.linspace(0, pp.L0, pp.n_points)
x_points = [x for x in x_points]
phi_1 = [Phi1(x) for x in x_points]
phi_2 = [Phi2(x) for x in x_points]
plt.axis('equal')
plt.plot(x_points, phi_1, color='k', linewidth=2)
plt.plot(x_points, phi_2, color='k', linewidth=2)

# Generate random positions
randoms_horizontal_guesses = np.random.uniform(0, 1.0 , pp.n_particles)
objective_function = HorizontalDistributionMinusObjective(pp.L0, pp.H0, pp.A1, pp.A2, pp.DeltaX)
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
eta_values = np.linspace(-1, 1, pp.n_streamlines)
eta_values = [value for value in eta_values]
eta_values = eta_values[1:-1]

def GetPositionAndFlowVariables(i):
    X = randoms_horizontal[i]
    Z = randoms_vertical[i] 
    UX, UZ, DUX, DUZ, D = GetFlowVariables(X, Z)
    return X, Z, UX, UZ, DUX, DUZ, D

for value in eta_values:
    streamline_Z_values = [pp.H0 * eta_to_z(x / pp.L0, value)[1] for x in x_points]
    plt.plot(x_points, streamline_Z_values, color='b', linestyle='dashed')
plt.scatter(randoms_horizontal, randoms_vertical)
print('delta_0 = ', pp.delta_0)
print('alpha = ', pp.alpha)
print('gamma = ', pp.gamma)
derivatives_horizontal = [0.0 for i in range(pp.n_particles)]
derivatives_vertical = [0.0 for i in range(pp.n_particles)]
for i in range(pp.n_particles):
    X = randoms_horizontal[i]
    Z = randoms_vertical[i]      
    UX, UZ, DUX, DUZ, D = GetFlowVariables(X, Z)
    derivatives_horizontal[i] = DUX
    derivatives_vertical[i] = DUZ

moduli_inv = [1.0 / math.sqrt(derivatives_horizontal[i] ** 2 + derivatives_vertical[i] ** 2) for i in range(pp.n_particles)]
max_modul_inv = min(moduli_inv)
print(moduli_inv)
print('max_modul_inv',max_modul_inv)
size_coeff = max_modul_inv * 0.25 * pp.H0
def PrintResult():
    for i in range(pp.n_particles): 
        X = randoms_horizontal[i]
        Z = randoms_vertical[i]  
        DUX = derivatives_horizontal[i]
        DUZ = derivatives_vertical[i]    
        pylab.arrow(X, Z, DUX * size_coeff, DUZ * size_coeff, fc = "k", ec = "k", width = 0.001 * pp.H0, head_width = 0.02 * pp.H0, head_length = 0.04 * pp.H0)
    plt.show()
