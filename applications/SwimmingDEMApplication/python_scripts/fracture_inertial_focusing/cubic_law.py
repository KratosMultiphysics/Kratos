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
        self.U0 = 0.31297739987328533717
        self.epsilon = self.H0 / self.L0
        self.A1 = 0.001
        self.A2 = 0.001
        self.DeltaX = self.L0 / 2
        self.n_points = 1000
        self.n_streamlines = 12
        self.n_periods = 1
        self.Omega = 2 * math.pi / self.L0 * self.n_periods
        self.n_particles = 1000

        self.delta_0 = (self.A1 + self.A2) / (2 * self.H0)
        self.alpha = 2 * math.pi * self.DeltaX / self.L0
        #self.alpha = self.DeltaX
        self.gamma = (self.A2 - self.A1) / (self.A2 + self.A1)
        self.cos_semi_alpha = math.cos(0.5 * self.alpha)
        self.sin_semi_alpha = math.sin(0.5 * self.alpha)
        self.eta_values = []
        self.randoms_horizontal = []
        self.randoms_vertical = []
        self.vels_horizontal = []
        self.vels_vertical = []        
        self.accelerations_horizontal = []
        self.accelerations_vertical = []
        self.x_points = []
        self.phi_1 = []
        self.phi_2 = []
        self.delta_alpha = self.delta_0 * math.sin(self.alpha / 2)
        self.h3 = 8 * (1 + self.delta_alpha ** 2) / (1 - 4 * self.delta_alpha ** 2) ** 2.5
        self.J_h = 16 * math.pi ** 2 * self.delta_alpha ** 2 / (1 - 4 * self.delta_alpha ** 2) ** 1.5
        self.J_phi = 16 * math.pi ** 2 * (self.delta_0 ** 2 - self.delta_alpha ** 2) * (1 + 8 * self.delta_alpha ** 2) / (1 - 4 * self.delta_alpha ** 2) ** 2.5
        self.J_phi_h = self.gamma * self.J_h
        self.beta = self.J_phi_h / self.J_h
        
        # Plotting options
        self.plot_arrows = False


    def SetNonGeometricParameters(self, a, rho_p, rho_f, nu, gz):
        self.a = a
        self.rho = rho_p / rho_f
        self.nu = nu
        self.R = 2.0 / (2 * self.rho + 1)
        self.tau = 2.0 / 9 * self.a ** 2 / self.nu / self.L0       
        self.Fr = self.U0 ** 2 / (self.L0 * gz)
        self.Gz = 1. / self.Fr        
        self.gamma_z = 8. * self.Gz / (9 * self.epsilon * self.J_h)
        self.h_and_phi_function = h_and_phi_function(self)
        
    def PrintParameters(self):
        print('tau = ', self.tau)
        print('R = ', self.R)
        print('epsilon = ', self.epsilon)
        print('delta_0 = ', self.delta_0)
        print('h3 = ', self.h3)
        print('J_h = ', self.J_h)
        print('J_phi = ', self.J_phi)
        print('J_phi = ', self.J_phi)
        print('J_phi_h = ', self.J_phi_h)
        print('beta = ', self.beta)
        print('gamma_z = ', self.gamma_z)
        

def Phi1(pp, X):
    return - 0.5 * pp.H0 + pp.A1 * math.sin(pp.Omega * X - 0.5 * pp.alpha)

def Phi2(pp, X):
    return   0.5 * pp.H0 + pp.A2 * math.sin(pp.Omega * X + 0.5 * pp.alpha) 

class HorizontalDistributionMinusObjective:
    def __init__(self, pp, L0, H0, A1, A2, DeltaX):
        self.Omega = pp.Omega
        self.H0 = pp.H0
        self.DeltaX = DeltaX
        self.A1 = pp.A1
        self.A2 = pp.A2
        self.C_phase = (A2 - A1) * math.cos(0.5 * DeltaX)
        self.C = 1.0 / (H0 * L0 + 1.0 / pp.Omega * (A1 * math.cos(pp.Omega * L0 - 0.5 * DeltaX) - A2 * math.cos(pp.Omega * L0 + 0.5 * DeltaX) + self.C_phase))    
    def SetObjective(self, x):
        self.x = x
    def f(self, y):
        return self.C * (self.H0 * y + 1.0 / self.Omega * (self.A1 * math.cos(self.Omega * y - 0.5 * self.DeltaX) - self.A2 * math.cos(self.Omega * y + 0.5 * self.DeltaX) + self.C_phase)) - self.x
    def df(self, y):
        return self.C * (self.H0 - self.A1 * math.sin(self.Omega * y - 0.5 * self.DeltaX) + self.A2 * math.sin(self.Omega * y + 0.5 * self.DeltaX))

class h_and_phi_function:
    def __init__(self, pp):
        self.Omega = pp.Omega
        self.L0 = pp.L0
        self.sin_semi_alpha = pp.sin_semi_alpha
        self.cos_semi_alpha = pp.cos_semi_alpha
        self.delta_0 = pp.delta_0
        self.gamma = pp.gamma
        self.OmegaL0 = self.Omega * self.L0
        
    def f(self, x):
        arg = self.OmegaL0 * x
        c_arg = math.cos(arg)
        s_arg = math.sin(arg)
        phi =                        self.delta_0 * (  s_arg * self.cos_semi_alpha + self.gamma * c_arg * self.sin_semi_alpha)
        phidx =  self.OmegaL0 *      self.delta_0 * (  c_arg * self.cos_semi_alpha - self.gamma * s_arg * self.sin_semi_alpha)
        phidx2 = self.OmegaL0 ** 2 * self.delta_0 * (- s_arg * self.cos_semi_alpha - self.gamma * c_arg * self.sin_semi_alpha)
        phidx3 = self.OmegaL0 ** 3 * self.delta_0 * (- c_arg * self.cos_semi_alpha + self.gamma * s_arg * self.sin_semi_alpha)
        
        h =                    0.5 + self.delta_0 * (  c_arg * self.sin_semi_alpha + self.gamma * s_arg * self.cos_semi_alpha)
        hdx =    self.OmegaL0      * self.delta_0 * (- s_arg * self.sin_semi_alpha + self.gamma * c_arg * self.cos_semi_alpha)
        hdx2 =   self.OmegaL0 ** 2 * self.delta_0 * (- c_arg * self.sin_semi_alpha - self.gamma * s_arg * self.cos_semi_alpha)
        hdx3 =   self.OmegaL0 ** 3 * self.delta_0 * (  s_arg * self.sin_semi_alpha - self.gamma * c_arg * self.cos_semi_alpha)
        return phi, phidx, phidx2, phidx3, h, hdx, hdx2, hdx3
    
    def f_phi_h(self, x):
        arg = self.OmegaL0 * x
        c_arg = math.cos(arg)
        s_arg = math.sin(arg)
        phi =     self.delta_0 * (  s_arg * self.cos_semi_alpha + self.gamma * c_arg * self.sin_semi_alpha)
        h = 0.5 + self.delta_0 * (  c_arg * self.sin_semi_alpha + self.gamma * s_arg * self.cos_semi_alpha)
        return phi, h
    
def z_to_eta(pp, x, z):
    phi, h = pp.h_and_phi_function.f_phi_h(x)
    eta = (z - phi) / h
    return x, eta

def eta_to_z(pp, x, eta):
    phi, h = pp.h_and_phi_function.f_phi_h(x)
    z = eta * h + phi
    return x, z

def mod_eta_to_z(pp, x, eta):
    phi, h = pp.h_and_phi_function.f_phi_h(x)
    z = pp.epsilon * (eta * h + phi)
    return x, z


def velocity_order_1(pp, x, eta):    
    phi, phidx, phidx2, phidx3, h, hdx, hdx2, hdx3 = pp.h_and_phi_function.f(x)
    h_inv = 1.0 / h
    etadx = - (phidx + eta * hdx) * h_inv
    etadx2 = h_inv ** 2 * (2 * eta * hdx ** 2 + 2 * hdx * phidx - eta * h * hdx2 - h * phidx2)
    etadx3 = h_inv ** 3 * (- 6 * eta * hdx ** 3 - 6 * hdx ** 2 * phidx + 6 * eta * h * hdx * hdx2 + 3 * phidx * h * hdx2 + 3 * h * hdx * phidx2 - h ** 2 * eta * hdx3 - h ** 2 * phidx3)
    one_minus_eta_2 = (1.0 - eta ** 2)
    
    UX =   0.75 * pp.U0              * h_inv * one_minus_eta_2
    UZ = - 0.75 * pp.U0 * pp.epsilon * etadx * one_minus_eta_2
    
    uxdx = - 0.75 * h_inv ** 2 * (2 * eta * etadx * h + one_minus_eta_2 * hdx)
    uxdz =  - 1.5 * eta * h_inv ** 2
    uzdx =   0.75 * (2 * eta * etadx ** 2 - etadx2 * one_minus_eta_2)
    uxdx2 =  0.75 * h_inv ** 3 * ((2 - 12 * eta ** 2) * hdx ** 2 - 12 * eta * hdx * phidx - 2 * phidx ** 2 + h * ((3 * eta ** 2 - 1) * hdx2 + 2 * eta * phidx2))
    uzdz2 = - 1.5 * h_inv ** 2 * (2 * hdx * h_inv * eta - etadx)
    
    uxdx3 = 0.75 * h_inv ** 4 * (hdx ** 3 * (60 * eta ** 2 - 6) + 72 * eta * hdx ** 2 * phidx + 6 * hdx * (3 * phidx ** 2 + h * (hdx2 * (1 - 6 * eta ** 2) - 3 * eta * phidx2)) + h * (- 6 * phidx * phidx2 - h * hdx3 + 3 * h * eta ** 2 * hdx3 + 2 * eta * (- 9 * phidx * etadx2 + h * phidx3)))
    uxdx2z = 1.5 * h_inv ** 4 * (- 6 * hdx * phidx + 3 * eta * (- 4 * hdx ** 2 + h * hdx2) + h * phidx2)
    uxdxz2 = 4.5 * h_inv ** 4 * hdx
    
    UXdX = pp.U0 / pp.L0 * uxdx
    UZdX = pp.U0 * pp.epsilon / pp.L0 * uzdx
    UXdZ = pp.U0 / pp.H0 * uxdz
    UZdZ = - UXdX
    
    DUX = UX * UXdX + UZ * UXdZ
    DUZ = UX * UZdX + UZ * UZdZ
    
    DUX2 =   pp.U0 / pp.L0 ** 2      * uxdx2
    DUZ2 = - pp.U0 / (pp.L0 * pp.H0) * uzdz2
    
    DUX3  =   pp.U0 / pp.L0 ** 3 * uxdx3
    DUX2Z =   pp.U0 / (pp.L0 ** 2 * pp.H0) * uxdx2z
    DUZ2X = - DUX2Z
    DUZ3  = - pp.U0 / (pp.L0 * pp.H0 ** 2) * uxdxz2
    
    CONV_DERIV_LAPL_X = UX * DUX3 + UZ * DUX2Z
    CONV_DERIV_LAPL_Z = UX * DUZ2X + UZ * DUZ3
    
    FAXEN_DRAG_X = pp.a ** 2 / 6 * DUX2
    FAXEN_DRAG_Z = pp.a ** 2 / 6 * DUZ2
    
    FAXEN_ADDED_MASS_X = pp.a ** 2 / 30 * CONV_DERIV_LAPL_X
    FAXEN_ADDED_MASS_Z = pp.a ** 2 / 30 * CONV_DERIV_LAPL_Z
    
    UX += FAXEN_DRAG_X
    UZ += FAXEN_DRAG_Z
    
    DUX += FAXEN_DRAG_X
    DUZ += FAXEN_DRAG_Z    
    
    return UX, UZ, DUX, DUZ, DUX2, DUZ2

def GetFlowVariables(pp, i, X, Z):
    x, eta = z_to_eta(pp, X / pp.L0, Z / pp.H0)
    UX, UZ, DUX, DUZ, DUX2, DUZ2 = velocity_order_1(pp, x, eta)
    D = min(abs(Phi1(pp, X) - Z), abs(Phi2(pp, X) - Z))
    pp.randoms_horizontal[i]       = X
    pp.randoms_vertical[i]         = Z
    pp.vels_horizontal[i]          = UX
    pp.vels_vertical[i]            = UZ      
    pp.accelerations_horizontal[i] = DUX
    pp.accelerations_vertical[i]   = DUZ   
    return UX, UZ, DUX, DUZ, DUX2, DUZ2, D    

def GenerateRandomPositions(pp, n_particles):
    pp.x_points = np.linspace(0, pp.L0, pp.n_points)
    pp.x_points = [x for x in pp.x_points]
    pp.phi_1 = [Phi1(pp, x) for x in pp.x_points]
    pp.phi_2 = [Phi2(pp, x) for x in pp.x_points]
    pp.n_particles = n_particles
    randoms_horizontal_guesses = np.random.uniform(0, 1.0 , n_particles)
    objective_function = HorizontalDistributionMinusObjective(pp, pp.L0, pp.H0, pp.A1, pp.A2, pp.DeltaX)
    pp.randoms_horizontal       = [0.0 for i in range(n_particles)]
    pp.randoms_vertical         = [0.0 for i in range(n_particles)]
    pp.vels_horizontal          = [0.0 for i in range(n_particles)]
    pp.vels_vertical            = [0.0 for i in range(n_particles)]    
    pp.accelerations_horizontal = [0.0 for i in range(n_particles)]
    pp.accelerations_vertical   = [0.0 for i in range(n_particles)]
    
    i_value = 0

    for value in randoms_horizontal_guesses:
        objective_function.SetObjective(value)
        corrected_value = root_finder.FindRootBisection(objective_function, value)
        pp.randoms_horizontal[i_value] = corrected_value
        pp.randoms_vertical[i_value] = np.random.uniform(Phi1(pp, corrected_value), Phi2(pp, corrected_value))
        i_value += 1

    # Streamlines
    pp.eta_values = np.linspace(-1, 1, pp.n_streamlines)
    pp.eta_values = [value for value in pp.eta_values]
    pp.eta_values = pp.eta_values[1:-1]

def GetPositionAndFlowVariables(pp, i):
    X = pp.randoms_horizontal[i]
    Z = pp.randoms_vertical[i] 
    UX, UZ, DUX, DUZ, DUX2, DUZ2, D = GetFlowVariables(pp, i, X, Z)
    return X, Z, UX, UZ, DUX, DUZ, D

def PrintResult(pp, time):

    for value in pp.eta_values:
        streamline_Z_values = [pp.H0 * eta_to_z(pp, x / pp.L0, value)[1] for x in pp.x_points]
        plt.plot(pp.x_points, streamline_Z_values, color = 'b', linestyle = 'dashed')
    
    eta_critical = 0.245806
    streamline_Z_values = [-pp.L0 * mod_eta_to_z(pp, x / pp.L0, eta_critical)[1] for x in pp.x_points]
    plt.plot(pp.x_points, streamline_Z_values, color = 'r')
    plt.scatter(pp.randoms_horizontal, pp.randoms_vertical)


    for i in range(pp.n_particles):
        X = pp.randoms_horizontal[i]
        Z = pp.randoms_vertical[i]      
        UX, UZ, DUX, DUZ, DUX2, DUZ2, D = GetFlowVariables(pp, i, X, Z)
        pp.vels_horizontal[i]          = UX
        pp.vels_vertical[i]            = UZ        
        pp.accelerations_horizontal[i] = DUX
        pp.accelerations_vertical[i]   = DUZ
    if pp.plot_arrows:
        acc_moduli_inv = [1.0 / math.sqrt(pp.accelerations_horizontal[i] ** 2 + pp.accelerations_vertical[i] ** 2) for i in range(pp.n_particles)]
        acc_max_modul_inv = min(acc_moduli_inv)
        acc_size_coeff = acc_max_modul_inv * 0.25 * pp.H0
        
        for i in range(pp.n_particles): 
            X = pp.randoms_horizontal[i]
            Z = pp.randoms_vertical[i]  
            DUX = pp.accelerations_horizontal[i]
            DUZ = pp.accelerations_vertical[i]    
            UX = pp.vels_horizontal[i]
            UZ = pp.vels_vertical[i]    
            pylab.arrow(X, Z, UX * acc_size_coeff * 10, UZ * acc_size_coeff * 10, fc = "k", ec = "k", width = 0.001 * pp.H0, head_width = 0.02 * pp.H0, head_length = 0.04 * pp.H0)
    plt.axis('equal')
    plt.plot(pp.x_points, pp.phi_1, color='k', linewidth=2)
    plt.plot(pp.x_points, pp.phi_2, color='k', linewidth=2)
    plt.savefig('fracture_' + str(time) + '.png', bbox_inches='tight')
    plt.close()
