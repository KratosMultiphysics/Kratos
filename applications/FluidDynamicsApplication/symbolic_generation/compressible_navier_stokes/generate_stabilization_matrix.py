from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

from sympy import *
import pprint

def computeTau(params):
    print("\nCompute Stabilization Matrix\n")
    dim = params["dim"]      # Spatial dimensions
    Tau = zeros(dim+2,dim+2) # Stabilization Matrix

    tau1 = Symbol('tau1')
    tau2 = Symbol('tau2')
    tau3 = Symbol('tau3')

    Tau[0,0] = tau1
    for i in range (0,dim):
        Tau[i+1,i+1] = tau2
    Tau[dim+1,dim+1] = tau3
    return(Tau)

def computeTauOnGaussPoint(params, U_gauss):
    print("\t- Compute Stabilization Matrix on Gauss pt.")

    # Calculate auxiliary values
    rho_g = U_gauss[0]
    e_t_g = U_gauss[params["dim"] + 1]
    norm_v_squared = 0.0
    for d in range(params["dim"]):
        norm_v_squared += (U_gauss[d + 1] * U_gauss[d + 1]) / (rho_g * rho_g)
    norm_v = sqrt(norm_v_squared)
    nu = params["mu"] / rho_g
    alpha = params["lambda"] / (rho_g * params["gamma"] * params["c_v"])

    # Calculate sound speed
    c = sqrt(params["gamma"] * (params["gamma"] -1) * ((e_t_g / rho_g) - ((1.0 / 2.0) * norm_v_squared)))

    # Calculate stabilization constants
    tau1_inv = (params["stab_c2"] * (norm_v + c)) / params["h"]
    tau2_inv = ((params["stab_c1"] / (params["h"] * params["h"])) * (4.0 * nu / 3.0)) + tau1_inv
    tau3_inv = (params["stab_c1"] * alpha / (params["h"] * params["h"])) + tau1_inv

    # Save the obtained values in the stabilization matrix
    Tau = zeros(params["dim"] + 2, params["dim"] + 2)
    Tau[0,0] = 1.0 / tau1_inv
    for i in range (params["dim"]):
        Tau[i + 1, i + 1] = 1.0 / tau2_inv
    Tau[params["dim"] + 1, params["dim"] + 1] = 1.0 / tau3_inv

    return(Tau)

def printTau(Tau, params):
    print("The Stabilization term matrix is:\n")
    dim = params["dim"]
    for i in range (0,dim+2):
        for j in range (0,dim+2):
            print("Tau[",i,",",j,"]=",Tau[i,j],"\n")

    return 0
