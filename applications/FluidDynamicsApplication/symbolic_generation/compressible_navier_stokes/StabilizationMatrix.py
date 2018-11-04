from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

def computeTau(params):
    print("\nCompute Stabilization Matrix\n")
    dim = params["dim"]				# Spatial dimensions
    Tau = zeros(dim+2,dim+2)        # Stabilization Matrix
    
    tau1 = Symbol('tau1')
    tau2 = Symbol('tau2')
    tau3 = Symbol('tau3')
 
    Tau[0,0] = tau1
    for i in range (0,dim):
        Tau[i+1,i+1] = tau2
    Tau[dim+1,dim+1] = tau3
    return(Tau)



def printTau(Tau, params):
    dim = params["dim"]				# Spatial dimensions
    print("The Stabilization term matrix is:\n")
    for i in range (0,dim+2):
        for j in range (0,dim+2):
            print("Tau[",i,",",j,"]=",Tau[i,j],"\n")

    return 0
    
    
    
    
    
    
    
    
    
    
