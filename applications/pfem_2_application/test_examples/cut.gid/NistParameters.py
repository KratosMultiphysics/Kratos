from Kratos import *
from KratosConvectionDiffusionApplication import *

from Kratos import *
from KratosConvectionDiffusionApplication import *

def AuxFunction2(T):
    
    nu = 1.0*(10**13)*(2.718282**(-0.046*T))
    mu=nu 		
    if(mu < 0.0):
        print ("attention negative viscosity!!!!")
        mu = 0.0
    return mu

def CalculateViscosity(model_part):
    for node in model_part.Nodes:
        rho = node.GetSolutionStepValue(DENSITY)
        T = node.GetSolutionStepValue(YCH4)
        mu = AuxFunction(T)
        node.SetSolutionStepValue(VISCOSITY,0,mu)
import math
def AuxFunction(T):
    Tc = T - 273.0
    e=0.0
    mu=0.0
    if( Tc <= 25):
        mu = 1e6
    if (Tc > 25.0 and Tc <= 200.0):
        e=14.48 - 0.13858*200.0 + 5.5960e-4*200.0*200.0 - 7.8665e-7 * 200.0 * 200.0 * 200.0
        mu = 1e6*(200.0-Tc)/175 + pow(10,e)
    if(Tc > 200.0 and Tc <= 350.0):
        e=14.48 - 0.13858*Tc + 5.5960e-4*Tc*Tc - 7.8665e-7*Tc*Tc*Tc
        mu = pow(10,e)
    if(Tc > 350.0 and Tc < 425.0):
        e=53.19 - 0.2542*Tc + 2.9879e-4*Tc*Tc
        mu = pow(10,e)
    if (Tc >= 425.0):
        e=53.19 - 0.2542*425.0 + 2.9879e-4*425.0*425.0
        mu = pow(10,e)
    return mu


