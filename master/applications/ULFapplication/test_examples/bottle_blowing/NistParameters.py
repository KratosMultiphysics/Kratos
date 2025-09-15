from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
#from KratosConvectionDiffusionApplication import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.MeshingApplication import *
#from KratosMultiphysics.PFEMApplication import PfemUtils
#from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

import math
#def AuxFunction(T):
#    nu = 1.0*(10**13)*(2.718282**(-0.046*T))
#    mu=nu 		
#    if(mu < 0.0):
#        print "attention negative viscosity!!!!"
#        mu = 0.0
#    return mu

#def CalculateViscosity(model_part):
#    for node in model_part.Nodes:
#        rho = node.GetSolutionStepValue(DENSITY)
#        T = node.GetSolutionStepValue(TEMPERATURE)
#        mu = AuxFunction(T)
#        node.SetSolutionStepValue(VISCOSITY,0,mu)



def AuxFunction(T):
    #print "Computing viscosity"
    
    #if (T>800):
      #nu=25834736472.9798*math.exp(-0.0140237333*T)/2400
      #nu=441203477*math.exp(-0.01965*T)
      #265677693762693*math.exp(-0.0233569026*T)/2400
      #divided by 2400 to get kinematic viscosiyu. 2400=density
     # nu=math.exp(1.8+4700/(T-220))/2500
    #elif (T<=800):
      #nu=57083
    #  nu=100

    #just to try
    #nu=0.3
    #if (T>250):
    #    nu=-0.05145*T+29.2913
    #if (T>494):
    #    nu=-0.01977*T+13.64138
    #if (T>690.004):
    #    nu=0

    #mu=10**nu


    #nu = 6*(10**13)*(2.718282**(-0.044*T))
    ##nu = 1.0*(10**13)*(2.718282**(-0.046*T)) no sirve
    #mu=nu 		
    #if(mu < 0.0):
        #print "attention negative viscosity!!!!"
        #mu = 0.0
    #print "Finished computing viscosity"
    nu=math.exp(1.8+4700/(T-220))/2500
    return nu

def CalculateViscosity(model_part):
    for node in model_part.Nodes:
        #rho = node.GetSolutionStepValue(DENSITY)
        T = node.GetSolutionStepValue(TEMPERATURE)
        nu = AuxFunction(T)
        node.SetSolutionStepValue(VISCOSITY,0,nu)






def FixVelocities(nodes):
    for node in nodes:
        if( node.GetSolutionStepValue(IS_STRUCTURE) == 1):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)
        else:
            node.Free(VELOCITY_X)
            node.Free(VELOCITY_Y)
            node.Free(VELOCITY_Z)
    for node in nodes:
        if( node.GetSolutionStepValue(VISCOSITY) > 100.0):
            if( node.GetSolutionStepValue(IS_FREE_SURFACE) != 1):
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)
            
def InitialConditions(model_part):
    model_part.Properties[1].SetValue(EMISSIVITY,1.0);
    model_part.Properties[1].SetValue(AMBIENT_TEMPERATURE,298.0);
    model_part.Properties[1].SetValue(CONVECTION_COEFFICIENT,0.0);
    gravity = Array3()
    gravity[0] = 0.00; gravity[1] = -9.81; gravity[2] = 0.0;
    for node in model_part.Nodes:
        node.SetSolutionStepValue(DENSITY,0,1130.0)
        node.SetSolutionStepValue(TEMPERATURE,0,298.0)
        node.SetSolutionStepValue(BODY_FORCE,0,gravity)
        node.SetSolutionStepValue(SPECIFIC_HEAT,0,1740.0)
        node.SetSolutionStepValue(CONDUCTIVITY,0,0.25)
        node.Free(TEMPERATURE);

        
    

    
        
