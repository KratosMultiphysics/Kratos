from Kratos import *
from KratosConvectionDiffusionApplication import *

def f1(Tc):
    return 10**(14.48 - 0.13858*Tc + 5.5960e-4*Tc*Tc - 7.8665e-7*Tc*Tc*Tc)

def f2(Tc):
    return 10**(53.19 - 0.2542*Tc + 2.9879e-4*Tc*Tc)
    
def AuxFunction(T):

    Tc = T - 273
    if( Tc <= 25):
        mu = 1e6
    elif (Tc > 25 and  Tc <= 200):
        mu = 1e6*(200-Tc)/175 + f1(200)
    elif(Tc > 200 and Tc <= 350):
        mu = f1(Tc)
    elif(Tc > 350 and Tc < 425 ):
        mu = f2(Tc)
    else:
        mu = f2(425)

    #mu = 100

    if(mu < 0.00):
        print "attention negative viscosity!!!!"
        mu = 0.00

    return mu
    
def CalculateViscosity(nodes):
    for node in nodes:
        rho = node.GetSolutionStepValue(DENSITY)
        T = node.GetSolutionStepValue(TEMPERATURE)
        mu = AuxFunction(T)
        node.SetSolutionStepValue(VISCOSITY,0,mu/rho)

def ApplyBoundaryconditions (nodes):
    for node in nodes:
        node.Free(TEMPERATURE);
        if( node.GetSolutionStepValue(IS_FREE_SURFACE) == True):
            node.SetSolutionStepValue(FACE_HEAT_FLUX,0,30000.0)
        else:
            node.SetSolutionStepValue(FACE_HEAT_FLUX,0,0.0)


def FixVelocities(nodes):
    for node in nodes:
        if( node.GetSolutionStepValue(IS_STRUCTURE) == 1):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)
##            node.Fix(PRESSURE)
##            node.Fix(FRACT_VEL_X)
##            node.Fix(FRACT_VEL_Y)
##            node.Fix(FRACT_VEL_Z)
        else:
            node.Free(VELOCITY_X)
            node.Free(VELOCITY_Y)
            node.Free(VELOCITY_Z)
##            node.Free(PRESSURE)
##            node.Free(FRACT_VEL_X)
##            node.Free(FRACT_VEL_Y)
##            node.Free(FRACT_VEL_Z)                
    for node in nodes:
        if( node.GetSolutionStepValue(VISCOSITY) > 100.0):
            if( node.GetSolutionStepValue(IS_FREE_SURFACE) != 1):
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)
##               node.Fix(PRESSURE)
##                node.Fix(FRACT_VEL_X)
##                node.Fix(FRACT_VEL_Y)
##                node.Fix(FRACT_VEL_Z)

            
def InitialConditions(model_part):
    model_part.Properties[1].SetValue(EMISSIVITY,1.0);
    model_part.Properties[1].SetValue(AMBIENT_TEMPERATURE,298.0);
    model_part.Properties[1].SetValue(CONVECTION_COEFFICIENT,8.0);
##    model_part.Properties[1][EMISSIVITY] = 1.0;
##    model_part.Properties[1][AMBIENT_TEMPERATURE] = 298.0;
##    model_part.Properties[1][CONVECTION_COEFFICIENT] = 8.0;
   
    gravity = Array3()
    gravity[0] = 0.00; gravity[1] = -9.81; gravity[2] = 0.0;
    for node in model_part.Nodes:
        node.SetSolutionStepValue(DENSITY,0,900.0)
        node.SetSolutionStepValue(TEMPERATURE,0,298.0)
        node.SetSolutionStepValue(BODY_FORCE,0,gravity)
        node.SetSolutionStepValue(SPECIFIC_HEAT,0,2400.0)
        node.SetSolutionStepValue(CONDUCTIVITY,0,0.25)
        node.Free(TEMPERATURE);

        
        
    

    
        
