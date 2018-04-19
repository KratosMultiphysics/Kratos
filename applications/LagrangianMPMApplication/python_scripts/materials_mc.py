import os, sys
import math

# Importing the Kratos Library
from KratosMultiphysics import *
#from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.LagrangianMPMApplication import *
from Kratos import *
from KratosMultiphysics.StructuralApplication import *
import math



#Definicion de las propiedades de los compuestos
def AssignMaterial(Properties):
  
    Cohe      =   CohesionSoftening()#PerfectCohesion() 
    #Friction  =  FrictionSoftening()   
    #Dilatancy =  DilatancySoftening()   
    #Ft        =  ExponentialSoftening()  #LinearSoftening() #
    #Hard      =  LinearHardening()       

    #fluency_1  =  MorhCoulombYieldFunction(Cohe,Friction, Dilatancy, State.Plane_Strain, PotencialPlastic.Not_Associated)
    #fluency_2  =  IsotropicRankineYieldFunction(Ft, State.Plane_Strain)
    #fluency_3  =  ModifiedMorhCoulombYieldFunction(State.Plane_Strain, fluency_1, fluency_2) 
    #fluency_4  = EnergyYieldFunction(State.Plane_Stress);
    fluency_5  =  StandardMorhCoulombYieldFunction(Cohe, State.Plane_Strain)
    #fluency_6  =  VonMisesYieldFunction(State.Plane_Strain, Hard)

    gf   = 10.0 
    ft   = 0.35E6;
    fric = 35;
    dil  = 13;
    pi   = 3.1415926535898;
    Nphi = math.tan(pi/4 + 0.5*(pi*fric/180));
    fc   = ft * Nphi * Nphi; 
    cohe = fc/(2.00 *Nphi) 
    gc   = (fc/ft)*(fc/ft)*gf
     
    #Properties[1].SetValue(DENSITY, 100.0);
    #Properties[1].SetValue(POISSON_RATIO, 0.20);
    #Properties[1].SetValue(YOUNG_MODULUS, 28000.00);
    Properties[1].SetValue(THICKNESS, 1.00);
    #Properties[1].SetValue(DAMPING_RATIO, 0.05)

    Properties[1].SetValue(FRACTURE_ENERGY,  0.5);
    Properties[1].SetValue(CRUSHING_ENERGY,  5.0)
    Properties[1].SetValue(INTERNAL_FRICTION_ANGLE, 20.0);
    Properties[1].SetValue(DILATANCY_ANGLE, 20.0);
    Properties[1].SetValue(COHESION, 100.00);
    Properties[1].SetValue(FT, 1E15);
    Properties[1].SetValue(YIELD_STRESS, 848.7);
    Properties[1].SetValue(FC, 1E15);
    Properties[1].SetValue(ISOTROPIC_HARDENING_MODULUS, 0.00);
    
    #Properties[1].SetValue(SHEAR_STRENGTH, 1E9);    
    #Mat_1 = PlaneStrain()  
    #Mat_1 = IsotropicRankineDamage2D()
    #Mat_1 = IsotropicDamage2D(fluency_1, Ft ,Properties[1])
    #Mat_1 =  Isotropic3D() ;
    #Mat_1 = Plasticity2D(fluency_6, Properties[1])
    #Mat_1 = PlaneStrain()  
    Mat_1 = BrittleMaterial2D(fluency_5, Properties[1])
    Properties[1].SetValue(CONSTITUTIVE_LAW, Mat_1)
    
    

    
    
