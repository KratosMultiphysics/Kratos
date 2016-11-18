from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.FSIApplication import *

CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(Process, KratosUnittest.TestCase):
  
    def __init__(self,model_part,params):

        self.interface_model_part = model_part[params["model_part_name"].GetString()]
        self.params = params

        
    def ExecuteInitialize(self):
        # Parabolic velocity distribution
        for node in self.interface_model_part.Nodes:
            node.SetSolutionStepValue(VELOCITY_X, 0, node.Y*10)             # Linear VELOCITY_X distribution
            node.SetSolutionStepValue(VELOCITY_Y, 0, node.Y*node.Y*10000)   # Parabolic VELOCITY_Y distribution
            node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)

        
    def ExecuteBeforeSolutionLoop(self):
        pass

    
    def ExecuteInitializeSolutionStep(self):
        pass


    def ExecuteFinalizeSolutionStep(self):
        # Mapped PRESSURE check
        for node in self.interface_model_part.Nodes:
            obtained_pressure_value = node.GetSolutionStepValue(PRESSURE,0)
            expected_pressure_value = node.Y*20
            self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=0.1)
            
        # Force equilibrium check (mapped POINT_LOAD)        
        expected_sum_fx = 765.0
        expected_sum_fy = 858.375
        expected_sum_fz = 0.0
        
        obtained_sum_fx = 0.0
        obtained_sum_fy = 0.0
        obtained_sum_fz = 0.0

        for node in self.interface_model_part.Nodes:
            f_solid = node.GetSolutionStepValue(POINT_LOAD)
            obtained_sum_fx += f_solid[0]
            obtained_sum_fy += f_solid[1]
            obtained_sum_fz += f_solid[2]
            
        self.assertAlmostEqual(obtained_sum_fx, expected_sum_fx, delta=0.1)
        self.assertAlmostEqual(obtained_sum_fy, expected_sum_fy, delta=0.1)
        self.assertAlmostEqual(obtained_sum_fz, expected_sum_fz, delta=0.1)
        
              
    def ExecuteBeforeOutputStep(self):
        pass


    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        pass

