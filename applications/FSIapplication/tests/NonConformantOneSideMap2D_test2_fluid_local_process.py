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
        
        # Constant pressure distribution
        for node in self.interface_model_part.Nodes:
            node.SetSolutionStepValue(PRESSURE, 0, node.Y*20)               # Linear PRESSURE variation
            node.SetSolutionStepValue(REACTION_X, 0, node.Y*30)             # Linear REACTION_X variation
            node.SetSolutionStepValue(REACTION_Y, 0, node.Y*node.Y*200)     # Parabolic REACTION_Y variation
        
        
    def ExecuteBeforeSolutionLoop(self):
        pass

    
    def ExecuteInitializeSolutionStep(self):
        pass


    def ExecuteFinalizeSolutionStep(self):
        
        for node in self.interface_model_part.Nodes:
            # Mapped VELOCITY check
            obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
            expected_velocity_value = [node.Y*10, node.Y*node.Y*10000, 0.0]
            for i in range(0,3):
                self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=0.15)
        
              
    def ExecuteBeforeOutputStep(self):
        pass


    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        pass

