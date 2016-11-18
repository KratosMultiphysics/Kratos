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
            node.SetSolutionStepValue(PRESSURE, 0, 2.0)
        
        # Linear x-reaction distribution
        self.interface_model_part.Nodes[1].SetSolutionStepValue(REACTION_X, 0, 2.0)
        self.interface_model_part.Nodes[3].SetSolutionStepValue(REACTION_X, 0, 4.7)
        self.interface_model_part.Nodes[5].SetSolutionStepValue(REACTION_X, 0, 5.6)
        self.interface_model_part.Nodes[7].SetSolutionStepValue(REACTION_X, 0, 6.6)
        self.interface_model_part.Nodes[9].SetSolutionStepValue(REACTION_X, 0, 3.6)
        
        # Linear y-reaction distribution
        self.interface_model_part.Nodes[1].SetSolutionStepValue(REACTION_Y, 0, 3.6)
        self.interface_model_part.Nodes[3].SetSolutionStepValue(REACTION_Y, 0, 6.6)
        self.interface_model_part.Nodes[5].SetSolutionStepValue(REACTION_Y, 0, 5.6)
        self.interface_model_part.Nodes[7].SetSolutionStepValue(REACTION_Y, 0, 4.7)
        self.interface_model_part.Nodes[9].SetSolutionStepValue(REACTION_Y, 0, 2.0)
        
        
    def ExecuteBeforeSolutionLoop(self):
        pass

    
    def ExecuteInitializeSolutionStep(self):
        pass


    def ExecuteFinalizeSolutionStep(self):
        
        for node in self.interface_model_part.Nodes:
            # Mapped VELOCITY check
            obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
            print(obtained_velocity_value)
            expected_velocity_value = [2*node.Y - node.Y*node.Y/2.0, 2*node.Y - node.Y*node.Y/2.0, 0.0]
            print(expected_velocity_value)
            for i in range(0,3):
                self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=0.1)
        
              
    def ExecuteBeforeOutputStep(self):
        pass


    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        pass

