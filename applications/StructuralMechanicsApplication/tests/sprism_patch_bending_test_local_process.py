from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(Process, KratosUnittest.TestCase):
  
    def __init__(self,model_part,params):

        self.model_part = model_part[params["model_part_name"].GetString()]
        self.params = params
        
    def ExecuteInitialize(self):
        pass
        
    def ExecuteBeforeSolutionLoop(self):
        # Add BC
        for node in self.model_part.Nodes:
            if (node.X >2.40000e-01 -1.0e-5) | (node.Y > 1.20000e-01 -1.0e-5) | (node.X < 1.0e-5) | (node.Y < 1.0e-5):
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0, -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0, -1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
                node.Fix(DISPLACEMENT_Z)
                node.SetSolutionStepValue(DISPLACEMENT_Z, 0, 0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y ** 2))

    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        for node in self.model_part.Nodes:
            value = node.GetSolutionStepValue(DISPLACEMENT_X,0)
            self.assertAlmostEqual(value, -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
            self.assertAlmostEqual(value,-1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Z,0)
            self.assertAlmostEqual(value, 0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y **2 ))
              
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

