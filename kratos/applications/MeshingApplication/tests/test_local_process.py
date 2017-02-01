from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *

CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(Process, KratosUnittest.TestCase):
  
    def __init__(self,model_part,params):
        
        self.model_part = model_part[params["model_part_name"].GetString()]
        self.test_name = params["test_name"].GetString()
        self.params = params
       
    def ExecuteInitialize(self):
        if self.test_name == "2D_hessian_test":
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(DISTANCE, 0, math.tanh(-100.0 * (node.Y - 0.5 - 0.25 * math.sin(2.0 * node.X * math.pi))) + math.tanh(100.0*(node.Y - node.X)))
        elif self.test_name == "3D_hessian_test":
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(DISTANCE, 0, (node.Z - 0.05) * (math.tanh(-100.0 * (node.Y - 0.5 - 0.25 * math.sin(2.0 * node.X * math.pi))) + math.tanh(100.0*(node.Y - node.X))))
        
    def ExecuteBeforeSolutionLoop(self): 
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass
              
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass
