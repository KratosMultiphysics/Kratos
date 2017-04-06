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
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        
        #comparing load and reactions
        sum_loadfy = 0
        sum_loadfx = 0
        for node in self.model_part.Nodes:
            sum_loadfy = sum_loadfy + node.GetSolutionStepValue(POINT_LOAD_Y,0)                   
        self.assertAlmostEqual(sum_loadfy,-100000)
        
        for node in self.model_part.Nodes:
            sum_loadfy = sum_loadfy + node.GetSolutionStepValue(REACTION_Y,0)
        self.assertAlmostEqual(sum_loadfy,0, places = 7)

        for node in self.model_part.Nodes:
            sum_loadfx = sum_loadfx + node.GetSolutionStepValue(REACTION_X,0)
        self.assertAlmostEqual(sum_loadfx,0, places = 7)
        
        #comparing nodal displacement to another FE-solution
        disp_y_N2_FE = -0.000266199
        disp_y_N2_model = self.model_part.Nodes[2].GetSolutionStepValue(DISPLACEMENT_Y,0)
        self.assertAlmostEqual(disp_y_N2_model,disp_y_N2_FE, places = 6)
        
        #comparing KRATOS internal displacement values
        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
            y = node.Y
            Y = node.Y0
            self.assertAlmostEqual(y - (Y + v), 0.0, places = 16)
              
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

