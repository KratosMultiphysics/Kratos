from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.FSIApplication as KratosFSI

KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
  
    def __init__(self,model_part,params):

        self.punctual_load_submodelpart = model_part[params["model_part_name"].GetString()]
        self.params = params
        self.updater = 1.5

        
    def ExecuteInitialize(self):
        
        # Set the initial punctual load value
        for node in self.punctual_load_submodelpart.Nodes:
            node.SetSolutionStepValue(KratosSolid.POINT_LOAD_X, 0, self.params["point_load_value"].GetDouble())
                
        
    def ExecuteBeforeSolutionLoop(self):
        pass

    
    def ExecuteInitializeSolutionStep(self):
                
        # Set the initial punctual load value
        for node in self.punctual_load_submodelpart.Nodes:
            node.SetSolutionStepValue(KratosSolid.POINT_LOAD_X, 0, self.updater*(self.params["point_load_value"].GetDouble()))


    def ExecuteFinalizeSolutionStep(self):
        pass
        
              
    def ExecuteBeforeOutputStep(self):
        pass


    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        pass

