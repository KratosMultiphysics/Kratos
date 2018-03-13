from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
      
def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SPRISMProcess(Model, settings["Parameters"])

class SPRISMProcess(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):

        self.model_part = model_part[params["model_part_name"].GetString()]
        self.params = params
        
    def ExecuteInitialize(self):
        sprism_neighbour_search = KratosMultiphysics.StructuralMechanicsApplication.SprismNeighbours(self.model_part)
        sprism_neighbour_search.Execute()
    
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