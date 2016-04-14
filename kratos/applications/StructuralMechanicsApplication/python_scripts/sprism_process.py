from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

def Factory(settings, Model):
	
    if settings["process_name"] == "SPRISM_process":
        my_process_parameters = settings["parameters"]
        model_part = Model.get(my_process_parameters.get("model_part_name","MODEL PART WAS NOT FOUND")," TEST ") 
	    
        return SPRISM_process(model_part)

class SPRISM_process:
    def __init__(self,model_part):
        self.model_part =  model_part
        
    def ExecuteInitialize(self):
        # Find neighbours
        sprism_neighbour_search = SprismNeighbours(self.model_part)
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
    
    def Clear(self):
        pass
 
