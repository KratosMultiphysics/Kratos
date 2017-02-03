import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication 

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")     
    
    return ConstantRotationProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"

class ConstantRotationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)        
        
        
        name_of_submodelpart = settings["model_part_name"].GetString()
        self.model_part = Model[name_of_submodelpart]        
        self.cplusplus_version_process = KratosMultiphysics.PfemBaseApplication.ConstantRotationProcess(self.model_part, settings)
        
    def ExecuteInitializeSolutionStep(self):
        self.cplusplus_version_process.ExecuteInitializeSolutionStep()
        

