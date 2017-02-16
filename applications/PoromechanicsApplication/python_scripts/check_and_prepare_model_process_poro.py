import KratosMultiphysics 

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

##all the python processes should be derived from "python_process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    
    def __init__(self, main_model_part, Parameters ):
        
        self.main_model_part = main_model_part
        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        
    def Execute(self):
        
        #construct the computing model part: a model part which contains the mesh to compute
        self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        poro_computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        poro_computing_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        poro_computing_model_part.Properties  = self.main_model_part.Properties
        
        #set flag to identify the computing model part
        poro_computing_model_part.Set(KratosMultiphysics.ACTIVE)
        
        print("Adding nodes and elements to computing_model_part")
        for node in self.main_model_part.Nodes:
            poro_computing_model_part.AddNode(node,0)
        for elem in self.main_model_part.Elements:
            poro_computing_model_part.AddElement(elem,0)
        for cond in self.main_model_part.Conditions:
            poro_computing_model_part.AddCondition(cond,0)
        
        print(poro_computing_model_part)
