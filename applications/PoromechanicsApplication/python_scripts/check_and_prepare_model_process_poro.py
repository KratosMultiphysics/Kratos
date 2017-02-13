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
        #~ self.solid_model_part_names     = Parameters["problem_domain_sub_model_part_list"]
        #~ self.processes_model_part_names = Parameters["processes_sub_model_part_list"]
        
        
    def Execute(self):
        
        #construct the computing model part:
        #~ solid_parts = []
        #~ for i in range(self.solid_model_part_names.size()):
            #~ solid_parts.append(self.main_model_part.GetSubModelPart(self.solid_model_part_names[i].GetString()))
        #~ processes_parts = []
        #~ for i in range(self.processes_model_part_names.size()):
            #~ processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        #construct a model part which contains the mesh to compute
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
            
        #~ list_of_ids = set()
        #~ for part in solid_parts:
            #~ for elem in part.Elements:
                #~ list_of_ids.add(elem.Id)
        #~ poro_computing_model_part.AddElements(list(list_of_ids))
        #~ list_of_ids = set()
        #~ for part in processes_parts:
            #~ for cond in part.Conditions:
                #~ list_of_ids.add(cond.Id)
        #~ poro_computing_model_part.AddConditions(list(list_of_ids))
                
        print(poro_computing_model_part)
