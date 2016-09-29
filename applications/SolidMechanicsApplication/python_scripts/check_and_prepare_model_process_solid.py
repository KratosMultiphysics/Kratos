import KratosMultiphysics 

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    def __init__(self, main_model_part, Parameters ):
        self.main_model_part = main_model_part
        
        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        self.solid_model_part_names     = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

        self.bodies_parts_list = []
        if( Parameters.Has("bodies_list") ):
            self.bodies_parts_list = Parameters["bodies_list"]
        
        
    def Execute(self):

        
        #construct body model parts:
        if( len(self.bodies_parts_list) != 0 ):
            for i in range(len(self.bodies_parts_list)):
                #create body model part
                body_model_part_name = self.bodies_parts_list[i]["body_name"].GetString()
                self.main_model_part.CreateSubModelPart(body_model_part_name)
                body_model_part = self.main_model_part.GetSubModelPart(body_model_part_name)
                
                print("[Model_Prepare]::Body Creation", body_model_part_name)
                body_model_part.ProcessInfo = self.main_model_part.ProcessInfo
                body_model_part.Properties  = self.main_model_part.Properties

                #build body from their parts
                body_parts_name_list = self.bodies_parts_list[i]["parts_list"]
                body_parts_list = []
                for j in range(len(body_parts_name_list)):
                    body_parts_list.append(self.main_model_part.GetSubModelPart(body_parts_name_list[j].GetString()))
                    
                for part in body_parts_list:
                    for node in part.Nodes:
                        body_model_part.Nodes.append(node)
                    for elem in part.Elements:
                        body_model_part.AddElement(elem,0)
                    for cond in part.Conditions:
                        body_model_part.AddCondition(cond,0) 

        #construct the computing model part:
        solid_parts = []
        for i in range(self.solid_model_part_names.size()):
            solid_parts.append(self.main_model_part.GetSubModelPart(self.solid_model_part_names[i].GetString()))
        
        processes_parts = []
        for i in range(self.processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        #construct a model part which contains the mesh to compute
        self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        solid_computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        solid_computing_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        solid_computing_model_part.Properties  = self.main_model_part.Properties

        #set flag to identify the solid model part
        solid_computing_model_part.Set(KratosMultiphysics.SOLID)
        #set flag to identify the computing model part
        solid_computing_model_part.Set(KratosMultiphysics.ACTIVE)
        
        for node in self.main_model_part.Nodes:
            solid_computing_model_part.AddNode(node,0)
 
        for part in solid_parts:
            part.Set(KratosMultiphysics.SOLID)
            for elem in part.Elements:
                solid_computing_model_part.AddElement(elem,0)
            #for node in part.Nodes:
            #    solid_computing_model_part.Nodes.append(node)
            
        for part in processes_parts:
            part.Set(KratosMultiphysics.BOUNDARY)
            for cond in part.Conditions:
                solid_computing_model_part.AddCondition(cond,0)  
                
        print(solid_computing_model_part)
