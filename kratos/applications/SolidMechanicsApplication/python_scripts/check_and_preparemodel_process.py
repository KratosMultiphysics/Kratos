import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication 

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    def __init__(self, main_model_part, Parameters ):
        self.main_model_part = main_model_part
        
        self.solid_model_part_names = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

    def Execute(self):
        
        solid_parts = []
        for i in range(self.solid_model_part_names.size()):
            solid_parts.append(self.main_model_part.GetSubModelPart(self.solid_model_part_names[i].GetString()))
        
        processes_parts = []
        for i in range(self.processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        #construct a model part which contains both the skin and the volume
        #temporarily we call it "solid_computational_model_part"
        self.main_model_part.CreateSubModelPart("solid_computational_model_part")
        solid_computational_model_part = self.main_model_part.GetSubModelPart("solid_computational_model_part")
        solid_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        
        for part in solid_parts:
            for node in part.Nodes:
                solid_computational_model_part.Nodes.append(node)
            for elem in part.Elements:
                solid_computational_model_part.Elements.append(elem)
            
        for part in processes_parts:
            for cond in part.Conditions:
                solid_computational_model_part.Conditions.append(cond)  
                
        print(solid_computational_model_part)
        
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            #verify that the skin is correct (no gaps and overlaps)
            KratosMultiphysics.CheckSkinProcess(solid_computational_model_part, KratosMultiphysics.Flags()).Execute()

        #~ #verify the orientation of the skin
        #~ throw_errors = False
        #~ KratosMultiphysics.TetrahedralMeshOrientationCheck(solid_computational_model_part,throw_errors).Execute()
