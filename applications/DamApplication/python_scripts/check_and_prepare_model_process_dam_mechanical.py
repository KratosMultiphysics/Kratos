import KratosMultiphysics
import KratosMultiphysics.DamApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcessDamMechanical(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class CheckAndPrepareModelProcessDamMechanical(KratosMultiphysics.Process):
    def __init__(self, main_model_part, Parameters ):
        self.main_model_part = main_model_part

        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        self.problem_domain_sub_model_part_list = Parameters["problem_domain_sub_model_part_list"]
        self.processes_sub_model_part_list = Parameters["processes_sub_model_part_list"]
        self.body_domain_sub_model_part_list = Parameters["body_domain_sub_model_part_list"]
        self.body_domain_sub_sub_model_part_list = Parameters["body_domain_sub_sub_model_part_list"]
        self.loads_sub_model_part_list = Parameters["loads_sub_model_part_list"]
        self.loads_sub_sub_model_part_list = Parameters["loads_sub_sub_model_part_list"]

    def Execute(self):
        # Construct the computing model part: a model part which contains the mesh to compute
        self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        computing_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        computing_model_part.Properties = self.main_model_part.Properties
        computing_model_part.Set(KratosMultiphysics.ACTIVE)

        domain_parts = []
        for i in range(self.problem_domain_sub_model_part_list.size()):
            domain_parts.append(self.main_model_part.GetSubModelPart(self.problem_domain_sub_model_part_list[i].GetString()))
        # Adding Nodes to Computing Model Part
        list_of_ids = set()
        for part in domain_parts:
            for node in part.Nodes:
                list_of_ids.add(node.Id)
        computing_model_part.AddNodes(list(list_of_ids))
        # Adding Elements to Computing Model Part
        list_of_ids = set()
        for part in domain_parts:
            for elem in part.Elements:
                list_of_ids.add(elem.Id)
        computing_model_part.AddElements(list(list_of_ids))
        # Adding Conditions to Computing Model Part
        domain_conditions = []
        for i in range(self.processes_sub_model_part_list.size()):
            domain_conditions.append(self.main_model_part.GetSubModelPart(self.processes_sub_model_part_list[i].GetString()))
        list_of_ids = set()
        for part in domain_conditions:
            for cond in part.Conditions:
                list_of_ids.add(cond.Id)
        computing_model_part.AddConditions(list(list_of_ids))

        # Adding Computing Sub Sub Model Parts
        # Body - Joints
        for i in range(self.body_domain_sub_model_part_list.size()):
            body_sub_model_part = self.main_model_part.GetSubModelPart(self.body_domain_sub_model_part_list[i].GetString())
            computing_model_part.CreateSubModelPart(self.body_domain_sub_sub_model_part_list[i].GetString())
            body_sub_sub_model_part = computing_model_part.GetSubModelPart(self.body_domain_sub_sub_model_part_list[i].GetString())
            list_of_ids = set()
            for node in body_sub_model_part.Nodes:
                list_of_ids.add(node.Id)
            body_sub_sub_model_part.AddNodes(list(list_of_ids))
            list_of_ids = set()
            for elem in body_sub_model_part.Elements:
                list_of_ids.add(elem.Id)
            body_sub_sub_model_part.AddElements(list(list_of_ids))
        # Arc-Length
        for i in range(self.loads_sub_model_part_list.size()):
            load_sub_model_part = self.main_model_part.GetSubModelPart(self.loads_sub_model_part_list[i].GetString())
            computing_model_part.CreateSubModelPart(self.loads_sub_sub_model_part_list[i].GetString())
            load_sub_sub_model_part = computing_model_part.GetSubModelPart(self.loads_sub_sub_model_part_list[i].GetString())
            list_of_ids = set()
            for node in load_sub_model_part.Nodes:
                list_of_ids.add(node.Id)
            load_sub_sub_model_part.AddNodes(list(list_of_ids))
