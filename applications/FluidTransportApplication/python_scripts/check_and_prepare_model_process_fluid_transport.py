import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):

    def __init__(self, main_model_part, Parameters ):
        KratosMultiphysics.Process.__init__(self)
        self.main_model_part = main_model_part

        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        self.problem_domain_sub_model_part_list = Parameters["problem_domain_sub_model_part_list"]
        self.processes_sub_model_part_list = Parameters["processes_sub_model_part_list"]

    def Execute(self):

        #construct the computing model part: a model part which contains the mesh to compute
        self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        fluid_transport_computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        fluid_transport_computing_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        fluid_transport_computing_model_part.Properties  = self.main_model_part.Properties

        #set flag to identify the computing model part
        fluid_transport_computing_model_part.Set(KratosMultiphysics.ACTIVE)

        domain_parts = []
        for i in range(self.problem_domain_sub_model_part_list.size()):
            domain_parts.append(self.main_model_part.GetSubModelPart(self.problem_domain_sub_model_part_list[i].GetString()))
        # Adding Nodes to Computing Model Part
        list_of_ids = set()
        for part in domain_parts:
            for node in part.Nodes:
                list_of_ids.add(node.Id)
        fluid_transport_computing_model_part.AddNodes(list(list_of_ids))
        # Adding Elements to Computing Model Part
        list_of_ids = set()
        for part in domain_parts:
            for elem in part.Elements:
                list_of_ids.add(elem.Id)
        fluid_transport_computing_model_part.AddElements(list(list_of_ids))
        # Adding Conditions to Computing Model Part
        domain_conditions = []
        for i in range(self.processes_sub_model_part_list.size()):
            domain_conditions.append(self.main_model_part.GetSubModelPart(self.processes_sub_model_part_list[i].GetString()))
        list_of_ids = set()
        for part in domain_conditions:
            for cond in part.Conditions:
                list_of_ids.add(cond.Id)
        fluid_transport_computing_model_part.AddConditions(list(list_of_ids))