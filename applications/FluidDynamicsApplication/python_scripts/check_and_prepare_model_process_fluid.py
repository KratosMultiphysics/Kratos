import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    def __init__(self, main_model_part, Parameters ):
        self.main_model_part = main_model_part

        self.volume_model_part_name = Parameters["volume_model_part_name"].GetString()
        self.skin_name_list = Parameters["skin_parts"]


        #self.volume_model_part_name = Parameters["volume_model_part_name"].GetString()
        #self.list_of_inlets = Parameters["list_of_inlets"]
        #self.list_of_slip = Parameters["list_of_inlets"]
        #self.list_of_inlets = Parameters["list_of_inlets"]



    def Execute(self):
        self.volume_model_part = self.main_model_part.GetSubModelPart(self.volume_model_part_name)

        skin_parts = []
        for i in range(self.skin_name_list.size()):
            skin_parts.append(self.main_model_part.GetSubModelPart(self.skin_name_list[i].GetString()))

        #construct a model part which contains both the skin and the volume
        #temporarily we call it "fluid_computational_model_part"
        self.main_model_part.CreateSubModelPart("fluid_computational_model_part")
        fluid_computational_model_part= self.main_model_part.GetSubModelPart("fluid_computational_model_part")
        fluid_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo

        for node in self.volume_model_part.Nodes:
            fluid_computational_model_part.AddNode(node,0)
        for elem in self.volume_model_part.Elements:
            fluid_computational_model_part.AddElement(elem,0)

        #do some gymnastics to have this done fast. - create an ordered list to be added
        list_of_ids = set()
        for part in skin_parts:
            for cond in part.Conditions:
                list_of_ids.add(cond.Id)

        fluid_computational_model_part.AddConditions(list(list_of_ids))

        #verify the orientation of the skin
        throw_errors = False
        KratosMultiphysics.TetrahedralMeshOrientationCheck(fluid_computational_model_part,throw_errors).Execute()

