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

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "volume_model_part_name" : "",
            "skin_parts" : [],
            "assign_neighbour_elements_to_conditions" : false
        }''')
        Parameters.ValidateAndAssignDefaults(default_parameters)
        if Parameters["volume_model_part_name"].GetString() == "":
            raise Exception("Please define the \"volume_model_part_name\" (string) argument.")

        self.volume_model_part_name = Parameters["volume_model_part_name"].GetString()
        self.skin_name_list = Parameters["skin_parts"]

        self.assign_neighbour_elements = Parameters["assign_neighbour_elements_to_conditions"].GetBool()


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
        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = tmoc.NOT_COMPUTE_NODAL_NORMALS | tmoc.NOT_COMPUTE_CONDITION_NORMALS
        if self.assign_neighbour_elements:
            flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        else:
            flags |= tmoc.NOT_ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        KratosMultiphysics.TetrahedralMeshOrientationCheck(fluid_computational_model_part,throw_errors, flags).Execute()

