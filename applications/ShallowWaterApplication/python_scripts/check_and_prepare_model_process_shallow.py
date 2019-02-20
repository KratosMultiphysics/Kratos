import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    def __init__(self, model, Parameters):

        default_parameters = KratosMultiphysics.Parameters('''{
            "topographic_model_part_name" : "",
            "computing_model_part_name" : "",
            "build_structured_mesh" : false
        }''')
        Parameters.ValidateAndAssignDefaults(default_parameters)

        # Topographic model part: the model part where whe have import the terrain model
        # Computing model part: the model part where to perform the simulation

        self.topographic_model_part = model[Parameters["topographic_model_part_name"].GetString()]
        self.computing_model_part = model[Parameters["computing_model_part_name"].GetString()]
        self.build_structured_mesh = Parameters["build_structured_mesh"].GetBool()

        self.model = model

    def Execute(self):
        skin_parts_names = []
        for sub_model_part in self.topographic_model_part.SubModelParts:
            skin_parts_names.append(sub_model_part.Name)

        # Create the sub model parts on the empty computing model part
        for name in skin_parts_names:
            self.computing_model_part.CreateSubModelPart(name)
        
        # Create the new mesh or create the existing elements
        if self.build_structured_mesh:
            raise Exception("Not yet implemented")
        else:
            # Populate the sub model parts
            Shallow.PrepareShallowModelUtility(self.topographic_model_part, self.computing_model_part).CopyTopographicModelPart()
