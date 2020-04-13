import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyWallFunctionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyWallFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "avoid_recomputing_normals": false
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.avoid_recomputing_normals = settings["avoid_recomputing_normals"].GetBool()

        for node in self.model_part.Nodes:
            node.Set(KratosMultiphysics.SLIP, True)
            node.Set(KratosMultiphysics.STRUCTURE, True)

        for condition in self.model_part.Conditions:
            condition.Set(KratosMultiphysics.SLIP, True)
            condition.Set(KratosMultiphysics.STRUCTURE, True)

