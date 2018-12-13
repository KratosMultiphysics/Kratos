import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VisualizationMeshProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class VisualizationMeshProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name"
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        self.process = KratosShallow.ShallowWaterVariablesUtility(Model[settings["model_part_name"].GetString()])

        # The DefineDryProperties methods duplicates the current number of properties:
        # For each property, it creates another one, which means dry stateself.
        # It should be called only once, otherwise, the number of properties will increase without limits
        self.process.DefineDryProperties()

    def ExecuteBeforeOutputStep(self):
        # The elements should be active to be included on the GidOutputProcess
        self.process.SetElementsActive()
        self.process.AssignDryWetProperties()
        self.process.SetMeshPosition()

    def ExecuteAfterOutputStep(self):
        self.process.ResetMeshPosition()
