import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetTopographyProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class SetTopographyProcess(KM.Process):

    def __init__(self, Model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "interval"             : [0.0, 1e30],
                "constrained"          : false,
                "value"                : "z"
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)
        settings.AddEmptyValue("variable_name").SetString("TOPOGRAPHY")

        self.model_part = Model[settings["model_part_name"].GetString()]

        from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
        self.process = AssignScalarVariableProcess(Model, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
        SW.ShallowWaterUtilities().FlipScalarVariable(SW.TOPOGRAPHY, SW.BATHYMETRY, self.model_part)
