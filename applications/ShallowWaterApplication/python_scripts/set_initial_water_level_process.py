import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInitialWaterLevelProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class SetInitialWaterLevelProcess(KM.Process):

    def __init__(self, Model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "interval"             : [0.0, 1e30],
                "variable_name"        : "HEIGHT",
                "constrained"          : false,
                "value"                : "1.0"
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = settings["variable_name"].GetString()
        self.model_part = Model[settings["model_part_name"].GetString()]

        from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
        self.process = AssignScalarVariableProcess(Model, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
        if self.variable == "HEIGHT":
            SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.model_part)
        elif self.variable == "FREE_SURFACE_ELEVATION":
            SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.model_part)
