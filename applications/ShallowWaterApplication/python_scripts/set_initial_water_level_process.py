import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInitialWaterLevelProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class SetInitialWaterLevelProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
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

        import assign_scalar_variable_process
        self.process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, settings)
        self.variables_utility = KratosShallow.ShallowWaterVariablesUtility(self.model_part)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
        if self.variable == "HEIGHT":
            self.variables_utility.ComputeFreeSurfaceElevation()
        elif self.variable == "FREE_SURFACE_ELEVATION":
            self.variables_utility.ComputeHeightFromFreeSurface()
