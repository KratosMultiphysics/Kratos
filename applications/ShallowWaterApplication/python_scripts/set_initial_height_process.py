import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInitialHeightProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class SetInitialHeightProcess(KratosMultiphysics.Process):

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

        if settings["variable_name"].GetString() == "FREE_SURFACE_ELEVATION":
            time_unit_converter = Model["main_model_part"].ProcessInfo.GetValue(KratosShallow.TIME_UNIT_CONVERTER)
            free_surface = settings["value"].GetString()
            settings["value"].SetString(free_surface + '-z*' + str(time_unit_converter))
            settings["variable_name"].SetString("HEIGHT")

        import assign_scalar_variable_process
        self.process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
