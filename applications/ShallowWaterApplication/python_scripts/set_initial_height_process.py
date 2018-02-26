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

        self.variable = settings["variable_name"].GetString()
        self.model_part = Model[settings["model_part_name"].GetString()]
        #~ if settings["variable_name"].GetString() == "FREE_SURFACE_ELEVATION":
            #~ height_unit_converter = Model["main_model_part"].ProcessInfo.GetValue(KratosShallow.WATER_HEIGHT_UNIT_CONVERTER)
            #~ free_surface = settings["value"].GetString()
            #~ settings["value"].SetString(free_surface + '-z*' + str(height_unit_converter))
            #~ settings["variable_name"].SetString("HEIGHT")

        import assign_scalar_variable_process
        self.process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
        if self.variable == "HEIGHT":
            for node in self.model_part.Nodes:
                free_surface = node.GetSolutionStepValue(KratosShallow.HEIGHT) + node.GetSolutionStepValue(KratosShallow.BATHYMETRY)
                node.SetSolutionStepValue(KratosShallow.FREE_SURFACE_ELEVATION,0,free_surface)
        elif self.variable == "FREE_SURFACE_ELEVATION":
            for node in self.model_part.Nodes:
                height = node.GetSolutionStepValue(KratosShallow.FREE_SURFACE_ELEVATION) - node.GetSolutionStepValue(KratosShallow.BATHYMETRY)
                node.SetSolutionStepValue(KratosShallow.HEIGHT,0,height)
