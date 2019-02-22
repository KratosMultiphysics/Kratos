import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetBathymetryProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class SetBathymetryProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "",
                "interval"             : [0.0, 1e30],
                "variable_name"        : "BATHYMETRY",
                "constrained"          : false,
                "value"                : "z"
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        import assign_scalar_variable_process
        # Data process
        self.variable_process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, settings)
        # Z-coordinate process
        self.coordinate_process = Shallow.ShallowWaterVariablesUtility(Model[settings["model_part_name"].GetString()])

    def ExecuteInitialize(self):
        # Data process: define the topography variable
        self.variable_process.ExecuteInitializeSolutionStep()
        # Z-coordinate process: assign the topography variable to the z coordinate
        self.coordinate_process.SetMeshPosition(True)
