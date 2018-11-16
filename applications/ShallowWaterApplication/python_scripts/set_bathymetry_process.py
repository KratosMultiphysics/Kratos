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
                "model_part_name"      : "please_specify_model_part_name",
                "interval"             : [0.0, 1e30],
                "variable_name"        : "BATYMETRY",
                "constrained"          : false,
                "value"                : "z"
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        import assign_scalar_variable_process
        self.process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess2D(
            self.model_part, Shallow.BATHYMETRY, Shallow.TOPOGRAPHY_GRADIENT, KratosMultiphysics.NODAL_AREA
        )
