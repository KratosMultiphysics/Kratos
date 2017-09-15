import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetVectorProcess(Model, settings["Parameters"])

## This process sets the value of a vector variable using the ApplyConstantVectorValueProcess.
class SetVectorProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "VELOCITY",
                "is_fixed_x": false,
                "is_fixed_y": false,
                "is_fixed_z": false,
                "modulus" : 1.0,
                "direction": [1.0, 0.0, 0.0]
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        import assign_scalar_variable_process
        self.process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
