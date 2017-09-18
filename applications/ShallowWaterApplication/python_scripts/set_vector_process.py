import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetVectorProcess(Model, settings["Parameters"])

## This process sets the value of a vector variable using the AssignVectorVariableProcess.
class SetVectorProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "VELOCITY",
                "interval"             : [0.0, 1e30],
                "value"                : [10.0, 0.0, 0.0],
                "constrained"          : [false,false,false],
                "local_axes"           : {}
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        import assign_vector_variable_process
        self.process = assign_vector_variable_process.AssignVectorVariableProcess(Model, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
