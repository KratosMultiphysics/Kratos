import KratosMultiphysics
from KratosMultiphysics.assign_vector_variable_process import AssignVectorVariableProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyNoSlipProcess(Model, settings["Parameters"])


class ApplyNoSlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "VELOCITY",
                "interval"             : [0.0, 1e30],
                "value"                : [0.0, 0.0, 0.0],
                "constrained"          : [true,true,true],
                "local_axes"           : {}
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.vector_process = AssignVectorVariableProcess(Model, settings)


    def ExecuteInitializeSolutionStep(self):
        self.vector_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        self.vector_process.ExecuteFinalizeSolutionStep()
