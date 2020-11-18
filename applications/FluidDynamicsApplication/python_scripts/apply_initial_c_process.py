import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication 
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyInitialCProcess(Model, settings["Parameters"])


class ApplyInitialCProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "SOUND_VELOCITY",
                "interval"             : [0.0, 1e30],
                "value"                : 1e12,
                "constrained"          : false,
                "local_axes"           : {}
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.scalar_process = AssignScalarVariableProcess(Model, settings)
        self.scalar_process.ExecuteInitializeSolutionStep()
        self.scalar_process.ExecuteFinalizeSolutionStep()

        
