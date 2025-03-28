import KratosMultiphysics
from KratosMultiphysics.assign_vector_variable_process import AssignVectorVariableProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFluidTopologyOptimizationNoSlipProcess(Model, settings["Parameters"])


class ApplyFluidTopologyOptimizationNoSlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)
        
        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "interval"        : [0.0,1e30]
            }
            """)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()

        velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "VELOCITY",
                "interval"             : """ + interval + """,
                "value"                : [0.0, 0.0, 0.0],
                "constrained"          : [true,true,true],
                "local_axes"           : {}
            }
            """
            )

        adjoint_velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "VELOCITY_ADJ",
                "interval"             : """ + interval + """,
                "value"                : [0.0, 0.0, 0.0],
                "constrained"          : [true,true,true],
                "local_axes"           : {}
            }
            """
            )
    
        self.velocity_vector_process = AssignVectorVariableProcess(Model, velocity_settings)
        self.adjoint_velocity_vector_process = AssignVectorVariableProcess(Model, adjoint_velocity_settings)


    def ExecuteInitializeSolutionStep(self):
        self.velocity_vector_process.ExecuteInitializeSolutionStep()
        self.adjoint_velocity_vector_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        self.velocity_vector_process.ExecuteFinalizeSolutionStep()
        self.adjoint_velocity_vector_process.ExecuteFinalizeSolutionStep()
