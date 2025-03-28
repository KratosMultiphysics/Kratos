import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTransportTopologyOptimizationDirichletProcess(Model, settings["Parameters"])


class ApplyTransportTopologyOptimizationDirichletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        KratosMultiphysics.Process.__init__(self)
        
        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "interval"        : [0.0,1e30],
                "value"           : 0.0
            }
            """)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()
        value = settings["value"].GetString().upper()

        physics_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "TEMPERATURE",
                "interval"             : """ + interval + """,
                "value"                : """ + value + """,
                "constrained"          : true,
                "local_axes"           : {}
            }
            """
            )

        adjoint_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "TEMPERATURE_ADJ",
                "interval"             : """ + interval + """,
                "value"                : """ + value + """,
                "constrained"          : true,
                "local_axes"           : {}
            }
            """
            )
    
        self.physics_process = AssignScalarVariableProcess(Model, physics_settings)
        self.adjoint_process = AssignScalarVariableProcess(Model, adjoint_settings)


    def ExecuteInitializeSolutionStep(self):
        self.physics_process.ExecuteInitializeSolutionStep()
        self.adjoint_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        self.physics_process.ExecuteFinalizeSolutionStep()
        self.adjoint_process.ExecuteFinalizeSolutionStep()