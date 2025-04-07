import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFluidTopologyOptimizationFixVelocityComponentProcess(Model, settings["Parameters"])


class ApplyFluidTopologyOptimizationFixVelocityComponentProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)
        
        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "interval"        : [0.0,1e30],
                "component"       : "Z",
                "value_physics"   : 0.0,
                "value_adjoint"   : 0.0,
                "constrained_physics" : true,
                "constrained_adjoint" : true
            }
            """)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        if (settings["component"].GetString().upper() not in ["X", "Y", "Z"]):
            raise Exception("NOt valid component to fix.")
        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()
        component = settings["component"].GetString().upper()
        value_physics = settings["value_physics"].GetDouble()
        value_adjoint = settings["value_adjoint"].GetDouble()
        constrained_physics = settings["constrained_physics"].GetBool()
        constrained_adjoint = settings["constrained_adjoint"].GetBool()

        velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "VELOCITY_""" + component + """\",
                "interval"             : """ + interval + """,
                "value"                : """ + str(value_physics) + """,
                "constrained"          : """ + str(constrained_physics).lower() + """,
                "local_axes"           : {}
            }
            """
            )

        adjoint_velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "VELOCITY_ADJ_""" + component + """\",
                "interval"             : """ + interval + """,
                "value"                : """ + str(value_adjoint) + """,
                "constrained"          : """ + str(constrained_adjoint).lower() + """,
                "local_axes"           : {}
            }
            """
            )
    
        self.velocity_component_process = AssignScalarVariableProcess(Model, velocity_settings)
        self.adjoint_velocity_component_process = AssignScalarVariableProcess(Model, adjoint_velocity_settings)


    def ExecuteInitializeSolutionStep(self):
        self.velocity_component_process.ExecuteInitializeSolutionStep()
        self.adjoint_velocity_component_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        self.velocity_component_process.ExecuteFinalizeSolutionStep()
        self.adjoint_velocity_component_process.ExecuteFinalizeSolutionStep()
