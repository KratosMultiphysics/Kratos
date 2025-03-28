import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

from KratosMultiphysics import assign_vector_by_direction_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFluidTopologyOptimizationInletProcess(Model, settings["Parameters"])


class ApplyFluidTopologyOptimizationInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "",
                "velocity_modulus"         : "0.0",
                "adjoint_velocity_modulus" : "0.0",
                "direction"       : "automatic_inwards_normal",
                "interval"        : [0.0,"End"]
            }
            """)

        if (settings.Has("direction")):
            if (not settings["direction"].IsString()):
                raise Exception("Not string_type direction in fluid topology optimization inlet.")
        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        else:
            if (settings["velocity_modulus"].GetString() == ""):
                raise Exception("Fluid topology optimization Inlet velocity modulus is empty.")
            elif (settings["adjoint_velocity_modulus"].GetString() == ""):
                raise Exception("Fluid topology optimization Inlet adjoint velocity modulus is empty.")
        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        velocity_modulus = settings["velocity_modulus"].GetString()
        adjoint_velocity_modulus = settings["adjoint_velocity_modulus"].GetString()
        direction = settings["direction"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()

        velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : """ + str(mesh_id) +  """,
                "model_part_name" : \"""" + model_part_name + """\",
                "variable_name"   : "VELOCITY",
                "modulus"         : \"""" + velocity_modulus + """\",
                "constrained"     : true,
                "direction"       : \"""" + direction + """\",
                "interval"        : """ + interval + """
            }
            """)

        adjoint_velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : """ + str(mesh_id) + """,
                "model_part_name" : \"""" + model_part_name + """\",
                "variable_name"   : "VELOCITY_ADJ",
                "modulus"         : \"""" + adjoint_velocity_modulus + """\",
                "constrained"     : true,
                "direction"       : \"""" + direction + """\",
                "interval"        : """ + interval + """
            }
            """)
            
        # Set the INLET flag in the inlet model part nodes and conditions
        self.inlet_model_part = Model[settings["model_part_name"].GetString()]
        for node in self.inlet_model_part.Nodes:
            node.Set(KratosMultiphysics.INLET, True)
        for condition in self.inlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.INLET, True)

        # Construct the base process AssignVectorByDirectionProcess
        self.aux_process         = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, velocity_settings)
        self.aux_process_adjoint = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, adjoint_velocity_settings)

    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_process.ExecuteInitializeSolutionStep()
        self.aux_process_adjoint.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_process.ExecuteFinalizeSolutionStep()
        self.aux_process_adjoint.ExecuteFinalizeSolutionStep()
