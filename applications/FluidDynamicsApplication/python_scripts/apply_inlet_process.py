import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

from KratosMultiphysics import assign_vector_by_direction_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyInletProcess(Model, settings["Parameters"])


class ApplyInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"         : 0,
            "model_part_name" : "",
            "variable_name"   : "VELOCITY",
            "modulus"         : 0.0,
            "constrained"     : true,
            "direction"       : [1.0,0.0,0.0],
            "interval"        : [0.0,"End"]
        }
        """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if (settings.Has("modulus")):
            if (settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")

        if (settings.Has("direction")):
            if (settings["direction"].IsString()):
                default_settings["direction"].SetString("automatic_inwards_normal")

        settings.ValidateAndAssignDefaults(default_settings)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        elif (settings["variable_name"].GetString() != "VELOCITY"):
            raise Exception("Inlet variable_name is not VELOCITY.")
        else:
            if (settings["modulus"].IsDouble()):
                if (settings["modulus"].GetDouble == 0.0):
                    raise Exception("Inlet scalar value equal to 0.")
            elif (settings["modulus"].IsString()):
                if (settings["modulus"].GetString == ""):
                    raise Exception("Inlet function sting is empty.")

        # Set the INLET flag in the inlet model part nodes and conditions
        self.inlet_model_part = Model[settings["model_part_name"].GetString()]
        for node in self.inlet_model_part.Nodes:
            node.Set(KratosMultiphysics.INLET, True)
        for condition in self.inlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.INLET, True)

        # Construct the base process AssignVectorByDirectionProcess
        self.aux_process = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_process.ExecuteFinalizeSolutionStep()
