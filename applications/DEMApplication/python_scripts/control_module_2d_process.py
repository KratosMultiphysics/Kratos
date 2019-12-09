import KratosMultiphysics
import KratosMultiphysics.DEMApplication as Dem

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ControlModule2DProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ControlModule2DProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        # Control module process acting on the imposed direction: 0 (X), 1 (Y), 2 (Z) or 3 (radial)
        # The radial direction is valid only for the vertical walls of a right cylinder with the base
        # on the 'X-Y' plane centered on (0,0). Negative target_stress means compression.

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.params = KratosMultiphysics.Parameters("{}")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.params.AddValue("imposed_direction",settings["imposed_direction"])
        self.params.AddValue("target_stress_table_id",settings["target_stress_table_id"])
        self.params.AddValue("initial_velocity",settings["initial_velocity"])
        self.params.AddValue("limit_velocity",settings["limit_velocity"])
        self.params.AddValue("velocity_factor",settings["velocity_factor"])
        self.params.AddValue("compression_length",settings["compression_length"])
        self.params.AddValue("young_modulus",settings["young_modulus"])
        self.params.AddValue("stress_increment_tolerance",settings["stress_increment_tolerance"])
        self.params.AddValue("update_stiffness",settings["update_stiffness"])
        self.params.AddValue("start_time",settings["start_time"])

        self.control_module_process = Dem.ControlModule2DProcess(self.model_part, self.params)

    def ExecuteInitialize(self):
        self.control_module_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.control_module_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.control_module_process.ExecuteFinalizeSolutionStep()