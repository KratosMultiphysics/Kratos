import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyAbsorbingBoundaryProcess(model, settings["Parameters"])

class ApplyAbsorbingBoundaryProcess(KM.Process):
    """ApplyAbsorbingBoundaryProcess

    This process sets the DISTANCE variable from every
    node to the nearest boundary condition
    """

    def __init__(self, model, settings):
        """The constructor of the ApplyAbsorbingBoundaryProcess

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The model to be used
        settings -- The ProjectParameters used
        """

        KM.Process.__init__(self)

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.processes = []
        self.model_part = model[self.settings["computing_model_part_name"].GetString()]
        boundaries_names = self.settings["absorbing_boundaries_list"].GetStringArray()
        distance_calculator_settings = KM.Parameters()
        distance_calculator_settings.AddValue("r_squared_threshold", self.settings["r_squared_threshold"])
        for name in boundaries_names:
            boundary_part = self.model_part.GetSubModelPart(name)
            self.processes.append(SW.CalculateDistanceToBoundaryProcess(self.model_part, boundary_part, distance_calculator_settings))

        for process in self.processes:
            process.Check()

    def Check(self):
        direction = self.settings["velocity_direction"].GetVector()
        if not direction.Size() == 3:
            raise Exception("The boundary velocity direction must be specified with a three dimensional array")
        if direction.norm_2() == 0.0:
            raise Exception("The boundary velocity direction has zero norm")

    def ExecuteInitialize(self):
        KM.VariableUtils().SetVariable(KM.DISTANCE, 1e+38, self.model_part.Nodes)

        if self.settings["apply_fixity"].GetBool():
            boundaries_names = self.settings["absorbing_boundaries_list"].GetStringArray()
            for name in boundaries_names:
                boundary_part = self.model_part.GetSubModelPart(name)
                KM.VariableUtils().ApplyFixity(KM.MOMENTUM_X, True, boundary_part.Nodes)
                KM.VariableUtils().ApplyFixity(KM.MOMENTUM_Y, True, boundary_part.Nodes)

    def ExecuteBeforeSolutionLoop(self):
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.model_part.ProcessInfo.SetValue(SW.ABSORBING_DISTANCE, self.settings["absorbing_distance"].GetDouble())
        self.model_part.ProcessInfo.SetValue(SW.DISSIPATION, self.settings["dissipation_factor"].GetDouble())
        direction = self.settings["velocity_direction"].GetVector()
        direction /= direction.norm_2()
        boundary_velocity = self.settings["velocity_modulus"].GetDouble() * direction
        self.model_part.ProcessInfo.SetValue(SW.BOUNDARY_VELOCITY, boundary_velocity)

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
                "computing_model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "absorbing_boundaries_list" : [],
                "r_squared_threshold"       : 0.99,
                "absorbing_distance"        : 0.0,
                "dissipation_factor"        : 0.0,
                "velocity_modulus"          : 0.0,
                "velocity_direction"        : [1.0, 0.0, 0.0],
                "apply_fixity"              : false
            }""")
