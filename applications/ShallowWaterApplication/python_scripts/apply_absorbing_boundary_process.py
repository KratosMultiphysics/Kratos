import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from math import pi, sqrt

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyAbsorbingBoundaryProcess(model, settings["Parameters"])

class ApplyAbsorbingBoundaryProcess(KM.Process):
    """ApplyAbsorbingBoundaryProcess

    This process sets the DISTANCE variable from every
    node to the nearest boundary condition
    """

    __formulation = {
        # Json input
        "primitive_variables" : SW.Formulation.PrimitiveVariables,
        "conserved_variables" : SW.Formulation.ConservativeVariables
    }

    def __init__(self, model, settings):
        """The constructor of the ApplyAbsorbingBoundaryProcess"""

        KM.Process.__init__(self)

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.distance_processes = []
        self.model_part = model[self.settings["computing_model_part_name"].GetString()]
        boundaries_names = self.settings["absorbing_boundaries_list"].GetStringArray()
        distance_calculator_settings = KM.Parameters()
        distance_calculator_settings.AddValue("r_squared_threshold", self.settings["r_squared_threshold"])
        for name in boundaries_names:
            boundary_part = self.model_part.GetSubModelPart(name)
            self.distance_processes.append(SW.CalculateDistanceToBoundaryProcess(self.model_part, boundary_part, distance_calculator_settings))

        for process in self.distance_processes:
            process.Check()

        self.formulation = self.__formulation[settings["formulation"].GetString()]

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
                if self.formulation == SW.Formulation.PrimitiveVariables:
                        KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, True, boundary_part.Nodes)
                        KM.VariableUtils().ApplyFixity(KM.VELOCITY_Y, True, boundary_part.Nodes)
                elif self.formulation == SW.Formulation.ConservativeVariables:
                        KM.VariableUtils().ApplyFixity(KM.MOMENTUM_X, True, boundary_part.Nodes)
                        KM.VariableUtils().ApplyFixity(KM.MOMENTUM_Y, True, boundary_part.Nodes)

    def ExecuteBeforeSolutionLoop(self):
        for process in self.distance_processes:
            process.ExecuteBeforeSolutionLoop()
        
        gravity = self.model_part.ProcessInfo.GetValue(KM.GRAVITY_Z)
        height = 0
        num_nodes = 0
        boundaries_names = self.settings["absorbing_boundaries_list"].GetStringArray()
        for name in boundaries_names:
            boundary_part = self.model_part.GetSubModelPart(name)
            height += KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.HEIGHT, boundary_part, 0)
            num_nodes += boundary_part.NumberOfNodes()
        height /= num_nodes
        self.wave_speed = sqrt(gravity * height)
        self.wave_period = self.settings["wave_period"].GetDouble()

    def ExecuteInitializeSolutionStep(self):
        absorbing_distance = self.__Wavelength() * self.settings["relative_distance"].GetDouble()
        dissipation_factor = self.__Frequency() * self.settings["relative_damping"].GetDouble()
        self.model_part.ProcessInfo.SetValue(SW.ABSORBING_DISTANCE, absorbing_distance)
        self.model_part.ProcessInfo.SetValue(SW.DISSIPATION, dissipation_factor)
        direction = self.settings["velocity_direction"].GetVector()
        direction /= direction.norm_2()
        boundary_velocity = self.settings["velocity_modulus"].GetDouble() * direction
        self.model_part.ProcessInfo.SetValue(SW.BOUNDARY_VELOCITY, boundary_velocity)

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "computing_model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
            "absorbing_boundaries_list" : [],
            "interval"                  : [0.0,"End"]
            "r_squared_threshold"       : 0.99,
            "relative_distance"         : 2.0,
            "relative_damping"          : 2.0,
            "wave_period"               : 1.0,
            "velocity_modulus"          : 0.0,
            "velocity_direction"        : [1.0, 0.0, 0.0],
            "formulation"               : "primitive_variables",
            "apply_fixity"              : false
        }""")

    def __Frequency(self):
        return 2 * pi / self.wave_period

    def __Wavelength(self):
        return self.wave_speed * self.wave_period
