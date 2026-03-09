import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.ShallowWaterApplication.utilities.wave_factory import WaveTheoryFactory

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyAbsorbingBoundaryProcess(model, settings["Parameters"])

class ApplyAbsorbingBoundaryProcess(KM.Process):
    """ApplyAbsorbingBoundaryProcess

    This process sets the DISTANCE variable from every
    node to the nearest boundary condition
    """

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "computing_model_part_name" : "",
            "absorbing_boundary_name"   : "",
            "r_squared_threshold"       : 0.99,
            "relative_distance"         : 2.0,
            "relative_damping"          : 2.0,
            "wave_specifications"       : {}
        }""")

    def __init__(self, model, settings):
        """The constructor of the ApplyAbsorbingBoundaryProcess"""

        KM.Process.__init__(self)

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model.GetModelPart(self.settings["computing_model_part_name"].GetString())
        self.boundary_part = model.GetModelPart(self.settings["absorbing_boundary_name"].GetString())
        self.distance_process = SW.CalculateDistanceToBoundaryProcess(self.model_part, self.boundary_part, self.settings["r_squared_threshold"].GetDouble())

        variables_names = KM.SpecificationsUtilities.GetDofsListFromConditionsSpecifications(self.boundary_part)
        self.variables_to_fix = []
        for variable_name in variables_names:
            if variable_name.startswith("VELOCITY") or variable_name.startswith("MOMENTUM"):
                variable = KM.KratosGlobals.GetVariable(variable_name)
                self.variables_to_fix.append(variable)

    def Check(self):
        """Check the correctness of the input."""
        self.distance_process.Check()

    def ExecuteBeforeSolutionLoop(self):
        """Calculate the distances and the damping parameters."""
        self.distance_process.ExecuteBeforeSolutionLoop()

        wave = WaveTheoryFactory(self.boundary_part, self.settings["wave_specifications"])
        absorbing_distance = wave.wavelength * self.settings["relative_distance"].GetDouble()
        dissipation_factor = wave.frequency * self.settings["relative_damping"].GetDouble()

        self.model_part.ProcessInfo.SetValue(SW.ABSORBING_DISTANCE, absorbing_distance)
        self.model_part.ProcessInfo.SetValue(SW.DISSIPATION, dissipation_factor)

        for variable in self.variables_to_fix:
            KM.VariableUtils().ApplyFixity(variable, True, self.boundary_part.Nodes)
