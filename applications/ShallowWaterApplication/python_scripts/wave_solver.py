# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return WaveSolver(model, custom_settings)

class WaveSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.element_name, self.condition_name, self.min_buffer_size = self._GetFormulationSettings()

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)
        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water primitive DOFs added correctly.")

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.main_model_part)

    def _SetProcessInfo(self):
        super()._SetProcessInfo()
        self.main_model_part.ProcessInfo.SetValue(KM.STABILIZATION_FACTOR, self.settings["stabilization_factor"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(SW.RELATIVE_DRY_HEIGHT, self.settings["relative_dry_height"].GetDouble())

    def _GetFormulationSettings(self):
        scheme = self.settings["time_integration_scheme"].GetString()
        order = self.settings["time_integration_order"].GetInt()
        if scheme == "bdf":
            element_name = "WaveElement"
            condition_name = "WaveCondition"
            buffer_size = order + 1
        elif scheme == "crank_nicolson":
            element_name = "CrankNicolsonWaveElement"
            condition_name = "WaveCondition"
            buffer_size = 2
            if not order == 2:
                KM.Logger.PrintWarning('WaveSolver', 'Setting the "time_integration_order" to 2')
        else:
            raise Exception('The possible "time_integration_scheme" are "bdf" and "crank_nicolson"')
        return element_name, condition_name, buffer_size

    def _CreateScheme(self):
        scheme = self.settings["time_integration_scheme"].GetString()
        if scheme == "bdf":
            scheme_settings = KM.Parameters()
            scheme_settings.AddStringArray("solution_variables", ["VELOCITY","HEIGHT"])
            scheme_settings.AddValue("integration_order", self.settings["time_integration_order"])
            time_scheme = SW.ShallowWaterResidualBasedBDFScheme(scheme_settings)
        else:
            time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        return time_scheme

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""{
            "time_integration_scheme"    : "bdf",
            "time_integration_order"     : 2,
            "relative_dry_height"        : 0.1,
            "stabilization_factor"       : 0.01
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings
