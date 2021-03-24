# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return FreeSurfaceShallowWaterSolver(model, custom_settings)

class FreeSurfaceShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "SWE"
        self.condition_name = "LineCondition"
        self.min_buffer_size = 2

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.FREE_SURFACE_ELEVATION, self.main_model_part)

        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water solver DOFs added correctly.")

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.GetComputingModelPart())
        SW.ShallowWaterUtilities().ComputeVelocity(self.GetComputingModelPart(), True)
        SW.ShallowWaterUtilities().ComputeAccelerations(self.GetComputingModelPart())

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
            "permeability"          : 1.0e-4,
            "dry_height_threshold"  : 1e-3,
            "dry_discharge_penalty" : 1.0e+2,
            "stabilization_factor"  : 0.005,
            "relative_dry_height"   : 0.1
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def PrepareModelPart(self):
        super().PrepareModelPart()
        permeability = self.settings["permeability"].GetDouble()
        dry_height = self.settings["dry_height_threshold"].GetDouble()
        relative_dry_height = self.settings["relative_dry_height"].GetDouble()
        discharge_penalty = self.settings["dry_discharge_penalty"].GetDouble()
        stabilization_factor = self.settings["stabilization_factor"].GetDouble()

        # Checks
        if permeability == 0.0:
            KM.Logger.PrintWarning(self.__class__.__name__, "Detected permeability == 0.0")
        if discharge_penalty == 0.0:
            KM.Logger.PrintWarning(self.__class__.__name__, "Detected dry_discharge_penalty == 0.0")

        self.main_model_part.ProcessInfo.SetValue(SW.PERMEABILITY, permeability)
        self.main_model_part.ProcessInfo.SetValue(SW.DRY_HEIGHT, dry_height)
        self.main_model_part.ProcessInfo.SetValue(SW.RELATIVE_DRY_HEIGHT, relative_dry_height)
        self.main_model_part.ProcessInfo.SetValue(SW.DRY_DISCHARGE_PENALTY, discharge_penalty)
        self.main_model_part.ProcessInfo.SetValue(KM.STABILIZATION_FACTOR, stabilization_factor)

    def _CreateScheme(self):
        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        return time_scheme
