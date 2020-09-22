# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return ShallowWaterSolver(model, custom_settings)

class ShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "SWE"
        self.condition_name = "LineCondition"
        self.min_buffer_size = 2
        self.advection_epsilon = self.settings["advection_epsilon"].GetDouble()

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.FREE_SURFACE_ELEVATION, self.main_model_part)

        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water solver DOFs added correctly.")

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        epsilon = max(self.advection_epsilon, self.GetComputingModelPart().ProcessInfo[SW.DRY_HEIGHT])
        SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.GetComputingModelPart())
        SW.ComputeVelocityProcess(self.GetComputingModelPart(), epsilon).Execute()
        SW.ShallowWaterUtilities().ComputeAccelerations(self.GetComputingModelPart())

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
            "advection_epsilon"     : 1.0e-2,
            "permeability"          : 1.0e-4,
            "dry_discharge_penalty" : 1.0e+2
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        default_settings["wetting_drying_model"] = KM.Parameters("""
        {
            "model_name" : "negative_height",
            "beta" : 1e4
        }""")
        return default_settings

    def PrepareModelPart(self):
        super().PrepareModelPart()
        permeability = self.settings["permeability"].GetDouble()
        discharge_penalty = self.settings["dry_discharge_penalty"].GetDouble()
        if permeability == 0.0:
            KM.Logger.PrintWarning(self.__class__.__name__, "Detected permeability == 0.0")
        if discharge_penalty == 0.0:
            KM.Logger.PrintWarning(self.__class__.__name__, "Detected dry_discharge_penalty == 0.0")
        self.main_model_part.ProcessInfo.SetValue(SW.PERMEABILITY, permeability)
        self.main_model_part.ProcessInfo.SetValue(SW.DRY_DISCHARGE_PENALTY, discharge_penalty)
