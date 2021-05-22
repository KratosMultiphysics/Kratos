# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return LagrangianShallowWaterSolver(model, custom_settings)

class LagrangianShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.min_buffer_size = 2
        self.element_name = "WaveElement"
        self.condition_name = "LineCondition"

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)
        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water primitive DOFs added correctly.")

    def Initialize(self):
        super().Initialize()
        self.main_model_part.ProcessInfo.SetValue(KM.STABILIZATION_FACTOR, self.settings["stabilization_factor"].GetDouble())

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def _CreateScheme(self):
        order = self.settings["time_integration_order"].GetInt()
        time_scheme = SW.VelocityHeightResidualBasedBDFScheme(order)
        return time_scheme

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""{
            "time_integration_order"     : 2,
            "stabilization_factor"       : 0.01,
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings
