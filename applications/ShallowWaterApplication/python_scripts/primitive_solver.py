# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.wave_solver import WaveSolver

def CreateSolver(model, custom_settings):
    return PrimitiveSolver(model, custom_settings)

class PrimitiveSolver(WaveSolver):

    def _GetFormulationSettings(self):
        order = self.settings["time_integration_order"].GetInt()
        element_name = "PrimitiveElement"
        condition_name = "PrimitiveCondition"
        buffer_size = order + 1
        return element_name, condition_name, buffer_size
