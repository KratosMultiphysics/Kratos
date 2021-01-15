# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.stabilized_shallow_water_solver import StabilizedShallowWaterSolver

def CreateSolver(model, custom_settings):
    return FluxCorrectedShallowWaterSolver(model, custom_settings)

class FluxCorrectedShallowWaterSolver(StabilizedShallowWaterSolver):

    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.element_name = "MonotonicElement"

    def _CreateScheme(self):
        time_scheme = SW.FluxCorrectedShallowWaterScheme(self.settings["time_integration_order"].GetInt())
        return time_scheme
