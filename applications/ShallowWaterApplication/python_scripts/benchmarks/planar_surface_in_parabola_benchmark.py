# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess

# Other imports
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PlanarSurfaceInParabola(model, settings["Parameters"])

class PlanarSurfaceInParabolaBenchmark(BaseBenchmarkProcess):
    def __init__(self):
        super(PlanarSurfaceInParabolaBenchmark, self).__init__(model, settings)

        benchmark_default_settings = KM.Parameters("""
            {
            }
            """
            )
        self.benchmark_settings.ValidateAndAssignDefaults(benchmark_default_settings)
