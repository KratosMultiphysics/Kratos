# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.benchmarks.planar_surface_in_parabola_benchmark import PlanarSurfaceInParabolaBenchmark

def Factory(model, settings):
    return PlanarSurfaceInParabolaModeler(model, settings)

class PlanarSurfaceInParabolaModeler(KM.Modeler):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        process_settings = KM.Parameters("""{}""")
        process_settings.AddValue("Parameters", settings)
        self.process = PlanarSurfaceInParabolaBenchmark(model, process_settings)

    def PrepareGeometryModel(self):
        self.process.Check()
        self.process.ExecuteInitialize()
