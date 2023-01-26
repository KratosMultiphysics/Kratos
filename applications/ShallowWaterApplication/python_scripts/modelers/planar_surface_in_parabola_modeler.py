# importing the Kratos Library
import KratosMultiphysics as KM

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.benchmarks.planar_surface_in_parabola_benchmark import PlanarSurfaceInParabolaBenchmark

def Factory(model, settings):
    return PlanarSurfaceInParabolaModeler(model, settings)

class PlanarSurfaceInParabolaModeler(KM.Modeler):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.process = PlanarSurfaceInParabolaBenchmark(model, settings)

    def PrepareGeometryModel(self):
        self.process.ExecuteInitialize()
        self.process.Check()
