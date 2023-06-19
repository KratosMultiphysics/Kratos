# importing the Kratos Library
import KratosMultiphysics as KM

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.benchmarks.dam_break_benchmark import DamBreakBenchmark

def Factory(model, settings):
    return DamBreakModeler(model, settings)

class DamBreakModeler(KM.Modeler):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.process = DamBreakBenchmark(model, settings)

    def PrepareGeometryModel(self):
        self.process.ExecuteInitialize()
        self.process.Check()
