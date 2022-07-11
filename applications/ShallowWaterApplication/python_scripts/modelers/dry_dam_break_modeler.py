# importing the Kratos Library
import KratosMultiphysics as KM

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.benchmarks.dry_dam_break_benchmark import DryDamBreakBenchmark

def Factory(model, settings):
    return DryDamBreakModeler(model, settings)

class DryDamBreakModeler(KM.Modeler):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.process = DryDamBreakBenchmark(model, settings)

    def PrepareGeometryModel(self):
        self.process.ExecuteInitialize()
        self.process.Check()
