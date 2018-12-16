from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import ethier_benchmark_analysis
BaseAnalysis = ethier_benchmark_analysis.EthierBenchmarkAnalysis

class ProductOfSinesBenchmarkAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        Add = self.pp.CFD_DEM.AddEmptyValue
        Add("field_period").SetDouble(1.0)

    def GetFieldUtility(self):
        period = self.pp.CFD_DEM["field_period"].GetDouble()
        self.flow_field = ProductOfSines(period)
        space_time_set = SpaceTimeSet()
        self.field_utility = FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility
