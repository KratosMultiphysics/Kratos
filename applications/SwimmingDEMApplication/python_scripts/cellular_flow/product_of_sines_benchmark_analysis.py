from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import KratosMultiphysics.SwimmingDEMApplication.cellular_flow.ethier_benchmark_analysis as ethier_benchmark_analysis
BaseAnalysis = ethier_benchmark_analysis.EthierBenchmarkAnalysis

class ProductOfSinesBenchmarkAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        Add = self.project_parameters.AddEmptyValue
        Add("field_period").SetDouble(1.0)

    def GetFieldUtility(self):
        period = self.project_parameters["field_period"].GetDouble()
        self.flow_field = SDEM.ProductOfSines(period)
        space_time_set = SDEM.SpaceTimeSet()
        self.field_utility = SDEM.FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility
