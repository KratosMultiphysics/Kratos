from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication.cellular_flow.ethier_benchmark_analysis as ethier_benchmark_analysis
BaseAnalysis = ethier_benchmark_analysis.EthierBenchmarkAnalysis

class EthierBenchmarkMakeMeshAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        self.project_parameters.AddEmptyValue("pressure_grad_recovery_type")
        self.project_parameters.AddEmptyValue("size_parameter").SetInt(1)

    def ReadFluidModelParts(self):
        from meshing import meshing_utilities
        self.mesh_generator = meshing_utilities.ParallelepipedRegularMesher(
                                                model_part_to_be_filled = self._GetFluidAnalysis().fluid_model_part,
                                                lower_corner_coordinates = [0.0, 0.0, 0.0],
                                                higher_corner_coordinates = [0.1, 0.1, 0.1],
                                                number_of_divisions_per_dimension = 10)

        self.mesh_generator.FillModelPartWithNewMesh()
