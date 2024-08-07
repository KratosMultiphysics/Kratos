import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis
import numpy as np

class GeoMechanicsAnalysisGatheringNonConverged(GeoMechanicsAnalysis):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

    def Initialize(self):
        super().Initialize()
        self._GetSolver()._GetSolutionStrategy().SetUpNonconvergedSolutionsFlag(True)

    def Finalize(self):
        self.kratos_matrix_nonconverged, self.dofs_before_clear = self._GetSolver()._GetSolutionStrategy().GetNonconvergedSolutions()
        super().Finalize()

    def GetNonconvergedSolutions(self):
        return np.array(self.kratos_matrix_nonconverged, copy=False), self.dofs_before_clear


class TestNonconvergedSolutions(KratosUnittest.TestCase):

    def test_nonconverged_solutions_netwton_raphson_erosion_process(self):
        self.work_folder = "test_compare_sellmeijer/HeightAquiferD10L30.gid"
        parameters_filename = "ProjectParameters.json"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            with open(parameters_filename, 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = GeoMechanicsAnalysisGatheringNonConverged(model,parameters)

            # Run test case
            simulation.Run()

            # Check results
            np_array_with_nonconverged_solutions, dofs = simulation.GetNonconvergedSolutions()
            num_nl_iters = simulation._GetSolver().GetComputingModelPart().GetRootModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]
            self.assertEqual(np_array_with_nonconverged_solutions.shape[1], num_nl_iters+1)
            self.assertEqual(np_array_with_nonconverged_solutions.shape[0], len(dofs))


if __name__ == '__main__':
    KratosUnittest.main()