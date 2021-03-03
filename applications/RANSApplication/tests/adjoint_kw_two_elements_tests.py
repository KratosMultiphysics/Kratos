import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis

from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ReadParameters
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import DeleteH5Files
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ComputeAdjointSensitivity
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis

def SolvePrimalProblem(kratos_parameters):
    test = RANSAnalysis(Kratos.Model(), kratos_parameters)
    test.Run()

    return test

def SolveAdjointProblem(kratos_parameters):
    test = AdjointRANSAnalysis(Kratos.Model(), kratos_parameters)
    test.Run()

    return test

@KratosUnittest.skipIfApplicationsNotAvailable("HDF5Application")
class AdjointKOmegaTwoElementsTest(KratosUnittest.TestCase):
    @staticmethod
    def _ReadParameters(parameters_file_name):
        parameters = ReadParameters(parameters_file_name)
        parameters["problem_data"]["parallel_type"].SetString("OpenMP")
        return parameters

    def testQSVMSKOmegaSteady(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointKOmegaTwoElementsTest._ReadParameters('./TwoElementsTest/kw_steady_test_parameters.json')
            step_size = 1e-10
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                AdjointKOmegaTwoElementsTest._AddHDF5PrimalOutputProcess,
                True)

            # solve adjoint
            adjoint_parameters = AdjointKOmegaTwoElementsTest._ReadParameters('./TwoElementsTest/kw_steady_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 3)

    def testQSVMSKOmegaBossak(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointKOmegaTwoElementsTest._ReadParameters('./TwoElementsTest/kw_bossak_test_parameters.json')
            step_size = 1.1111e-2
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                AdjointKOmegaTwoElementsTest._AddHDF5PrimalOutputProcess)

            # solve adjoint
            adjoint_parameters = AdjointKOmegaTwoElementsTest._ReadParameters('./TwoElementsTest/kw_bossak_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 1)

    @staticmethod
    def _AddHDF5PrimalOutputProcess(parameters):
        process_parameters = Kratos.Parameters("""
            {
                "python_module": "primal_hdf5_output_process",
                "kratos_module": "KratosMultiphysics.RANSApplication",
                "Parameters": {
                    "model_part_name": "MainModelPart"
                }
            }
        """)
        parameters["output_processes"].AddEmptyList("hdf5_output")
        parameters["output_processes"]["hdf5_output"].Append(process_parameters)

    @classmethod
    def tearDownClass(_):
        with KratosUnittest.WorkFolderScope('.', __file__):
            DeleteH5Files()
            DeleteTimeFiles(".")

if __name__ == '__main__':
    KratosUnittest.main()