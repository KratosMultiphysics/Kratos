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
class AdjointOneElementTest(KratosUnittest.TestCase):
    @staticmethod
    def _ReadParameters(parameters_file_name):
        parameters = ReadParameters(parameters_file_name)
        parameters["problem_data"]["parallel_type"].SetString("OpenMP")
        return parameters

    def testQSVMS(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointOneElementTest._ReadParameters('./OneElementTest/one_element_qsvms_test_parameters.json')
            step_size = 1e-6
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.Structure', SolvePrimalProblem)

            # solve adjoint
            adjoint_parameters = AdjointOneElementTest._ReadParameters('./OneElementTest/one_element_qsvms_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            print(adjoint_sensitivities)
            print(fd_sensitivities)
            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 5)

    def testQSVMSKEpsilon(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointOneElementTest._ReadParameters('./OneElementTest/one_element_ke_test_parameters.json')
            step_size = 1e-4
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.Structure', SolvePrimalProblem)

            # solve adjoint
            adjoint_parameters = AdjointOneElementTest._ReadParameters('./OneElementTest/one_element_ke_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            print(adjoint_sensitivities)
            print(fd_sensitivities)
            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 3)

    @classmethod
    def tearDownClass(_):
        with KratosUnittest.WorkFolderScope('.', __file__):
            DeleteH5Files()
            DeleteTimeFiles(".")

if __name__ == '__main__':
    KratosUnittest.main()