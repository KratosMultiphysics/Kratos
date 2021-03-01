import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ReadParameters
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import DeleteH5Files
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolvePrimalProblem
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolveAdjointProblem
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ComputeAdjointSensitivity
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis


@KratosUnittest.skipUnless(CheckIfApplicationsAvailable("HDF5Application"), "Missing HDF5Application")
class AdjointVMSSensitivity2D(KratosUnittest.TestCase):
    @staticmethod
    def _ReadParameters(parameters_file_name):
        parameters = ReadParameters(parameters_file_name)
        parameters["problem_data"]["parallel_type"].SetString("OpenMP")
        return parameters

    def testOneElement(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/one_element_test_parameters.json')
            step_size = 0.00000001
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.Structure', SolvePrimalProblem)

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/one_element_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 4)

    def testCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/cylinder_test_parameters.json')
            step_size = 0.00000001
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.NoSlip2D_Cylinder', SolvePrimalProblem)

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/cylinder_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 5)

    def testSteadyCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/steady_cylinder_test_parameters.json')
            step_size = 1e-10
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.NoSlip2D_Cylinder', SolvePrimalProblem)

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/steady_cylinder_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 3)

    @classmethod
    def tearDownClass(_):
        with KratosUnittest.WorkFolderScope('.', __file__):
            DeleteH5Files()
            DeleteTimeFiles(".")


if __name__ == '__main__':
    KratosUnittest.main()
