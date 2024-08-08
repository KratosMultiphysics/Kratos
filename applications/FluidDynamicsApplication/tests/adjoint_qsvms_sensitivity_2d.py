import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ReadParameters
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import DeleteH5Files
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolvePrimalProblem
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolveAdjointProblem
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ComputeAdjointSensitivity
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceBodyFittedMomentShapeSensitivityAnalysis
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceBodyFittedDragFrequencyShapeSensitivityAnalysis

@KratosUnittest.skipIfApplicationsNotAvailable("HDF5Application")
class AdjointQSVMSSensitivity2D(KratosUnittest.TestCase):
    @staticmethod
    def _ReadParameters(parameters_file_name):
        parameters = ReadParameters(parameters_file_name)
        parameters["problem_data"]["parallel_type"].SetString("OpenMP")
        return parameters

    def testOneElement(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/one_element_test_parameters.json')
            step_size = 6.4e-5
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                AdjointQSVMSSensitivity2D.AddHDF5PrimalOutputProcess)

            # solve adjoint
            adjoint_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/one_element_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 4)

    def testOneElementMoment(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/one_element_test_parameters.json')
            step_size = 6.4e-5
            fd_sensitivities = FiniteDifferenceBodyFittedMomentShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [0.0, 0.0, 1.0],
                'MainModelPart.Structure',
                [1.0, 1.0, 0.0],
                SolvePrimalProblem,
                AdjointQSVMSSensitivity2D.AddHDF5PrimalOutputProcess)

            # solve adjoint
            adjoint_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/one_element_test_moment_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 4)

    def testCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/cylinder_test_parameters.json')
            step_size = 1e-7
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                AdjointQSVMSSensitivity2D.AddHDF5PrimalOutputProcess)

            # solve adjoint
            adjoint_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/cylinder_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 5)

    def testDragFrequencyRealCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/cylinder_test_parameters.json')
            step_size = 1e-9
            fd_sensitivities = FiniteDifferenceBodyFittedDragFrequencyShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                AdjointQSVMSSensitivity2D.AddHDF5PrimalOutputProcess, True, 0.05, 1)

            # solve adjoint
            adjoint_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/cylinder_test_adjoint_drag_frequency_real_coeff_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 5)

    def testDragFrequencyImagCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/cylinder_test_parameters.json')
            step_size = 1e-9
            fd_sensitivities = FiniteDifferenceBodyFittedDragFrequencyShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                AdjointQSVMSSensitivity2D.AddHDF5PrimalOutputProcess, False, 0.05, 1)

            # solve adjoint
            adjoint_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/cylinder_test_adjoint_drag_frequency_imag_coeff_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 5)

    def testSteadyCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/steady_cylinder_test_parameters.json')
            step_size = 1e-8
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                AdjointQSVMSSensitivity2D.AddHDF5PrimalOutputProcess)

            # solve adjoint
            adjoint_parameters = AdjointQSVMSSensitivity2D._ReadParameters('./AdjointQSVMSSensitivity2DTest/steady_cylinder_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 3)

    @staticmethod
    def AddHDF5PrimalOutputProcess(parameters):
        process_parameters = Kratos.Parameters(R'''
        {
            "kratos_module": "KratosMultiphysics.HDF5Application",
            "python_module": "single_mesh_temporal_output_process",
            "help": "",
            "process_name": "",
            "Parameters": {
                "model_part_name": "MainModelPart",
                "file_settings": {
                    "file_access_mode": "truncate"
                },
                "nodal_solution_step_data_settings": {
                    "list_of_variables": [
                        "VELOCITY",
                        "PRESSURE",
                        "REACTION"
                    ]
                },
                "nodal_data_value_settings": {
                    "list_of_variables": [
                        "RELAXED_ACCELERATION"
                    ]
                }
            }
        }
        ''')
        parameters["processes"]["auxiliar_process_list"].Append(process_parameters)

    def tearDown(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            DeleteH5Files()
            DeleteTimeFiles(".")

if __name__ == '__main__':
    KratosUnittest.main()
