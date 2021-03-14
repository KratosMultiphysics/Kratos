
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ReadParameters
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import DeleteH5Files
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolvePrimalProblem
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolveAdjointProblem
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ComputeAdjointSensitivity
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceVelocityPressureNormSquareShapeSensitivityAnalysis


@KratosUnittest.skipIfApplicationsNotAvailable("HDF5Application")
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
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                lambda params : AdjointVMSSensitivity2D.AddHDF5PrimalOutputProcess(params, False))

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/one_element_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 4)

    def testTwoElementsSlipSteady(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/two_elements_steady_test_parameters.json')
            step_size = 1e-9
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                lambda params : AdjointVMSSensitivity2D.AddHDF5PrimalOutputProcess(params, True))

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/two_elements_steady_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 1)

    def testTwoElementsSlipBossak(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/two_elements_bossak_test_parameters.json')
            step_size = 1e-9
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                lambda params : AdjointVMSSensitivity2D.AddHDF5PrimalOutputProcess(params, False))

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/two_elements_bossak_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 4)

    def testCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/cylinder_test_parameters.json')
            step_size = 0.00000001
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                lambda params : AdjointVMSSensitivity2D.AddHDF5PrimalOutputProcess(params, False))

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
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                lambda params : AdjointVMSSensitivity2D.AddHDF5PrimalOutputProcess(params, True))

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/steady_cylinder_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 3)

    def testSlipNormCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/cylinder_slip_test_parameters.json')
            step_size = 1e-10
            fd_sensitivities = FiniteDifferenceVelocityPressureNormSquareShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters,
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                lambda params : AdjointVMSSensitivity2D.AddHDF5PrimalSlipOutputProcess(params, False))

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/cylinder_slip_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 5)

    def testSlipSteadyNormCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1968]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/steady_cylinder_slip_test_parameters.json')
            step_size = 1e-9
            fd_sensitivities = FiniteDifferenceVelocityPressureNormSquareShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters,
                'MainModelPart.NoSlip2D_Cylinder',
                SolvePrimalProblem,
                lambda params : AdjointVMSSensitivity2D.AddHDF5PrimalSlipOutputProcess(params, True))

            # solve adjoint
            adjoint_parameters = AdjointVMSSensitivity2D._ReadParameters('./AdjointVMSSensitivity2DTest/steady_cylinder_slip_test_norm_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 6)

    @staticmethod
    def AddHDF5PrimalOutputProcess(parameters, is_steady):
        process_parameters = Kratos.Parameters(R'''
        {
            "kratos_module" : "KratosMultiphysics.HDF5Application",
            "python_module" : "single_mesh_primal_output_process",
            "help"          : "",
            "process_name"  : "",
            "Parameters" : {
                "model_part_name" : "MainModelPart",
                "file_settings" : {
                    "file_access_mode" : "truncate"
                },
                "nodal_solution_step_data_settings" : {
                    "list_of_variables": ["VELOCITY", "PRESSURE"]
                }
            }
        }
        ''')
        if (not is_steady):
            process_parameters["Parameters"]["nodal_solution_step_data_settings"]["list_of_variables"].Append("ACCELERATION")

        parameters["processes"]["auxiliar_process_list"].Append(process_parameters)

    @staticmethod
    def AddHDF5PrimalSlipOutputProcess(parameters, is_steady):
        process_parameters = Kratos.Parameters(R'''
        {
            "kratos_module": "KratosMultiphysics.HDF5Application",
            "python_module": "single_mesh_primal_output_process",
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
                        "NORMAL",
                        "BODY_FORCE"
                    ]
                },
                "nodal_data_value_settings": {
                    "list_of_variables": [
                        "Y_WALL"
                    ]
                },
                "nodal_flag_value_settings": {
                    "list_of_variables": [
                        "SLIP"
                    ]
                },
                "condition_flag_value_settings": {
                    "list_of_variables": [
                        "SLIP"
                    ]
                }
            }
        }
        ''')
        if (not is_steady):
            process_parameters["Parameters"]["nodal_solution_step_data_settings"]["list_of_variables"].Append("ACCELERATION")

        parameters["processes"]["auxiliar_process_list"].Append(process_parameters)

    def tearDown(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            DeleteH5Files()
            DeleteTimeFiles(".")


if __name__ == '__main__':
    KratosUnittest.main()
