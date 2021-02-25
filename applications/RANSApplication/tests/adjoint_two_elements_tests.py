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
class AdjointTwoElementsTest(KratosUnittest.TestCase):
    @staticmethod
    def _ReadParameters(parameters_file_name):
        parameters = ReadParameters(parameters_file_name)
        parameters["problem_data"]["parallel_type"].SetString("OpenMP")
        return parameters

    def testQSVMSSteadySlip(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/qsvms_test_parameters.json')
            step_size = 1e-9
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.Structure', SolvePrimalProblem)

            # solve adjoint
            adjoint_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/qsvms_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 3)

    # def testQSVMSSteadySlip(self):
    #     with KratosUnittest.WorkFolderScope('.', __file__):
    #         node_ids = [4]

    #         # calculate sensitivity by finite difference
    #         primal_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/one_element_qsvms_slip_steady_test_parameters.json')

    #         # total_steps = [1e-8, 5e-9, 4e-9, 3e-9, 2e-9, 1e-9, 8e-10, 5e-10, 2e-10, 1e-10, 5e-11]
    #         # sensitivities = []
    #         # for step_size in total_steps:
    #         #     # step_size = 1e-9
    #         #     fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
    #         #         node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.Structure', SolvePrimalProblem, True)
    #         #     sensitivities.append(fd_sensitivities)

    #         # for i, steps in enumerate(total_steps):
    #         #     print(steps, sensitivities[i])

    #         step_size = 2e-10
    #         fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
    #             node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0], 'MainModelPart.Structure', SolvePrimalProblem, True)

    #         # # solve unperturbed
    #         # SolvePrimalProblem(primal_parameters)

    #         # solve adjoint
    #         adjoint_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/one_element_qsvms_slip_steady_test_adjoint_parameters.json')
    #         adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

    #         # fd_sensitivities = Kratos.Matrix(1, 2)
    #         # fd_sensitivities[0, 0] = 0.0422357
    #         # fd_sensitivities[0, 1] = 0.033163
    #         print(adjoint_sensitivities)
    #         print(fd_sensitivities)
    #         self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 5)

    def testQSVMSKEpsilonSteady(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/ke_steady_test_parameters.json')
            step_size = 1e-4
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                AdjointTwoElementsTest._AddHDF5PrimalOutputProcess,
                True)

            # solve adjoint
            adjoint_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/ke_steady_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            print(adjoint_sensitivities)
            print(fd_sensitivities)
            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 2)

    def testQSVMSKEpsilonBossak(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/ke_bossak_test_parameters.json')
            step_size = 1e-4
            fd_sensitivities = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, [1.0, 0.0, 0.0],
                'MainModelPart.Structure',
                SolvePrimalProblem,
                AdjointTwoElementsTest._AddHDF5PrimalOutputProcess)

            # solve adjoint
            adjoint_parameters = AdjointTwoElementsTest._ReadParameters('./TwoElementsTest/ke_bossak_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            print(adjoint_sensitivities)
            print(fd_sensitivities)
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