import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis

from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ReadParameters
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import DeleteH5Files
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ComputeAdjointSensitivity
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import FiniteDifferenceDomainIntegratedShapeSensitivityAnalysis

def SolvePrimalProblem(kratos_parameters):
    test = RANSAnalysis(Kratos.Model(), kratos_parameters)
    test.Run()

    return test

def SolveAdjointProblem(kratos_parameters):
    test = AdjointRANSAnalysis(Kratos.Model(), kratos_parameters)
    test.Run()

    return test

@KratosUnittest.skipIfApplicationsNotAvailable("HDF5Application")
class AdjointCircularConvectionTwoElementsTest(KratosUnittest.TestCase):
    @staticmethod
    def _ReadParameters(parameters_file_name):
        parameters = ReadParameters(parameters_file_name)
        parameters["problem_data"]["parallel_type"].SetString("OpenMP")
        return parameters

    def testCircularConvectionSteady(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointCircularConvectionTwoElementsTest._ReadParameters('./TwoElementsTest/circular_convection_steady_test_parameters.json')
            step_size = 1e-4
            fd_sensitivities = FiniteDifferenceDomainIntegratedShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, 'MainModelPart',
                SolvePrimalProblem,
                AdjointCircularConvectionTwoElementsTest._AddHDF5PrimalOutputProcess,
                True)

            # solve adjoint
            adjoint_parameters = AdjointCircularConvectionTwoElementsTest._ReadParameters('./TwoElementsTest/circular_convection_steady_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 2)

    def testCircularConvectionBossak(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            node_ids = [1]

            # calculate sensitivity by finite difference
            primal_parameters = AdjointCircularConvectionTwoElementsTest._ReadParameters('./TwoElementsTest/circular_convection_bossak_test_parameters.json')
            step_size = 1e-8
            fd_sensitivities = FiniteDifferenceDomainIntegratedShapeSensitivityAnalysis.ComputeSensitivity(
                node_ids, step_size, primal_parameters, 'MainModelPart',
                SolvePrimalProblem,
                AdjointCircularConvectionTwoElementsTest._AddHDF5PrimalOutputProcess)

            # solve adjoint
            adjoint_parameters = AdjointCircularConvectionTwoElementsTest._ReadParameters('./TwoElementsTest/circular_convection_bossak_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 2)

    @staticmethod
    def _AddHDF5PrimalOutputProcess(parameters):
        process_parameters = Kratos.Parameters("""
            {
                "python_module": "single_mesh_temporal_output_process",
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "file_settings": {
                        "file_name": "<model_part_name>-<time>.h5",
                        "file_access_mode": "truncate"
                    },
                    "nodal_solution_step_data_settings": {
                        "list_of_variables": ["MESH_VELOCITY", "NORMAL", "RANS_AUXILIARY_VARIABLE_1", "VELOCITY_POTENTIAL", "VELOCITY_POTENTIAL_RATE"]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": ["RELAXED_ACCELERATION"]
                    },
                    "condition_data_value_settings": {
                        "list_of_variables": []
                    }
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