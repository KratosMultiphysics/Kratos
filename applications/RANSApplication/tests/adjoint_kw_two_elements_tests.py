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
                lambda x: AdjointKOmegaTwoElementsTest._AddHDF5PrimalOutputProcess(x, False),
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
                lambda x: AdjointKOmegaTwoElementsTest._AddHDF5PrimalOutputProcess(x, True))

            # solve adjoint
            adjoint_parameters = AdjointKOmegaTwoElementsTest._ReadParameters('./TwoElementsTest/kw_bossak_test_adjoint_parameters.json')
            adjoint_sensitivities = ComputeAdjointSensitivity(node_ids, adjoint_parameters, SolveAdjointProblem)

            self.assertMatrixAlmostEqual(adjoint_sensitivities, fd_sensitivities, 1)

    @staticmethod
    def _AddHDF5PrimalOutputProcess(parameters, use_time = True):
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
                        "list_of_variables": [
                            "ACCELERATION",
                            "PRESSURE",
                            "MESH_VELOCITY",
                            "VELOCITY",
                            "TURBULENT_KINETIC_ENERGY",
                            "RANS_AUXILIARY_VARIABLE_2",
                            "DISPLACEMENT",
                            "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE",
                            "TURBULENT_KINETIC_ENERGY_RATE",
                            "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2",
                            "RANS_AUXILIARY_VARIABLE_1"
                        ]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": ["RELAXED_ACCELERATION"]
                    },
                    "nodal_flag_value_settings": {
                        "list_of_variables": ["SLIP"]
                    },
                    "condition_data_value_settings": {
                        "list_of_variables": [
                            "GAUSS_RANS_Y_PLUS",
                            "DISTANCE",
                            "RANS_IS_WALL_FUNCTION_ACTIVE"
                        ]
                    },
                    "condition_flag_value_settings": {
                        "list_of_variables": ["SLIP"]
                    }
                }
            }
        """)
        if not use_time:
            process_parameters["Parameters"]["file_settings"]["file_name"].SetString("<model_part_name>")
        parameters["output_processes"].AddEmptyList("hdf5_output")
        parameters["output_processes"]["hdf5_output"].Append(process_parameters)

    @classmethod
    def tearDownClass(_):
        with KratosUnittest.WorkFolderScope('.', __file__):
            DeleteH5Files()
            DeleteTimeFiles(".")

if __name__ == '__main__':
    KratosUnittest.main()