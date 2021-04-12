import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles

from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import ReadParameters
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import DeleteH5Files
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolvePrimalProblem
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_sensitivity_utilities import SolveAdjointProblem


@UnitTest.skipIfApplicationsNotAvailable("HDF5Application")
class AdjointFluidTest(UnitTest.TestCase):

    def setUp(self):
        self.write_json_output = False

    @staticmethod
    def _ReadParameters(parameters_file_name):
        parameters = ReadParameters(parameters_file_name)
        if (Kratos.IsDistributedRun()):
            parameters["problem_data"]["parallel_type"].SetString("MPI")
        else:
            parameters["problem_data"]["parallel_type"].SetString("OpenMP")

        return parameters

    def testCylinder(self):
        with UnitTest.WorkFolderScope('.', __file__):
            primal_parameters = AdjointFluidTest._ReadParameters(
                './AdjointVMSSensitivity2DTest/cylinder_test_parameters.json')
            if (not self.write_json_output):
                AdjointFluidTest._AddJsonCheckProcess(primal_parameters, "./AdjointVMSSensitivity2DTest/cylinder_test_primal_results.json", 1e-9)
            else:
                AdjointFluidTest._AddJsonOutputProcess(primal_parameters, "./AdjointVMSSensitivity2DTest/cylinder_test_primal_results.json")
            AdjointFluidTest._AddHDF5PrimalOutputProcess(primal_parameters)
            SolvePrimalProblem(primal_parameters)

            adjoint_parameters = AdjointFluidTest._ReadParameters(
                './AdjointVMSSensitivity2DTest/cylinder_test_adjoint_parameters.json')
            AdjointFluidTest._AddAdjointProcesses(adjoint_parameters, "cylinder_test_adjoint")
            SolveAdjointProblem(adjoint_parameters)

    def testSlipCylinder(self):
        with UnitTest.WorkFolderScope('.', __file__):
            primal_parameters = AdjointFluidTest._ReadParameters(
                './AdjointVMSSensitivity2DTest/cylinder_slip_test_parameters.json')
            if (not self.write_json_output):
                AdjointFluidTest._AddJsonCheckProcess(primal_parameters, "./AdjointVMSSensitivity2DTest/cylinder_slip_test_primal_results.json", 1e-9)
            else:
                AdjointFluidTest._AddJsonOutputProcess(primal_parameters, "./AdjointVMSSensitivity2DTest/cylinder_slip_test_primal_results.json")
            AdjointFluidTest._AddHDF5PrimalSlipOutputProcess(primal_parameters)
            SolvePrimalProblem(primal_parameters)

            adjoint_parameters = AdjointFluidTest._ReadParameters(
                './AdjointVMSSensitivity2DTest/cylinder_slip_test_adjoint_parameters.json')
            AdjointFluidTest._AddAdjointProcesses(adjoint_parameters, "cylinder_slip_test_adjoint")
            SolveAdjointProblem(adjoint_parameters)

    @staticmethod
    def _AddJsonCheckProcess(kratos_parameters, input_file_name, tolerance):
        process_settings = Kratos.Parameters(R'''
        {
            "python_module": "from_json_check_result_process",
            "kratos_module": "KratosMultiphysics",
            "help": "",
            "process_name": "FromJsonCheckResultProcess",
            "Parameters": {
                "check_variables": [
                        "VELOCITY_X",
                        "VELOCITY_Y",
                        "PRESSURE",
                        "ACCELERATION_X",
                        "ACCELERATION_Y"
                ],
                "input_file_name": "<INPUT_FILE_NAME>",
                "model_part_name": "MainModelPart.Parts_Fluid",
                "tolerance": 1e-9,
                "relative_tolerance": 1e-12,
                "time_frequency": -2
            }
        }
        ''')
        process_settings["Parameters"]["input_file_name"].SetString(
            input_file_name)
        process_settings["Parameters"]["tolerance"].SetDouble(tolerance)
        kratos_parameters["processes"]["auxiliar_process_list"].Append(
            process_settings)

    @staticmethod
    def _AddJsonOutputProcess(kratos_parameters, output_file_name):
        process_settings = Kratos.Parameters(R'''
        {
            "python_module": "json_output_process",
            "kratos_module": "KratosMultiphysics",
            "help": "",
            "process_name": "JsonOutputProcess",
            "Parameters": {
                "output_variables": [
                        "VELOCITY_X",
                        "VELOCITY_Y",
                        "PRESSURE",
                        "ACCELERATION_X",
                        "ACCELERATION_Y"
                ],
                "output_file_name": "<OUTPUT_FILE_NAME>",
                "model_part_name": "MainModelPart.Parts_Fluid",
                "time_frequency": -2
            }
        }''')
        process_settings["Parameters"]["output_file_name"].SetString(
            output_file_name)
        kratos_parameters["processes"]["auxiliar_process_list"].Append(
            process_settings)

    @staticmethod
    def _AddAdjointProcesses(kratos_parameters, adjoint_prefix):
        AdjointFluidTest._AddPointOutputProcesses(
            kratos_parameters,
            [0.020957, 0.0055272, 0.0],
            adjoint_prefix + "_probe1.dat",
            ["ADJOINT_FLUID_VECTOR_1_X", "ADJOINT_FLUID_VECTOR_1_Y"])

        AdjointFluidTest._AddPointOutputProcesses(
            kratos_parameters,
            [0.014931,-0.0034173, 0.0],
            adjoint_prefix + "_probe2.dat",
            ["ADJOINT_FLUID_SCALAR_1"])

        AdjointFluidTest._AddPointOutputProcesses(
            kratos_parameters,
            [0.023303,-0.0037623, 0.0],
            adjoint_prefix + "_probe3.dat",
            ["SHAPE_SENSITIVITY_X", "SHAPE_SENSITIVITY_Y"])

        AdjointFluidTest._AddCompareTwoFilesCheckProcess(
            kratos_parameters,
            './AdjointVMSSensitivity2DTest/{:s}_probe1.dat'.format(adjoint_prefix),
            './AdjointVMSSensitivity2DTest/{:s}_probe1_ref.dat'.format(adjoint_prefix),
            1e-5)

        AdjointFluidTest._AddCompareTwoFilesCheckProcess(
            kratos_parameters,
            './AdjointVMSSensitivity2DTest/{:s}_probe2.dat'.format(adjoint_prefix),
            './AdjointVMSSensitivity2DTest/{:s}_probe2_ref.dat'.format(adjoint_prefix),
            1e-2)

        AdjointFluidTest._AddCompareTwoFilesCheckProcess(
            kratos_parameters,
            './AdjointVMSSensitivity2DTest/{:s}_probe3.dat'.format(adjoint_prefix),
            './AdjointVMSSensitivity2DTest/{:s}_probe3_ref.dat'.format(adjoint_prefix),
            1e-10)

    @staticmethod
    def _AddPointOutputProcesses(kratos_parameters, point_coordinates, output_file_name, output_variables_list):
        point_output_process = Kratos.Parameters(R'''
        {
            "kratos_module"   : "KratosMultiphysics",
            "python_module"   : "point_output_process",
            "help"            : "",
            "process_name"    : "PointOutputProcess",
            "Parameters" : {
                "position"         : [0.0, 0.0, 0.0],
                "model_part_name"  : "MainModelPart.Parts_Fluid",
                "output_file_settings": {
                    "output_path": "AdjointVMSSensitivity2DTest",
                    "file_name"  : "<OUTPUT_FILE_NAME>"
                },
                "output_variables" : []
            }
        }
        ''')

        position_vector = Kratos.Vector(3)
        position_vector[0] = point_coordinates[0]
        position_vector[1] = point_coordinates[1]
        position_vector[2] = point_coordinates[2]

        point_output_process["Parameters"]["position"].SetVector(
            position_vector)
        point_output_process["Parameters"]["output_file_settings"]["file_name"].SetString(
            output_file_name)
        point_output_process["Parameters"]["output_variables"].SetStringArray(
            output_variables_list)
        kratos_parameters["processes"]["auxiliar_process_list"].Append(
            point_output_process)


    @staticmethod
    def _AddCompareTwoFilesCheckProcess(kratos_parameters, output_file_name, reference_file_name, tolerance):
        check_process = Kratos.Parameters(R'''
        {
            "python_module"   : "compare_two_files_check_process",
            "kratos_module"   : "KratosMultiphysics",
            "help"            : "",
            "process_name"    : "CompareTwoFilesCheckProcess",
            "Parameters" :{
                "output_file_name"    : "<OUTPUT_FILE_NAME>",
                "reference_file_name" : "<REFERENCE_FILE_NAME>",
                "comparison_type"     : "dat_file",
                "remove_output_file"  : true,
                "tolerance"           : 1e-5
            }
        }
        ''')
        check_process["Parameters"]["output_file_name"].SetString(
            output_file_name)
        check_process["Parameters"]["reference_file_name"].SetString(
            reference_file_name)
        check_process["Parameters"]["tolerance"].SetDouble(tolerance)
        kratos_parameters["processes"]["auxiliar_process_list"].Append(
            check_process)

    @staticmethod
    def _AddHDF5PrimalOutputProcess(parameters):
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
                    "list_of_variables": ["VELOCITY", "PRESSURE", "ACCELERATION"]
                }
            }
        }
        ''')

        parameters["processes"]["auxiliar_process_list"].Append(process_parameters)

    @staticmethod
    def _AddHDF5PrimalSlipOutputProcess(parameters):
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
                        "BODY_FORCE",
                        "ACCELERATION"
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

        parameters["processes"]["auxiliar_process_list"].Append(process_parameters)

    @classmethod
    def tearDownClass(_):
        with UnitTest.WorkFolderScope('.', __file__):
            DeleteH5Files()
            DeleteTimeFiles(".")


if __name__ == '__main__':
    UnitTest.main()
