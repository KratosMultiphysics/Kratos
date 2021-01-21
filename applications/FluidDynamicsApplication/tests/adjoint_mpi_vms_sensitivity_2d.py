import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis


@KratosUnittest.skipUnless(CheckIfApplicationsAvailable("HDF5Application"), "Missing HDF5Application")
class AdjointMPIVMSSensitivity(KratosUnittest.TestCase):

    @staticmethod
    def _remove_h5_files(model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                DeleteFileIfExisting(name)

    @staticmethod
    def _AddJsonProcess(parameters_string, output_file_name):
        json_check_results_process = R'''
            {
                "python_module": "from_json_check_result_process",
                "kratos_module": "KratosMultiphysics",
                "help": "",
                "process_name": "FromJsonCheckResultProcess",
                "Parameters": {
                    "check_variables": [
                        "VELOCITY",
                        "PRESSURE",
                        "ACCELERATION",
                        "NORMAL"
                    ],
                    "input_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_JSON_CHECK_FILE_NAME>",
                    "model_part_name": "MainModelPart.Parts_Fluid",
                    "tolerance": 1e-9,
                    "relative_tolerance": 1e-12,
                    "time_frequency": 0.04
                }
            }
        '''

        json_output_process = R'''
            {
                "python_module": "json_output_process",
                "kratos_module": "KratosMultiphysics",
                "help": "",
                "process_name": "JsonOutputProcess",
                "Parameters": {
                    "output_variables": [
                        "VELOCITY",
                        "PRESSURE",
                        "ACCELERATION",
                        "NORMAL"
                    ],
                    "output_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_JSON_CHECK_FILE_NAME>",
                    "model_part_name": "MainModelPart.Parts_Fluid",
                        "time_frequency": 0.04
                }
            }
        '''

        parameters_string = parameters_string.replace("<OUTPUT_PROCESSES>", json_check_results_process)
        # parameters_string = parameters_string.replace("<OUTPUT_PROCESSES>", json_output_process)

        parameters_string = parameters_string.replace("<OUTPUT_JSON_CHECK_FILE_NAME>", output_file_name)

        return parameters_string

    @staticmethod
    def _AddAdjointOutputProcesses(parameters_string, output_file_name_prefix):
        adjoint_output_processes = R'''
            ,
            {
                "kratos_module": "KratosMultiphysics",
                "python_module": "point_output_process",
                "help": "",
                "process_name": "PointOutputProcess",
                "Parameters": {
                    "position": [
                        0.020957,
                        0.0055272,
                        0.0
                    ],
                    "model_part_name": "MainModelPart.Parts_Fluid",
                    "output_file_settings": {
                        "file_name": "<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe1.dat",
                        "output_path": "AdjointVMSSensitivity2DTest",
                        "write_buffer_size": 1
                    },
                    "output_variables": [
                        "ADJOINT_FLUID_VECTOR_1_X",
                        "ADJOINT_FLUID_VECTOR_1_Y"
                    ]
                }
            },
            {
                "kratos_module": "KratosMultiphysics",
                "python_module": "point_output_process",
                "help": "",
                "process_name": "PointOutputProcess",
                "Parameters": {
                    "position": [
                        0.014931,
                        -0.0034173,
                        0.0
                    ],
                    "model_part_name": "MainModelPart.Parts_Fluid",
                    "output_file_settings": {
                        "file_name": "<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe2.dat",
                        "output_path": "AdjointVMSSensitivity2DTest",
                        "write_buffer_size": 1
                    },
                    "output_variables": [
                        "ADJOINT_FLUID_SCALAR_1"
                    ]
                }
            },
            {
                "kratos_module": "KratosMultiphysics",
                "python_module": "point_output_process",
                "help": "",
                "process_name": "PointOutputProcess",
                "Parameters": {
                    "position": [
                        0.023303,
                        -0.0037623,
                        0.0
                    ],
                    "model_part_name": "MainModelPart.Parts_Fluid",
                    "output_file_settings": {
                        "file_name": "<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe3.dat",
                        "output_path": "AdjointVMSSensitivity2DTest",
                        "write_buffer_size": 1
                    },
                    "output_variables": [
                        "SHAPE_SENSITIVITY_X",
                        "SHAPE_SENSITIVITY_Y"
                    ]
                }
            },
            {
                "python_module": "compare_two_files_check_process",
                "kratos_module": "KratosMultiphysics",
                "help": "",
                "process_name": "CompareTwoFilesCheckProcess",
                "Parameters": {
                    "output_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe1.dat",
                    "reference_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe1_ref.dat",
                    "comparison_type": "dat_file",
                    "tolerance": 0.01,
                    "remove_output_file": true
                }
            },
            {
                "python_module": "compare_two_files_check_process",
                "kratos_module": "KratosMultiphysics",
                "help": "",
                "process_name": "CompareTwoFilesCheckProcess",
                "Parameters": {
                    "output_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe2.dat",
                    "reference_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe2_ref.dat",
                    "comparison_type": "dat_file",
                    "tolerance": 0.01,
                    "remove_output_file": true
                }
            },
            {
                "python_module": "compare_two_files_check_process",
                "kratos_module": "KratosMultiphysics",
                "help": "",
                "process_name": "CompareTwoFilesCheckProcess",
                "Parameters": {
                    "output_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe3.dat",
                    "reference_file_name": "AdjointVMSSensitivity2DTest/<OUTPUT_FILE_NAME_PREFIX_ADJOINT>_probe3_ref.dat",
                    "comparison_type": "dat_file",
                    "tolerance": 1e-7,
                    "remove_output_file": true
                }
            }
        '''

        parameters_string = parameters_string.replace("<OUTPUT_PROCESSES>", adjoint_output_processes)
        parameters_string = parameters_string.replace("<OUTPUT_FILE_NAME_PREFIX_ADJOINT>", output_file_name_prefix)

        return parameters_string

    def setUp(self):
        if IsDistributedRun():
            self.parallel_type = "MPI"
        else:
            self.parallel_type = "OpenMP"

    def testCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            # solve fluid
            model = Kratos.Model()
            with open('AdjointVMSSensitivity2DTest/cylinder_test_parameters.json', 'r') as parameter_file:
                lines = parameter_file.read()

            lines = AdjointMPIVMSSensitivity._AddJsonProcess(lines, "cylinder_test_primal_results.json")
            lines = lines.replace("<PARALLEL_TYPE>", self.parallel_type)

            project_parameters = Kratos.Parameters(lines)
            primal_simulation = FluidDynamicsAnalysis(model,project_parameters)
            primal_simulation.Run()
            Kratos.DataCommunicator.GetDefault().Barrier()

            # solve adjoint
            with open('AdjointVMSSensitivity2DTest/cylinder_test_adjoint_parameters.json', 'r') as parameter_file:
                lines = parameter_file.read()

            lines = AdjointMPIVMSSensitivity._AddAdjointOutputProcesses(lines, "cylinder_test_adjoint")
            lines = lines.replace("<PARALLEL_TYPE>", self.parallel_type)

            project_parameters = Kratos.Parameters(lines)
            adjoint_model = Kratos.Model()
            adjoint_simulation = AdjointFluidAnalysis(adjoint_model,project_parameters)
            adjoint_simulation.Run()

            Kratos.DataCommunicator.GetDefault().Barrier()

    def testSlipCylinder(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            # solve fluid
            model = Kratos.Model()
            with open('AdjointVMSSensitivity2DTest/cylinder_slip_test_parameters.json', 'r') as parameter_file:
                lines = parameter_file.read()

            lines = AdjointMPIVMSSensitivity._AddJsonProcess(lines, "cylinder_slip_test_primal_results.json")
            lines = lines.replace("<PARALLEL_TYPE>", self.parallel_type)

            project_parameters = Kratos.Parameters(lines)
            primal_simulation = FluidDynamicsAnalysis(model,project_parameters)
            primal_simulation.Run()
            Kratos.DataCommunicator.GetDefault().Barrier()

            # solve adjoint
            with open('AdjointVMSSensitivity2DTest/cylinder_slip_test_adjoint_parameters.json', 'r') as parameter_file:
                lines = parameter_file.read()

            lines = AdjointMPIVMSSensitivity._AddAdjointOutputProcesses(lines, "cylinder_slip_test_adjoint")
            lines = lines.replace("<PARALLEL_TYPE>", self.parallel_type)

            project_parameters = Kratos.Parameters(lines)
            adjoint_model = Kratos.Model()
            adjoint_simulation = AdjointFluidAnalysis(adjoint_model,project_parameters)
            adjoint_simulation.Run()
            Kratos.DataCommunicator.GetDefault().Barrier()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope('.', __file__):
            DeleteTimeFiles('./AdjointVMSSensitivity2DTest/')
            DeleteDirectoryIfExisting('./AdjointVMSSensitivity2DTest/cylinder_test_partitioned')
            AdjointMPIVMSSensitivity._remove_h5_files("MainModelPart")

if __name__ == '__main__':
    KratosUnittest.main()
