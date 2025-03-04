import os
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
if kratos_utilities.CheckIfApplicationsAvailable("CompressiblePotentialFlowApplication"):
    import KratosMultiphysics.CompressiblePotentialFlowApplication

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
@KratosUnittest.skipIfApplicationsNotAvailable("CompressiblePotentialFlowApplication")
class TestSaveRomCoefficientsProcess(KratosUnittest.TestCase):

    def testSaveRomCoefficientsProcess(self):
        self.work_folder = "save_rom_coefficients_process_test_files/structural_newton_raphson/"
        parameters_filename = "ProjectParameters.json"

        with KratosUnittest.WorkFolderScope(self.work_folder_newton_raphson, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            process_settings = KratosMultiphysics.Parameters("""{
                "Parameters": {
                    "model_part_name": "Structure",
                    "rom_coefficients_output_folder": "",
                    "rom_coefficients_output_name": "ObtainedOutputSaveRomCoefficientsNewtonRaphson",
                    "snapshots_control_type": "step",
                    "snapshots_interval": 1.0,
                    "snapshot_solution_type": "incremental"
                },
                "kratos_module": "KratosMultiphysics.RomApplication",
                "process_name": "SaveROMCoefficientsProcess",
                "python_module": "save_rom_coefficients_process"
            }""")

            parameters["output_processes"]["rom_output"].Append(process_settings)
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            # Determine the output path and the expected output path
            output_folder = process_settings["Parameters"]["rom_coefficients_output_folder"].GetString()
            output_filename = process_settings["Parameters"]["rom_coefficients_output_name"].GetString()+".npy"
            self.output_path = os.path.join(output_folder, output_filename)
            expected_path_newton_raphson = "ExpectedOutputSaveRomCoefficientsNewtonRaphson.npy"

            # Check results
            self.__CheckResults(self.work_folder, self.output_path, expected_path_newton_raphson)

    def testSaveRomCoefficientsProcessWithLineSearch(self):
        self.work_folder = "save_rom_coefficients_process_test_files/compressible_potential_flow_line_search/"
        parameters_filename = "ProjectParameters.json"

        with KratosUnittest.WorkFolderScope(self.work_folder_line_search, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            process_settings = KratosMultiphysics.Parameters("""{
                "Parameters": {
                    "model_part_name": "FluidModelPart",
                    "rom_coefficients_output_folder": "",
                    "rom_coefficients_output_name": "ObtainedOutputSaveRomCoefficientsLineSearch",
                    "snapshots_control_type": "step",
                    "snapshots_interval": 1.0,
                    "snapshot_solution_type": "total"
                },
                "kratos_module": "KratosMultiphysics.RomApplication",
                "process_name": "SaveROMCoefficientsProcess",
                "python_module": "save_rom_coefficients_process"
            }""")

            parameters["output_processes"]["rom_output"].Append(process_settings)
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            # Determine the output path and the expected output path
            output_folder = process_settings["Parameters"]["rom_coefficients_output_folder"].GetString()
            output_filename = process_settings["Parameters"]["rom_coefficients_output_name"].GetString()+".npy"
            self.output_path = os.path.join(output_folder, output_filename)
            expected_path_line_search = "ExpectedOutputSaveRomCoefficientsLineSearch.npy"

            # Check results
            self.__CheckResults(self.work_folder, self.output_path, expected_path_line_search)

    def tearDown(self):
        # Construct the full path to the output file and delete it
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            output_file_path = os.path.join(os.getcwd(), self.output_path)
            kratos_utilities.DeleteFileIfExisting(output_file_path)

    def __CheckResults(self, work_folder, obtained_output_path, expected_output_path):
        with KratosUnittest.WorkFolderScope(work_folder, __file__):
            # Load ROM coefficients output file
            rom_coeffs_obtained_output = np.load(obtained_output_path)
            # Load reference file
            rom_coeffs_expected_output = np.load(expected_output_path)

            self.assertMatrixAlmostEqual(KratosMultiphysics.Matrix(rom_coeffs_obtained_output), KratosMultiphysics.Matrix(rom_coeffs_expected_output))


##########################################################################################

if __name__ == '__main__':
    # KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
