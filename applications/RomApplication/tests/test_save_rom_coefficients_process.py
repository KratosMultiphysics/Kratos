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

    # def testSaveRomCoefficientsProcess(self):
    #     self.work_folder = "structural_static_test_files/ROM/"
    #     parameters_filename = "../ProjectParameters.json"

    #     with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
    #         # Set up simulation
    #         with open(parameters_filename,'r') as parameter_file:
    #             parameters = KratosMultiphysics.Parameters(parameter_file.read())

    #         self.process_settings = KratosMultiphysics.Parameters("""{
    #             "Parameters": {
    #                 "model_part_name": "Structure",
    #                 "rom_coefficients_output_folder": "../../save_rom_coefficients_process_test_files/save_rom_coefficients_process_test",
    #                 "rom_coefficients_output_name": "ObtainedOutputSaveRomCoefficientsProcess",
    #                 "snapshots_control_type": "step",
    #                 "snapshots_interval": 1.0,
    #                 "snapshot_solution_type": "incremental"
    #             },
    #             "kratos_module": "KratosMultiphysics.RomApplication",
    #             "process_name": "SaveROMCoefficientsProcess",
    #             "python_module": "save_rom_coefficients_process"
    #         }""")

    #         parameters["output_processes"]["rom_output"].Append(self.process_settings)
    #         model = KratosMultiphysics.Model()
    #         self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

    #         # Run test case
    #         self.simulation.Run()

    #         # Check results
    #         self.__CheckResults()

    def testSaveRomCoefficientsProcessWithLineSearch(self):
        self.work_folder = "compressible_potential_test_files/"
        parameters_filename = "ProjectParameters_LineSearch.json"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.process_settings = KratosMultiphysics.Parameters("""{
                "Parameters": {
                    "model_part_name": "FluidModelPart",
                    "rom_coefficients_output_folder": "../save_rom_coefficients_process_test_files/save_rom_coefficients_process_test_line_search",
                    "rom_coefficients_output_name": "ObtainedOutputSaveRomCoefficientsProcessLineSearch",
                    "snapshots_control_type": "step",
                    "snapshots_interval": 1.0,
                    "snapshot_solution_type": "incremental"
                },
                "kratos_module": "KratosMultiphysics.RomApplication",
                "process_name": "SaveROMCoefficientsProcess",
                "python_module": "save_rom_coefficients_process"
            }""")

            parameters["output_processes"]["rom_output"].Append(self.process_settings)
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            exit()

            # Check results
            self.__CheckResults()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Construct the full path to the output directory
            output_directory_path = os.path.join(os.getcwd(), self.output_folder)
            # Delete the output directory and its contents
            kratos_utilities.DeleteDirectoryIfExisting(output_directory_path)

    def __CheckResults(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Load ROM coefficients output file
            self.output_folder = self.process_settings["Parameters"]["rom_coefficients_output_folder"].GetString()
            rom_coeffs_obtained_output_name = os.path.join(self.output_folder, "ObtainedOutputSaveRomCoefficientsProcess.npy")
            rom_coeffs_obtained_output = np.load(rom_coeffs_obtained_output_name)

            # Load reference file
            rom_coeffs_expected_output_name = "../../save_rom_coefficients_process_test_files/ExpectedOutputSaveRomCoefficientsProcess.npy"
            rom_coeffs_expected_output = np.load(rom_coeffs_expected_output_name)

            self.assertMatrixAlmostEqual(KratosMultiphysics.Matrix(rom_coeffs_obtained_output), KratosMultiphysics.Matrix(rom_coeffs_expected_output))


##########################################################################################

if __name__ == '__main__':
    # KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
