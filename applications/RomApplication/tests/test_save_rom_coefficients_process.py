import os
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestStructuralRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testSaveRomCoefficientsProcess(self):
        self.work_folder = "structural_static_test_files/ROM/"
        parameters_filename = "../ProjectParameters.json"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.process_settings = KratosMultiphysics.Parameters("""{
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "rom_coefficients_output_folder": "rom_data_test",
                    "rom_coefficients_output_name": "ObtainedOutputSaveRomCoefficientsProcess",
                    "snapshots_control_type": "step",
                    "snapshots_interval": 1.0,
                    "snapshot_solution_type": "incremental"
                },
                "kratos_module": "KratosMultiphysics.RomApplication",
                "process_name": "SaveROMCoefficientsProcess",
                "python_module": "save_rom_coefficients_process"
            }""")

            # parameters["output_processes"]["rom_output"].Append(self.process_settings)
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            # Check results
            # self.__CheckResults()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            for path in os.listdir():
                full_path = os.path.join(os.getcwd(), path)
                if not "Expected" in path:
                    if os.path.isfile(full_path):
                        kratos_utilities.DeleteFileIfExisting(full_path)
                    elif os.path.isdir(full_path):
                        kratos_utilities.DeleteDirectoryIfExisting(full_path)

    def __CheckResults(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Load ROM coefficients output file
            output_folder = self.process_settings["rom_coefficients_output_folder"].GetString()
            rom_coeffs_obtained_output_name = os.path.join(output_folder, "ObtainedOutputSaveRomCoefficientsProcess.npy")
            rom_coeffs_obtained_output = np.load(rom_coeffs_obtained_output_name)

            # Load reference file
            rom_coeffs_expected_output_name = "ExpectedOutputSaveRomCoefficientsProcess.npy"
            rom_coeffs_expected_output = np.load(rom_coeffs_expected_output_name)

            for i in range(rom_coeffs_obtained_output.shape[0]):
                for j in range(rom_coeffs_obtained_output.shape[1]):
                    self.assertAlmostEqual(rom_coeffs_obtained_output[i,j], rom_coeffs_expected_output[i,j])

            #JSON part tested in the JSON test already


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
