import os
import json

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess

class TestCalculateRomBasisOutputProcess(KratosUnittest.TestCase):

    def setUp(self):
        # Test data
        self.print_output = False
        self.work_folder = "calculate_rom_basis_output_process_test_files"

        # List containing the fake total times
        self.time_steps_list = [1.0,2.0,4.0]

        # Create a model part to perform the ROM basis calculation
        self.model = KratosMultiphysics.Model()
        model_part = self.model.CreateModelPart("MainModelPart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        n_nodes = 5
        for i in range(n_nodes):
            model_part.CreateNewNode(i+1,float(i),0.0,0.0)

    # @KratosUnittest.skipUnless(numpy_available, "numpy is required for RomApplication")
    def testCalculateRomBasisOutputProcess(self):
        # Create a CalculateROMBasisOutputProcess instance
        self.process_settings = KratosMultiphysics.Parameters("""{
            "model_part_name": "MainModelPart",
            "snapshots_control_type": "step",
            "snapshots_interval": 1.0,
            "nodal_unknowns": ["TEMPERATURE"],
            "rom_basis_output_format": "json",
            "rom_basis_output_name": "RomParameters",
            "svd_truncation_tolerance": 1.0e-6
        }""")

        # Run a "fake" simulation to calculate the ROM basis from its results
        self.__ExecuteTest(self.process_settings)

        # Check results
        self.__CheckResults()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            if not self.print_output:
                for file_name in os.listdir():
                    if not file_name.count("Results"):
                        kratos_utilities.DeleteFileIfExisting(file_name)

    def __ExecuteTest(self, process_settings):
        # Emulate a simulation to get the ROM basis output
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Create a calculate ROM basis process instance
            rom_basis_process = CalculateRomBasisOutputProcess(self.model, process_settings)

            # Emulate simulation to fill the database and calculate the ROM basis
            rom_basis_process.ExecuteInitialize()
            rom_basis_process.Check()
            rom_basis_process.ExecuteBeforeSolutionLoop()

            model_part = self.model.GetModelPart("MainModelPart")
            for time,step in zip(self.time_steps_list, range(1,len(self.time_steps_list)+1,1)):
                # Fake time advance
                model_part.CloneTimeStep(time)
                model_part.ProcessInfo[KratosMultiphysics.STEP] = step

                # Set nodal values
                if step % 2 == 0:
                    get_value = lambda node_id, time : node_id*time**2 + step
                else:
                    get_value = lambda node_id, time : -node_id/time**step
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, get_value(node.Id, time))

                # Call the process instance methods
                rom_basis_process.ExecuteInitializeSolutionStep()
                if rom_basis_process.IsOutputStep():
                    rom_basis_process.ExecuteBeforeOutputStep()
                    rom_basis_process.PrintOutput()
                    rom_basis_process.ExecuteAfterOutputStep()
                rom_basis_process.ExecuteFinalizeSolutionStep()

            rom_basis_process.ExecuteFinalize()

    def __CheckResults(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Load ROM basis output file
            output_filename = "{}.{}".format(self.process_settings["rom_basis_output_name"].GetString(), self.process_settings["rom_basis_output_format"].GetString())
            with open(output_filename) as f:
                output_data = json.load(f)

            # Load reference file
            reference_filename = "{}Results.{}".format(self.process_settings["rom_basis_output_name"].GetString(), self.process_settings["rom_basis_output_format"].GetString())
            with open(reference_filename) as f:
                reference_data = json.load(f)

            # Check files
            # Note that we need to do it manually as literal (deterministic) file comparison cannot be used because of the basis float values
            self.assertEqual(output_data["rom_settings"], reference_data["rom_settings"])
            for node_output, node_reference in zip(output_data["nodal_modes"],reference_data["nodal_modes"]):
                self.assertEqual(node_output, node_reference)
                for dof_basis_output, dof_basis_reference in zip(output_data["nodal_modes"][node_output], reference_data["nodal_modes"][node_reference]):
                    self.assertVectorAlmostEqual(dof_basis_output, dof_basis_reference)

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
