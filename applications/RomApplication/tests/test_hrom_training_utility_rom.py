import os
import types
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
class TestHromTrainingUtilityRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testHromTrainingUtilityRom(self):
        self.work_folder = "hrom_training_utility_test_files/"
        parameters_filename = "ProjectParameters.json"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            debug = True

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Cleaning
            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
