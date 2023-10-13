import os
import numpy as np
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities

if kratos_utilities.CheckIfApplicationsAvailable("CompressiblePotentialFlowApplication"):
    import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosPotentialFlow


@KratosUnittest.skipIfApplicationsNotAvailable("CompressiblePotentialFlowApplication")
class TestCompressiblePotentialRom(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "compressible_potential_test_files"
        self.relative_tolerance = 1.0e-12

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Cleaning
            files_to_remove = [f for f in os.listdir() if f.endswith(".time")]
            for file_name in files_to_remove:
                kratos_utilities.DeleteFileIfExisting(file_name)

    def testCompressiblePotentialRom(self):
        parameters_filename = "ProjectParameters.json"
        expected_output_filename = "ExpectedOutput.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename, 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            # Check results
            obtained_output = rom_testing_utilities.GetScalarNodalResults(
                self.simulation._GetSolver().GetComputingModelPart(), KratosPotentialFlow.VELOCITY_POTENTIAL)

            expected_output = np.load(expected_output_filename)
            nodal_area = rom_testing_utilities.GetNodalAreaVector(self.simulation._GetSolver().GetComputingModelPart())

            l2 = np.sqrt((sum(nodal_area*((1 - obtained_output/expected_output)**2)))/(sum(nodal_area)))*100
            self.assertLess(l2, self.relative_tolerance)


if __name__ == '__main__':
    KratosUnittest.main()
