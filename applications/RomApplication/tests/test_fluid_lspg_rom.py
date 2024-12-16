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
class TestFluidLSPGRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testFluidLSPGRom2D(self):
        self.work_folder = "fluid_dynamics_test_files/LSPGROM/"
        parameters_filename = "../ProjectParameters.json"
        expected_output_filename = "ExpectedOutputLSPGROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            dummy = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Patch the RomAnalysis class to save the selected time steps results
            class DummyAnalysis(type(dummy)):
                def Initialize(cls):
                    super().Initialize()
                    cls.selected_time_step_solution_container = []
                
                def FinalizeSolutionStep(cls):
                    super().FinalizeSolutionStep()
                    variables_array = [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.PRESSURE]
                    array_of_results = rom_testing_utilities.GetNodalResults(cls._solver.GetComputingModelPart(), variables_array)
                    cls.selected_time_step_solution_container.append(array_of_results)
            
            self.simulation = DummyAnalysis(model, parameters)

            # Run test case
            self.simulation.Run()

            # Check results
            expected_output = np.load(expected_output_filename)
            n_values = len(self.simulation.selected_time_step_solution_container[0])
            n_snapshots = len(self.simulation.selected_time_step_solution_container)
            obtained_snapshot_matrix = np.zeros((n_values, n_snapshots))
            for i in range(n_snapshots):
                snapshot_i= np.array(self.simulation.selected_time_step_solution_container[i])
                obtained_snapshot_matrix[:,i] = snapshot_i.transpose()

            for i in range (n_snapshots):
                up = sum((expected_output[:,i] - obtained_snapshot_matrix[:,i])**2)
                down = sum((expected_output[:,i])**2)
                l2 = np.sqrt(up/down)
                self.assertLess(l2, self.relative_tolerance)
        
    def testFluidLSPGHRom2D(self):
        self.work_folder = "fluid_dynamics_test_files/LSPGHROM/"
        parameters_filename = "ProjectParametersHROM.json"
        expected_output_filename = "ExpectedOutputLSPGHROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            dummy = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Patch the RomAnalysis class to save the selected time steps results
            class DummyAnalysis(type(dummy)):
                def Initialize(cls):
                    super().Initialize()
                    cls.selected_time_step_solution_container = []
                
                def FinalizeSolutionStep(cls):
                    super().FinalizeSolutionStep()
                    variables_array = [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.PRESSURE]
                    array_of_results = rom_testing_utilities.GetNodalResults(cls._solver.GetComputingModelPart(), variables_array)
                    cls.selected_time_step_solution_container.append(array_of_results)
            
            self.simulation = DummyAnalysis(model, parameters)

            # Run test case
            self.simulation.Run()

            # Check results
            expected_output = np.load(expected_output_filename)
            n_values = len(self.simulation.selected_time_step_solution_container[0])
            n_snapshots = len(self.simulation.selected_time_step_solution_container)
            obtained_snapshot_matrix = np.zeros((n_values, n_snapshots))
            for i in range(n_snapshots):
                snapshot_i= np.array(self.simulation.selected_time_step_solution_container[i])
                obtained_snapshot_matrix[:,i] = snapshot_i.transpose()

            for i in range (n_snapshots):
                up = sum((expected_output[:,i] - obtained_snapshot_matrix[:,i])**2)
                down = sum((expected_output[:,i])**2)
                l2 = np.sqrt(up/down)
                self.assertLess(l2, self.relative_tolerance)

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
