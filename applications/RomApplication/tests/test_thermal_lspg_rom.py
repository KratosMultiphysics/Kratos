import os
import types
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    import KratosMultiphysics.ConvectionDiffusionApplication

@KratosUnittest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
class TestThermalLSPGRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testConvDiffStationaryLSPGRom2D(self):
        self.work_folder = "thermal_static_test_files/LSPGROM/"
        parameters_filename = "../ProjectParameters.json"
        expected_output_filename = "ExpectedOutputLSPGROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            # Check results
            expected_output = np.load(expected_output_filename)
            obtained_output = rom_testing_utilities.GetScalarNodalResults(self.simulation._GetSolver().GetComputingModelPart(), KratosMultiphysics.TEMPERATURE)
            nodal_area = rom_testing_utilities.GetNodalAreaVector(self.simulation._GetSolver().GetComputingModelPart())

            l2 = np.sqrt((sum(nodal_area*((1 - obtained_output/expected_output )**2)))/(sum(nodal_area)))*100
            self.assertLess(l2, self.relative_tolerance)

    def testConvDiffDynamicRom2D(self):
        self.work_folder = "thermal_dynamic_test_files/LSPGROM/"
        parameters_filename = "../ProjectParameters.json"
        expected_output_filename = "ExpectedOutputLSPGROM.npy"

        time_snapshots = [500,1200,2500,3000,3600]

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up and run simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Patch the RomAnalysis class to save the selected time steps results
            def Initialize(cls):
                super(type(self.simulation), cls).Initialize()
                cls.selected_time_step_solution_container = []

            def FinalizeSolutionStep(cls):
                super(type(self.simulation), cls).FinalizeSolutionStep()

                array_of_temperatures = []
                if int(cls._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]) in time_snapshots:
                    for node in cls._solver.GetComputingModelPart().Nodes:
                        array_of_temperatures.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
                    cls.selected_time_step_solution_container.append(array_of_temperatures)

            self.simulation.Initialize  = types.MethodType(Initialize, self.simulation)
            self.simulation.FinalizeSolutionStep  = types.MethodType(FinalizeSolutionStep, self.simulation)

            # Run test case
            self.simulation.Run()

            # Check results
            expected_output = np.load(expected_output_filename)
            n_nodes = len(self.simulation.selected_time_step_solution_container[0])
            n_snapshots = len(self.simulation.selected_time_step_solution_container)
            obtained_snapshot_matrix = np.zeros((n_nodes, n_snapshots))
            for i in range(n_snapshots):
                snapshot_i= np.array(self.simulation.selected_time_step_solution_container[i])
                obtained_snapshot_matrix[:,i] = snapshot_i.transpose()
            nodal_area = rom_testing_utilities.GetNodalAreaVector(self.simulation._GetSolver().GetComputingModelPart())

            for i in range (n_snapshots):
                l2 = np.sqrt((sum(nodal_area*((1 - obtained_snapshot_matrix[:,i]/expected_output[:,i] )**2)))/(sum(nodal_area)))*100
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
