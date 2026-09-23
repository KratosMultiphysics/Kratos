import os
import types
import numpy as np
from pathlib import Path

import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_manager import RomManager
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
class TestFluidRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testFluidRom2D(self):
        self.work_folder = "fluid_dynamics_test_files/ROM/"
        parameters_filename = "../ProjectParameters.json"
        expected_output_filename = "ExpectedOutputROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
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

                variables_array = [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.PRESSURE]
                array_of_results = rom_testing_utilities.GetNodalResults(cls._solver.GetComputingModelPart(), variables_array)
                cls.selected_time_step_solution_container.append(array_of_results)

            self.simulation.Initialize  = types.MethodType(Initialize, self.simulation)
            self.simulation.FinalizeSolutionStep  = types.MethodType(FinalizeSolutionStep, self.simulation)

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


    def testFluidGalerkinRom2D_ANN(self):
        self.work_folder = "fluid_dynamics_test_files/GALERKIN_HROM_ANN/"
        expected_output_filename = "ExpectedOutputGalerkinHROM_ANN.npy"
        parameters_filename = "ProjectParametersROM.json"

        general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM","HROM"],
            "rom_stages_to_test" : [],
            "projection_strategy": "galerkin",
            "type_of_decoder" : "ann_enhanced",
            "assembling_strategy": "elemental",
            "save_gid_output": false,
            "save_vtk_output": false,
            "ROM":{
                "svd_truncation_tolerance": 0,
                "model_part_name": "FluidModelPart",
                "nodal_unknowns": ["VELOCITY_X","VELOCITY_Y", "PRESSURE"],
                "ann_enhanced_settings": {
                    "training": {
                        "retrain_if_exists": false
                    },
                    "modes": [3, 10],
                    "layers_size": [50, 50],
                    "batch_size": 2,
                    "epochs": 200,
                    "NN_gradient_regularisation_weight": 1.0,
                    "lr_strategy": {
                        "scheduler": "sgdr",
                        "base_lr": 0.001,
                        "additional_params": [0.0001, 10, 400]
                    },
                    "online": {
                        "model_number": 0
                    }
                }
            },
            "HROM":{
                "element_selection_svd_truncation_tolerance": 0
            }
        }""")

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(
                project_parameters_name=parameters_filename,
                general_rom_manager_parameters=general_rom_manager_parameters
            )
            rom_manager.RunHROM(use_full_model_part=True)
            expected_output = np.load(expected_output_filename)
            n_values = expected_output.shape[0]
            n_snapshots = expected_output.shape[1]

            obtained_snapshot_matrix = np.zeros((n_values, n_snapshots))
            snapshots_folder = Path("rom_data/Snapshots")

            for i in range(n_snapshots):
                step = i + 1  # Assuming step output starts at 1
                snapshot_path = snapshots_folder / f"solution_{step}.npy"

                self.assertTrue(snapshot_path.exists(), f"Snapshot file {snapshot_path} was not generated.")

                # Load the step data and place it in the corresponding column
                snapshot_data = np.load(snapshot_path)
                obtained_snapshot_matrix[:, i] = snapshot_data.flatten()

            # Check L2 error column by column (step by step)
            for i in range(n_snapshots):
                up = np.sum((expected_output[:, i] - obtained_snapshot_matrix[:, i])**2)
                down = np.sum((expected_output[:, i])**2)

                # Avoid division by zero if the expected solution is strictly 0
                if down > 1e-12:
                    l2 = np.sqrt(up / down)
                    self.assertLess(l2, self.relative_tolerance, msg=f"Tolerance failed at step {i+1}")
                else:
                    self.assertLess(np.sqrt(up), self.relative_tolerance, msg=f"Absolute tolerance failed at step {i+1}")

            self.tearDown()


    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Cleaning
            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
        kratos_utilities.DeleteDirectoryIfExisting("rom_data/Snapshots")

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
