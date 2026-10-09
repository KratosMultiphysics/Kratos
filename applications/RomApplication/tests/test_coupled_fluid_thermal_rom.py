import os
import types
import numpy as np
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
from KratosMultiphysics.RomApplication.rom_manager import RomManager


if kratos_utilities.CheckIfApplicationsAvailable("FluidDynamicsApplication", "ConvectionDiffusionApplication"):
    import KratosMultiphysics.FluidDynamicsApplication
    import KratosMultiphysics.ConvectionDiffusionApplication

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication", "ConvectionDiffusionApplication")
class TestCoupledFluidThermalRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testCoupledFluidThermalRom2D(self):
        self.work_folder = "coupled_fluid_thermal_test_files/ROM/"
        parameters_filename = "../ProjectParameters.json"
        expected_output_filename = "ExpectedOutputCoupledROM.npy"

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

                variables_array = [KratosMultiphysics.PRESSURE, KratosMultiphysics.TEMPERATURE, KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y]
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

    def testCoupledFluidThermalRomManager(self):
        # Train the ROM and HROM of the coupled problem with the RomManager, each coupled solver with its own settings
        self.work_folder = "coupled_fluid_thermal_test_files/RomManager/"
        with KratosUnittest.WorkFolderScope("coupled_fluid_thermal_test_files", __file__):
            os.makedirs("RomManager", exist_ok=True)
            # The ROM settings of the sub-solvers are set by the RomAnalysis, so they are removed from the FOM parameters
            with open("ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            for sub_solver_settings_name in ["fluid_solver_settings", "thermal_solver_settings"]:
                for key in ["rom_settings", "projection_strategy", "assembling_strategy"]:
                    parameters["solver_settings"][sub_solver_settings_name].RemoveValue(key)
            with open("RomManager/ProjectParameters.json",'w') as parameter_file:
                parameter_file.write(parameters.PrettyPrintJsonString())
        self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, str(Path(__file__).parent / self.work_folder))

        general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train": ["ROM", "HROM"],
            "projection_strategy": "galerkin",
            "assembling_strategy": "elemental",
            "ROM": {
                "model_part_name": "FluidModelPart",
                "nodal_unknowns": ["PRESSURE", "VELOCITY_X", "VELOCITY_Y"],
                "svd_truncation_tolerance": 1e-4
            },
            "HROM": {
                "element_selection_svd_truncation_tolerance": 1e-4
            },
            "coupled_solvers": [
                {"sub_solver_name": "fluid_solver"},
                {"sub_solver_name": "thermal_solver",
                 "ROM": {"model_part_name": "ThermalModelPart", "nodal_unknowns": ["TEMPERATURE"]},
                 "HROM": {"element_selection_svd_truncation_tolerance": 1e-5}}
            ]
        }""")

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(project_parameters_name="ProjectParameters.json", general_rom_manager_parameters=general_rom_manager_parameters)
            rom_manager.Fit()

            # A single (monolithic) basis is built with the nodal unknowns of all the coupled solvers
            with open("rom_data/RomParameters.json",'r') as rom_parameters_file:
                rom_parameters = KratosMultiphysics.Parameters(rom_parameters_file.read())
            self.assertEqual(rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray(), ["PRESSURE", "TEMPERATURE", "VELOCITY_X", "VELOCITY_Y"])
            self.assertEqual(rom_parameters["coupled_solvers"].GetStringArray(), ["fluid_solver", "thermal_solver"])

            # One column of HROM weights per coupled solver
            hrom_element_weights = np.load("rom_data/HROM_ElementWeights.npy")
            self.assertEqual(hrom_element_weights.shape[1], 2)

            self.assertLess(rom_manager.ROMvsFOM["Fit"], 1.0e-2)
            self.assertLess(rom_manager.ROMvsHROM["Fit"], 1.0e-3)

    def testCoupledSolversSettingsValidation(self):
        self.work_folder = "coupled_fluid_thermal_test_files"
        project_parameters_name = str(Path(__file__).parent / "coupled_fluid_thermal_test_files/ProjectParameters.json")
        with KratosUnittest.WorkFolderScope("coupled_fluid_thermal_test_files", __file__):
            self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, str(Path(__file__).parent / "coupled_fluid_thermal_test_files/rom_data"))

            # Settings shared by all the coupled solvers cannot be set per coupled solver
            general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
                "coupled_solvers": [{"sub_solver_name": "fluid_solver", "ROM": {"rom_basis_output_folder": "other_folder"}}]
            }""")
            with self.assertRaisesRegex(Exception, "'rom_basis_output_folder' cannot be set in the 'ROM' settings of the coupled solver 'fluid_solver'"):
                RomManager(project_parameters_name=project_parameters_name, general_rom_manager_parameters=general_rom_manager_parameters)

            # The monolithic basis uses the smallest tolerance and the nodal unknowns of all the coupled solvers
            general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
                "ROM": {"nodal_unknowns": ["PRESSURE"], "svd_truncation_tolerance": 1e-3},
                "coupled_solvers": [
                    {"sub_solver_name": "fluid_solver", "ROM": {"svd_truncation_tolerance": 1e-5}},
                    {"sub_solver_name": "thermal_solver", "ROM": {"nodal_unknowns": ["TEMPERATURE"]}}
                ]
            }""")
            rom_manager = RomManager(project_parameters_name=project_parameters_name, general_rom_manager_parameters=general_rom_manager_parameters)
            rom_settings = rom_manager.general_rom_manager_parameters["ROM"]
            self.assertEqual(rom_settings["svd_truncation_tolerance"].GetDouble(), 1e-5)
            self.assertEqual(rom_settings["nodal_unknowns"].GetStringArray(), ["PRESSURE", "TEMPERATURE"])

            # Settings not given per coupled solver are taken from the general ones
            thermal_solver = rom_manager._GetCoupledSolvers()[1]
            self.assertEqual(thermal_solver["ROM"]["svd_truncation_tolerance"].GetDouble(), 1e-3)
            self.assertEqual(thermal_solver["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble(),
                             rom_manager.general_rom_manager_parameters["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble())

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
