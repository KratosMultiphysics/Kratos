import os
import shutil
import types
import numpy as np
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication as KratosROM
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
from KratosMultiphysics.RomApplication.rom_manager import RomManager


if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication", "ConvectionDiffusionApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
    import KratosMultiphysics.ConvectionDiffusionApplication

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "ConvectionDiffusionApplication")
class TestCoupledThermoMechanicalRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-10

    def testCoupledThermoMechanicalRom2D(self):
        self.work_folder = "coupled_thermo_mechanical_test_files/ROM/"
        self._RunAndCheck("ExpectedOutputCoupledROM.npy")

    def testCoupledThermoMechanicalHRom2D(self):
        self.work_folder = "coupled_thermo_mechanical_test_files/HROM/"
        self._RunAndCheck("ExpectedOutputCoupledHROM.npy")
        self._CheckHromWeights()

    def testCoupledThermoMechanicalAnnEnhancedRom2D(self):
        # ANN-enhanced ROM run without the RomManager: each coupled solver has its own basis and network
        self.work_folder = "coupled_thermo_mechanical_test_files/ROM_ANN/"
        self._RunAndCheck("ExpectedOutputCoupledROM_ANN.npy")

    def testCoupledThermoMechanicalAnnEnhancedHRom2D(self):
        # ANN-enhanced HROM run without the RomManager: one set of HROM weights per coupled solver
        self.work_folder = "coupled_thermo_mechanical_test_files/HROM_ANN/"
        self._RunAndCheck("ExpectedOutputCoupledHROM_ANN.npy")
        self._CheckHromWeights()

    def testCoupledThermoMechanicalSegregatedBasesRomManager(self):
        # Train the ROM and HROM with the RomManager with a basis per coupled solver. A shared basis is near-singular in the
        # displacement rows (displacement ~1e-7, temperature ~1e2), so the reduced structural system would be set by the round-off
        self.work_folder = "coupled_thermo_mechanical_test_files/RomManagerSegregatedBases/"
        self._WriteRomManagerProjectParameters()

        general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train": ["ROM", "HROM"],
            "projection_strategy": "galerkin",
            "assembling_strategy": "elemental",
            "coupled_solvers_segregated_bases": true,
            "ROM": {
                "model_part_name": "Structure",
                "nodal_unknowns": ["DISPLACEMENT_X", "DISPLACEMENT_Y"],
                "svd_truncation_tolerance": 1e-6
            },
            "HROM": {
                "element_selection_svd_truncation_tolerance": 1e-6
            },
            "coupled_solvers": [
                {"sub_solver_name": "structural_solver"},
                {"sub_solver_name": "thermal_solver",
                 "ROM": {"model_part_name": "ThermalModelPart", "nodal_unknowns": ["TEMPERATURE"]}}
            ]
        }""")

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(project_parameters_name="ProjectParameters.json", general_rom_manager_parameters=general_rom_manager_parameters)
            rom_manager.Fit()

            # The bases of the coupled solvers are stored side by side from the first column, zero-padded up to the largest one.
            # Each solver uses its own modes. The rows store DISPLACEMENT_X, DISPLACEMENT_Y and TEMPERATURE node by node
            with open("rom_data/RomParameters.json",'r') as rom_parameters_file:
                rom_parameters = KratosMultiphysics.Parameters(rom_parameters_file.read())
            basis = np.load("rom_data/RightBasisMatrix.npy").reshape(-1, 3, rom_parameters["rom_settings"]["number_of_rom_dofs"].GetInt())
            for sub_solver_name, components in [("structural_solver", [0, 1]), ("thermal_solver", [2])]:
                number_of_rom_dofs = rom_parameters["coupled_solvers_number_of_rom_dofs"][sub_solver_name].GetInt()
                self.assertTrue(np.all(np.any(basis[:, components, :number_of_rom_dofs] != 0.0, axis=(0, 1))))
                self.assertTrue(np.all(basis[:, components, number_of_rom_dofs:] == 0.0))

            self.assertLess(rom_manager.ROMvsFOM["Fit"], 1.0e-5)
            self.assertLess(rom_manager.ROMvsHROM["Fit"], 1.0e-7)

    def _WriteRomManagerProjectParameters(self):
        # The mesh and materials are copied to the work folder, as the HROM training writes the HROM mesh next to the input mesh
        parameters = self._GetProjectParameters()
        work_folder = Path(__file__).parent / self.work_folder
        os.makedirs(work_folder, exist_ok=True)
        self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, str(work_folder))
        for sub_solver_settings_name in ["structural_solver_settings", "thermal_solver_settings"]:
            sub_solver_settings = parameters["solver_settings"][sub_solver_settings_name]
            for import_settings_name, file_name_key, extension in [("model_import_settings", "input_filename", ".mdpa"), ("material_import_settings", "materials_filename", "")]:
                file_name = sub_solver_settings[import_settings_name][file_name_key]
                source_file = Path(file_name.GetString() + extension)
                shutil.copy(source_file, work_folder / source_file.name)
                file_name.SetString(source_file.stem if extension else source_file.name)
        parameters.AddValue("output_processes", KratosMultiphysics.Parameters("""{}"""))
        with open(work_folder / "ProjectParameters.json",'w') as parameter_file:
            parameter_file.write(parameters.PrettyPrintJsonString())

    def _GetProjectParameters(self):
        # The problem (parameters, mdpa and materials) is the thermo-mechanical test of the ConvectionDiffusionApplication.
        # Its paths are relative to the tests folder of that application, so they are made absolute
        convection_diffusion_tests_folder = Path(__file__).resolve().parents[2] / "ConvectionDiffusionApplication" / "tests"
        with open(convection_diffusion_tests_folder / "thermo_mechanical_tests/thermo_mechanical/coupled_problem_test_parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters.AddString("analysis_stage", "KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis")
        for sub_solver_settings_name in ["structural_solver_settings", "thermal_solver_settings"]:
            sub_solver_settings = parameters["solver_settings"][sub_solver_settings_name]
            for file_name in [sub_solver_settings["model_import_settings"]["input_filename"], sub_solver_settings["material_import_settings"]["materials_filename"]]:
                file_name.SetString(str(convection_diffusion_tests_folder / file_name.GetString()))
        parameters["processes"].RemoveValue("json_check_process") # It checks against the FOM results
        return parameters

    def _RunSimulation(self):
        # Set up simulation
        model = KratosMultiphysics.Model()
        self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, self._GetProjectParameters())

        # Patch the RomAnalysis class to save the selected time steps results
        self.simulation.selected_time_step_solution_container = []

        def FinalizeSolutionStep(cls):
            super(type(self.simulation), cls).FinalizeSolutionStep()

            variables_array = [KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.TEMPERATURE]
            array_of_results = rom_testing_utilities.GetNodalResults(cls._solver.GetComputingModelPart(), variables_array)
            cls.selected_time_step_solution_container.append(array_of_results)

        self.simulation.FinalizeSolutionStep  = types.MethodType(FinalizeSolutionStep, self.simulation)

        # Run test case
        self.simulation.Run()

        return np.array(self.simulation.selected_time_step_solution_container).T

    def _RunAndCheck(self, expected_output_filename):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            obtained_snapshot_matrix = self._RunSimulation()

            # Check results
            expected_output = np.load(expected_output_filename)
            self.assertEqual(obtained_snapshot_matrix.shape, expected_output.shape)
            # Each field is checked on its own, as the temperature (~1e2) would hide the displacement (~1e-7) in a common norm.
            # The rows store DISPLACEMENT_X, DISPLACEMENT_Y and TEMPERATURE node by node
            number_of_steps = expected_output.shape[1]
            expected_fields = expected_output.reshape(-1, 3, number_of_steps)
            obtained_fields = obtained_snapshot_matrix.reshape(-1, 3, number_of_steps)
            for components in [[0, 1], [2]]:
                for i in range(number_of_steps):
                    up = np.sum((expected_fields[:,components,i] - obtained_fields[:,components,i])**2)
                    down = np.sum((expected_fields[:,components,i])**2)
                    l2 = np.sqrt(up/down)
                    self.assertLess(l2, self.relative_tolerance)

    def _CheckHromWeights(self):
        # The coupled solvers share the geometries, where the element and condition data is stored.
        # Each solver must only find HROM_WEIGHT in the HROM entities, not in the rest of the mesh
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            hrom_element_ids = set(np.load("rom_data/HROM_ElementIds.npy"))
            hrom_condition_ids = set(np.load("rom_data/HROM_ConditionIds.npy"))
        for solver in [self.simulation._GetSolver().thermal_solver, self.simulation._GetSolver().structural_solver]:
            model_part = solver.GetComputingModelPart()
            self.assertSetEqual({element.Id for element in model_part.Elements if element.Has(KratosROM.HROM_WEIGHT)}, hrom_element_ids)
            self.assertSetEqual({condition.Id for condition in model_part.Conditions if condition.Has(KratosROM.HROM_WEIGHT)}, hrom_condition_ids)

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
