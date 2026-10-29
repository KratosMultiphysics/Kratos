import os
import json
import numpy as np
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
from KratosMultiphysics.RomApplication.rom_manager import RomManager

if kratos_utilities.CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication


@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
class TestRbfEnhancedRom(KratosUnittest.TestCase):

    def setUp(self):
        # The FOM is the one of the ANN-enhanced ROM test. A temporary folder at the same level is used, since its paths are relative
        self.work_folder = "fluid_dynamics_test_files/RBF/"
        self.project_parameters_name = "../GALERKIN_HROM_ANN/ProjectParameters.json"
        os.makedirs(Path(__file__).parent / self.work_folder, exist_ok=True)
        self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, str(Path(__file__).parent / self.work_folder))

    def _GetRomManagerParameters(self, projection_strategy):
        general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM"],
            "rom_stages_to_test" : [],
            "type_of_decoder" : "rbf_enhanced",
            "assembling_strategy": "elemental",
            "ROM":{
                "svd_truncation_tolerance": 0,
                "model_part_name": "FluidModelPart",
                "nodal_unknowns": ["VELOCITY_X","VELOCITY_Y", "PRESSURE"],
                "rbf_enhanced_settings": {
                    "modes": [3, 10]
                }
            }
        }""")
        general_rom_manager_parameters.AddString("projection_strategy", projection_strategy)
        return general_rom_manager_parameters

    def _FitAndCheck(self, projection_strategy):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(project_parameters_name=self.project_parameters_name, general_rom_manager_parameters=self._GetRomManagerParameters(projection_strategy))
            rom_manager.Fit()
            self.assertLess(rom_manager.ROMvsFOM["Fit"], self.tolerance)

    def testFluidGalerkinRbfEnhancedRom2D(self):
        self.tolerance = 1.0e-6
        self._FitAndCheck("galerkin")

    def testFluidLSPGRbfEnhancedRom2D(self):
        self.tolerance = 1.0e-6
        self._FitAndCheck("lspg")

    def testFluidGalerkinRbfEnhancedRom2D_Export(self):
        # The exported model must run without the RomManager (no database) and give the results of the ROM run by the RomManager
        standalone_folder = Path(__file__).parent / "fluid_dynamics_test_files/RBF_STANDALONE"
        os.makedirs(standalone_folder, exist_ok=True)
        self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, str(standalone_folder))

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(project_parameters_name=self.project_parameters_name, general_rom_manager_parameters=self._GetRomManagerParameters("galerkin"))
            rom_manager.Fit()
            rom_manager.ExportRbfEnhancedRom(export_folder=standalone_folder / "rom_data")
            rom_manager_snapshots = rom_manager.data_base.get_snapshots_matrix_from_database([None], table_name='ROM')

        with KratosUnittest.WorkFolderScope(standalone_folder, __file__):
            for file_name in ["RightBasisMatrix.npy", "SingularValues.npy", "rbf_weights.npy", "rbf_centers.npy", "NodeIds.npy"]:
                self.assertTrue(Path("rom_data", file_name).exists(), msg=file_name)
            self.assertEqual(np.load("rom_data/RightBasisMatrix.npy").shape[1], 10)
            self.assertEqual(np.load("rom_data/rbf_centers.npy").shape[1], 3)
            self.assertEqual(np.load("rom_data/rbf_weights.npy").shape[1], 7)

            with open("rom_data/RomParameters.json") as f:
                rom_parameters = json.load(f)
            self.assertFalse(rom_parameters["rom_manager"])
            self.assertEqual(rom_parameters["projection_strategy"], "galerkin_rbf")
            self.assertIn(rom_parameters["rbf_enhanced_settings"]["kernel"], ["gaussian", "imq"])

            with open(self.project_parameters_name,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)
            simulation.Run()

            # The snapshots store the (alphabetically sorted) unknowns node by node
            variables = [KratosMultiphysics.PRESSURE, KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y]
            obtained = np.array(rom_testing_utilities.GetNodalResults(simulation._GetSolver().GetComputingModelPart(), variables))
            expected = rom_manager_snapshots[:, -1]
            self.assertLess(np.linalg.norm(obtained - expected) / np.linalg.norm(expected), 1.0e-10)

            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

    def _FitHromAndCheck(self, projection_strategy):
        # The HROM integrates the projected residuals of the RBF-enhanced ROM with a subset of the elements and conditions
        general_rom_manager_parameters = self._GetRomManagerParameters(projection_strategy)
        general_rom_manager_parameters["rom_stages_to_train"].SetStringArray(["ROM", "HROM"])
        general_rom_manager_parameters.AddValue("HROM", KratosMultiphysics.Parameters("""{"element_selection_svd_truncation_tolerance": 0}"""))
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # The HROM training writes the HROM model parts next to the input mesh (here in the FluidDynamicsApplication tests)
            with open(self.project_parameters_name, 'r') as parameter_file:
                input_filename = json.load(parameter_file)["solver_settings"]["model_import_settings"]["input_filename"]
            for suffix in ["HROM", "HROMVisualization"]:
                self.addCleanup(kratos_utilities.DeleteFileIfExisting, os.path.abspath(f"{input_filename}{suffix}.mdpa"))

            rom_manager = RomManager(project_parameters_name=self.project_parameters_name, general_rom_manager_parameters=general_rom_manager_parameters)
            rom_manager.Fit()
            self.assertLess(rom_manager.ROMvsFOM["Fit"], self.tolerance)
            self.assertLess(rom_manager.ROMvsHROM["Fit"], self.hrom_tolerance)

            self.assertGreater(np.load("rom_data/HROM_ElementIds.npy").size, 0)

    def testFluidGalerkinRbfEnhancedHRom2D_HRomModelPart(self):
        # The HROM run on the HROM model part (written during the HROM training) must give the solution of the HROM run on the full one
        general_rom_manager_parameters = self._GetRomManagerParameters("galerkin")
        general_rom_manager_parameters["rom_stages_to_train"].SetStringArray(["ROM", "HROM"])
        # The fluid solver requires the parent elements of the conditions in the HROM model part
        general_rom_manager_parameters.AddValue("HROM", KratosMultiphysics.Parameters("""{"element_selection_svd_truncation_tolerance": 0, "include_condition_parents": true}"""))

        results = {}
        def CustomizeSimulation(cls, global_model, parameters, mu=None):
            class CustomSimulation(cls):
                def Finalize(self):
                    model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
                    variables = [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.PRESSURE]
                    results[CustomizeSimulation.label] = {node.Id: [node.GetSolutionStepValue(variable) for variable in variables] for node in model_part.Nodes}
                    super().Finalize()
            return CustomSimulation(global_model, parameters)

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # The HROM training writes the HROM model parts next to the input mesh (here in the FluidDynamicsApplication tests)
            with open(self.project_parameters_name, 'r') as parameter_file:
                input_filename = json.load(parameter_file)["solver_settings"]["model_import_settings"]["input_filename"]
            for suffix in ["HROM", "HROMVisualization"]:
                self.addCleanup(kratos_utilities.DeleteFileIfExisting, os.path.abspath(f"{input_filename}{suffix}.mdpa"))

            rom_manager = RomManager(project_parameters_name=self.project_parameters_name, general_rom_manager_parameters=general_rom_manager_parameters, CustomizeSimulation=CustomizeSimulation)
            CustomizeSimulation.label = "Fit"
            rom_manager.Fit()
            for label, use_full_model_part in [("full model part", True), ("HROM model part", False)]:
                CustomizeSimulation.label = label
                rom_manager.RunHROM(use_full_model_part=use_full_model_part)

        # The HROM model part has only some of the nodes, on which both solutions must coincide
        full, hrom = results["full model part"], results["HROM model part"]
        self.assertLess(len(hrom), len(full))
        expected = np.array([full[node_id] for node_id in hrom])
        obtained = np.array([hrom[node_id] for node_id in hrom])
        self.assertLess(np.linalg.norm(obtained - expected) / np.linalg.norm(expected), 1.0e-9)

    def testFluidGalerkinRbfEnhancedHRom2D(self):
        # The HROM of the enhanced decoders is not exact (the tangent operator changes with the solution); the ANN-enhanced one gives ~1e-3 here
        self.tolerance = 1.0e-6
        # The round-off (e.g. the number of threads) can change the elements selected in the HROM training: ROMvsHROM from 4e-4 to 9e-4
        self.hrom_tolerance = 1.0e-2
        self._FitHromAndCheck("galerkin")

    def testFluidLSPGRbfEnhancedHRom2D(self):
        # The ANN-enhanced HROM gives ~7e-3 here
        self.tolerance = 1.0e-6
        self.hrom_tolerance = 1.0e-2
        self._FitHromAndCheck("lspg")

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
