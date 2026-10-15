import json
import shutil
import numpy as np
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
from KratosMultiphysics.RomApplication.rom_manager import RomManager

if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
if kratos_utilities.CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication


class TestElementalLSPGRom(KratosUnittest.TestCase):
    """The elemental LSPG ROM must give the same solution as the global one (with the normal equations)."""

    def setUp(self):
        self.tests_folder = Path(__file__).parent
        self.relative_tolerance = 1.0e-9

    def _RunLSPG(self, source_folder, parameters_filename, assembling_strategy, variables):
        # The ROM data of an existing LSPG test is copied to a temporary folder at the same level (its paths are relative),
        # where the assembling strategy is set
        work_folder = f"{source_folder}_{assembling_strategy}_tmp"
        kratos_utilities.DeleteDirectoryIfExisting(str(self.tests_folder / work_folder))
        shutil.copytree(self.tests_folder / source_folder, self.tests_folder / work_folder)
        self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, str(self.tests_folder / work_folder))

        rom_parameters_path = self.tests_folder / work_folder / "rom_data" / "RomParameters.json"
        with open(rom_parameters_path, 'r') as rom_parameters_file:
            rom_parameters = json.load(rom_parameters_file)
        rom_parameters["assembling_strategy"] = assembling_strategy
        rom_parameters["rom_settings"].setdefault("rom_bns_settings", {})["solving_technique"] = "normal_equations"
        with open(rom_parameters_path, 'w') as rom_parameters_file:
            json.dump(rom_parameters, rom_parameters_file)

        with KratosUnittest.WorkFolderScope(work_folder, __file__):
            with open(parameters_filename, 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            dummy = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Save the results of each time step
            class DummyAnalysis(type(dummy)):
                def FinalizeSolutionStep(cls):
                    super().FinalizeSolutionStep()
                    cls.snapshots.append(rom_testing_utilities.GetNodalResults(cls._GetSolver().GetComputingModelPart(), variables))

            simulation = DummyAnalysis(model, parameters)
            simulation.snapshots = []
            simulation.Run()
            return np.array(simulation.snapshots).T

    def _CheckElementalAgainstGlobal(self, source_folder, parameters_filename, variables):
        global_solution = self._RunLSPG(source_folder, parameters_filename, "global", variables)
        elemental_solution = self._RunLSPG(source_folder, parameters_filename, "elemental", variables)
        self.assertGreater(np.linalg.norm(global_solution), 0.0)
        self.assertLess(np.linalg.norm(elemental_solution - global_solution)/np.linalg.norm(global_solution), self.relative_tolerance)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testStructuralStaticElementalLSPGRom2D(self):
        variables = [KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y]
        self._CheckElementalAgainstGlobal("structural_static_test_files/LSPGROM", "../ProjectParameters.json", variables)

    @KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
    def testFluidElementalLSPGRom2D(self):
        variables = [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.PRESSURE]
        self._CheckElementalAgainstGlobal("fluid_dynamics_test_files/LSPGROM", "../ProjectParameters.json", variables)

    @KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
    def testFluidElementalLSPGHRom2D(self):
        variables = [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.PRESSURE]
        self._CheckElementalAgainstGlobal("fluid_dynamics_test_files/LSPGHROM", "ProjectParametersHROM.json", variables)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testElementalLSPGRomManager(self):
        # ROM and HROM trained with the RomManager: the projected residuals of the HROM are computed with the elemental builder and solver
        work_folder = "test_rom_manager"
        rom_folder = "rom_data_elemental_lspg"
        self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, str(self.tests_folder / work_folder / rom_folder))

        def UpdateProjectParameters(parameters, mu=None):
            for load_process, value in zip(parameters["processes"]["loads_process_list"].values(), mu):
                load_process["Parameters"]["value"].SetString(f"({value})")
            return parameters

        general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "projection_strategy": "lspg",
            "assembling_strategy": "elemental",
            "ROM": {}
        }""")
        general_rom_manager_parameters["ROM"].AddString("rom_basis_output_folder", rom_folder)

        with KratosUnittest.WorkFolderScope(work_folder, __file__):
            rom_manager = RomManager(general_rom_manager_parameters=general_rom_manager_parameters, UpdateProjectParameters=UpdateProjectParameters)
            rom_manager.Fit([[0.1, 0.1, 0.1]])
            self.assertLess(rom_manager.ROMvsFOM["Fit"], 1.0e-8)
            self.assertLess(rom_manager.ROMvsHROM["Fit"], 1.0e-8)
            for file_name in Path(".").iterdir():
                if file_name.suffix == ".time":
                    kratos_utilities.DeleteFileIfExisting(str(file_name))

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
