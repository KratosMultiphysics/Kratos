import os
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
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

    def testRbfEnhancedHromNotAvailable(self):
        general_rom_manager_parameters = self._GetRomManagerParameters("galerkin")
        general_rom_manager_parameters["rom_stages_to_train"].SetStringArray(["HROM"])
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(project_parameters_name=self.project_parameters_name, general_rom_manager_parameters=general_rom_manager_parameters)
            with self.assertRaisesRegex(Exception, "HROM is not available yet for rbf_enhanced decoders"):
                rom_manager.Fit()

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
