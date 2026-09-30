import copy
import json
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.RomApplication import rom_parameters


class TestRomParameters(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = Path("rom_parameters_test_output")
        # Version 1 files currently committed with the tests
        tests_folder = Path(__file__).parent
        self.version_1_files = sorted(tests_folder.glob("**/rom_data/RomParameters*.json"))
        with (tests_folder / "coupled_fluid_thermal_test_files/ROM/rom_data/RomParameters.json").open('r') as f:
            self.coupled_version_1 = json.load(f)

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            kratos_utilities.DeleteDirectoryIfExisting(str(self.work_folder))

    def test_committed_files_round_trip(self):
        self.assertGreater(len(self.version_1_files), 0)
        for file_path in self.version_1_files:
            with self.subTest(file=str(file_path)):
                with file_path.open('r') as f:
                    version_1 = json.load(f)
                version_2 = rom_parameters.UpgradeToVersion2(version_1)
                self.assertTrue(rom_parameters.IsVersion2(version_2))
                self.assertEqual(len(version_2["components"]), 1)
                self.assertEqual(rom_parameters.NormalizeVersion1(rom_parameters.DowngradeToVersion1(version_2)), rom_parameters.NormalizeVersion1(version_1))

    def test_coupled_solvers(self):
        # Content written by the RomManager for coupled solvers sharing a single basis
        version_1 = copy.deepcopy(self.coupled_version_1)
        version_1["coupled_solvers"] = ["fluid_solver", "thermal_solver"]
        version_2 = rom_parameters.UpgradeToVersion2(version_1)
        self.assertEqual(version_2["coupling"], {"basis_layout" : "monolithic", "solvers" : ["fluid_solver", "thermal_solver"]})
        self.assertNotIn("legacy", version_2)
        self.assertEqual(rom_parameters.DowngradeToVersion1(version_2), version_1)

    def test_ann_enhanced(self):
        for strategy in ["galerkin", "lspg"]:
            with self.subTest(strategy=strategy):
                version_1 = copy.deepcopy(self.coupled_version_1)
                version_1["projection_strategy"] = f"{strategy}_ann"
                version_2 = rom_parameters.UpgradeToVersion2(version_1)
                component = version_2["components"][0]
                self.assertEqual(component["projection"]["strategy"], strategy)
                self.assertEqual(component["manifold"]["type"], "global")
                self.assertEqual(component["manifold"]["decoder"]["type"], "ann_enhanced")
                self.assertEqual(rom_parameters.DowngradeToVersion1(version_2), version_1)

    def test_numpy_files_are_referenced(self):
        version_2 = rom_parameters.UpgradeToVersion2(self.coupled_version_1)
        manifold = version_2["components"][0]["manifold"]
        self.assertEqual(manifold["type"], "global")
        self.assertEqual(manifold["decoder"]["type"], "linear")
        basis = manifold["decoder"]["basis"]
        self.assertEqual(basis["file"], "RightBasisMatrix.npy")
        self.assertEqual(basis["node_ids"], "NodeIds.npy")
        self.assertEqual(basis["number_of_modes"], self.coupled_version_1["rom_settings"]["number_of_rom_dofs"])
        self.assertEqual(version_2["hyper_reduction"]["element_weights"], "HROM_ElementWeights.npy")

    def test_read_and_write(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            # A version 2 file is read by the version 1 readers
            version_2 = rom_parameters.UpgradeToVersion2(self.coupled_version_1)
            self.work_folder.mkdir(parents=True, exist_ok=True)
            with rom_parameters.GetRomParametersFilePath(self.work_folder, "RomParametersV2").open('w') as f:
                json.dump(version_2, f)
            self.assertEqual(rom_parameters.ReadRomParametersAsVersion1(self.work_folder, "RomParametersV2"), self.coupled_version_1)
            self.assertEqual(rom_parameters.LoadRomParameters(self.work_folder, "RomParametersV2"), version_2)

            # Files are written in the version 1 layout for the moment
            rom_parameters.WriteRomParameters(self.work_folder, "RomParametersV1", version_2)
            with rom_parameters.GetRomParametersFilePath(self.work_folder, "RomParametersV1").open('r') as f:
                self.assertEqual(json.load(f), self.coupled_version_1)

    def test_unsupported_downgrade(self):
        version_2 = rom_parameters.UpgradeToVersion2(self.coupled_version_1)
        version_2["components"].append(copy.deepcopy(version_2["components"][0]))
        with self.assertRaisesRegex(Exception, "2 components"):
            rom_parameters.DowngradeToVersion1(version_2)

        version_2 = rom_parameters.UpgradeToVersion2(self.coupled_version_1)
        version_2["components"][0]["manifold"]["decoder"]["basis"]["file"] = "OtherBasis.npy"
        with self.assertRaisesRegex(Exception, "RightBasisMatrix.npy"):
            rom_parameters.DowngradeToVersion1(version_2)

        # Local manifolds only exist in the version 2 layout
        version_2 = rom_parameters.UpgradeToVersion2(self.coupled_version_1)
        decoder = version_2["components"][0]["manifold"]["decoder"]
        version_2["components"][0]["manifold"] = {
            "type" : "local",
            "selector" : {"type" : "nearest_centroid", "centroids" : "centroids.npy"},
            "clusters" : [{"decoder" : decoder}, {"decoder" : decoder}]
        }
        with self.assertRaisesRegex(Exception, "'local' manifold"):
            rom_parameters.DowngradeToVersion1(version_2)

        # So do the nonlinear decoders other than the ANN-enhanced one
        version_2 = rom_parameters.UpgradeToVersion2(self.coupled_version_1)
        version_2["components"][0]["manifold"]["decoder"]["type"] = "rbf"
        with self.assertRaisesRegex(Exception, "'rbf' decoder"):
            rom_parameters.DowngradeToVersion1(version_2)


if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
