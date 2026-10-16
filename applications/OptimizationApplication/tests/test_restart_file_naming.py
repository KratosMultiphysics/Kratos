from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.optimization_restart import OptimizationProblemRestartOutput
from KratosMultiphysics.OptimizationApplication.utilities.restart_file_naming import ParseStepFromFileName, ListCheckpoints, ResolveRestartLoadStep

class TestRestartFileNaming(kratos_unittest.TestCase):
    restart_files_path = Path("restart_file_naming_test")
    restart_file_name = "restart_<step>.pkl"

    def setUp(self) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteDirectoryIfExisting(str(self.restart_files_path))
            self.restart_files_path.mkdir()
            for name in ["restart_2.pkl", "restart_10.pkl", "restart_7.pkl", "restart_5_old.pkl", "restart_.pkl", "restart_3.pkl.tmp", "restart_12.pkl.tmp", "restart_\u00b2.pkl"]:
                (self.restart_files_path / name).touch()
            # a directory with a checkpoint-like name is not a checkpoint
            (self.restart_files_path / "restart_11.pkl").mkdir()

    def tearDown(self) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteDirectoryIfExisting(str(self.restart_files_path))

    def test_ParseStepFromFileName(self):
        self.assertEqual(ParseStepFromFileName(self.restart_file_name, "restart_12.pkl"), 12)
        self.assertEqual(ParseStepFromFileName("ckpt_<step>", "ckpt_3"), 3)
        for name in ["restart_5_old.pkl", "restart_.pkl", "restart_-1.pkl", "other_5.pkl", "restart_5.pkl.tmp", "restart_+5.pkl", "restart_\u00b2.pkl"]:
            self.assertIsNone(ParseStepFromFileName(self.restart_file_name, name), msg=name)
        # temp files of a pattern without suffix
        self.assertIsNone(ParseStepFromFileName("ckpt_<step>", "ckpt_3.tmp"))

    def test_ListCheckpoints(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            checkpoints = ListCheckpoints(self.restart_files_path, self.restart_file_name)
            self.assertEqual([step for step, _ in checkpoints], [2, 7, 10])
            self.assertEqual(ListCheckpoints(Path("does_not_exist"), self.restart_file_name), [])

    def test_ResolveRestartLoadStep(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            # "latest" is the highest parsed step, not the lexicographically last name
            self.assertEqual(ResolveRestartLoadStep(self.restart_files_path, self.restart_file_name, Kratos.Parameters('"latest"')), 10)
            self.assertIsNone(ResolveRestartLoadStep(Path("does_not_exist"), self.restart_file_name, Kratos.Parameters('"latest"')))
            self.assertEqual(ResolveRestartLoadStep(self.restart_files_path, self.restart_file_name, Kratos.Parameters('4')), 4)
            with self.assertRaises(RuntimeError):
                ResolveRestartLoadStep(self.restart_files_path, self.restart_file_name, Kratos.Parameters('"last"'))

    def test_PruneWithMaxFilesToKeep(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            optimization_problem = OptimizationProblem()
            optimization_problem.SetStep(13)
            restart_output = OptimizationProblemRestartOutput(Kratos.Parameters("""{
                "restart_file_name": \"""" + str(self.restart_files_path) + "/" + self.restart_file_name + """\",
                "max_files_to_keep": 2
            }"""), optimization_problem, Kratos.Model())
            restart_output.SaveRestartFile()
            # the two newest checkpoints (by step) stay; temp files, directories and unparseable names
            # are neither counted nor pruned
            remaining = sorted(path.name for path in self.restart_files_path.iterdir())
            self.assertEqual(remaining, sorted(["restart_10.pkl", "restart_13.pkl", "restart_5_old.pkl", "restart_.pkl", "restart_3.pkl.tmp", "restart_12.pkl.tmp", "restart_\u00b2.pkl", "restart_11.pkl"]))

if __name__ == "__main__":
    kratos_unittest.main()
