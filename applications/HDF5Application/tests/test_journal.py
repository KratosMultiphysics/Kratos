# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.HDF5Application as HDF5
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

# --- STD Imports ---
import pathlib


class TestJournal(KratosUnittest.TestCase):

    @property
    def test_file_path(self) -> pathlib.Path:
        return pathlib.Path("test_journal.log")

    def setUp(self) -> None:
        DeleteFileIfExisting(str(self.test_file_path))

    def tearDown(self) -> None:
        DeleteFileIfExisting(str(self.test_file_path))

    def test_Journal(self) -> None:
        journal = HDF5.Journal(self.test_file_path)
        model = KratosMultiphysics.Model()

        extractor = lambda model: KratosMultiphysics.Parameters('["1st","2nd"]')
        journal.SetExtractor(extractor)

        journal.Clear()
        self.assertEqual(len(journal), 0)

        for i in range(1, 3):
            journal.Push(model)
            self.assertEqual(len(journal), i)
            for item in journal:
                self.assertTrue(isinstance(item, KratosMultiphysics.Parameters))
                self.assertEqual(item.WriteJsonString(), '["1st","2nd"]')

        journal.Clear()
        self.assertEqual(len(journal), 0)

    def test_Extractor(self) -> None:
        def extractor(model: KratosMultiphysics.Model) -> KratosMultiphysics.Parameters:
            model_part = model.GetModelPart("test")
            process_info = model_part.ProcessInfo
            return KratosMultiphysics.Parameters(f"""[
                "example/relative/directory/example_file_name_{process_info[KratosMultiphysics.STEP]}.h5",
                {process_info[KratosMultiphysics.STEP]},
                {process_info[KratosMultiphysics.TIME]},
                true
            ]""")

        journal = HDF5.Journal(self.test_file_path, extractor)
        journal.Clear()

        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("test")

        for i_step in range(10):
            model_part.CloneTimeStep()
            model_part.ProcessInfo[KratosMultiphysics.STEP] = i_step
            model_part.ProcessInfo[KratosMultiphysics.TIME] = i_step / 10.0
            journal.Push(model)

        self.assertEqual(len(journal), 10)
        for index, item in enumerate(journal):
            self.assertTrue(isinstance(item, KratosMultiphysics.Parameters))
            self.assertTrue(item.IsArray())
            self.assertEqual(item.size(), 4)

            self.assertTrue(item[0].IsString())
            self.assertEqual(item[0].GetString(), f"example/relative/directory/example_file_name_{index}.h5")

            self.assertTrue(item[1].IsInt())
            self.assertEqual(item[1].GetInt(), index)

            self.assertTrue(item[2].IsDouble())
            self.assertEqual(item[2].GetDouble(), index / 10.0)

            self.assertTrue(item[3].IsBool())
            self.assertEqual(item[3].GetBool(), True)

        # Check whether python objects are protected in a multithreaded environment
        journal.Clear()
        self.assertFalse(len(journal))
        KratosMultiphysics.HDF5Application.TestingUtilities.TestJournal(model, journal)
        self.assertTrue(len(journal))


if __name__ == "__main__":
    KratosUnittest.main()
