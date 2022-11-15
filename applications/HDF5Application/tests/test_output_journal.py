# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.HDF5Application.output_journal import OutputJournal

# --- STD Imports ---
import pathlib


class TestOutputJournal(KratosUnittest.TestCase):

    @property
    def test_directory_path(self) -> pathlib.Path:
        return pathlib.Path(__file__).absolute().parent / "test_output_journal"


    def setUp(self) -> None:
        DeleteDirectoryIfExisting(str(self.test_directory_path))
        self.test_directory_path.mkdir(parents = True)


    def tearDown(self) -> None:
        DeleteDirectoryIfExisting(str(self.test_directory_path))


    def test_FileManagement(self) -> None:
        file_path_generator = lambda step: self.test_directory_path / f"test_file{'_' * step}.suffix"

        def extractor(
                model: KratosMultiphysics.Model,
                output_base: KratosMultiphysics.Parameters = KratosMultiphysics.Parameters()
                ) -> KratosMultiphysics.Parameters:
            output = output_base.Clone()
            step = model.GetModelPart("test").ProcessInfo[KratosMultiphysics.STEP]
            output.AddString("file_path", str(file_path_generator(step)))
            output.AddInt("step", step)
            return output

        # Create a dummy model
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("test")

        # Construct and populate a journal
        journal_path = self.test_directory_path / "test_output_journal.txt"
        journal = OutputJournal(str(journal_path), extractor)

        for step in range(10):
            file_path = file_path_generator(step)
            file_path.touch() # <== create the file

            model_part.ProcessInfo[KratosMultiphysics.STEP] = step
            journal.Push(model)

        # Check whether all the required files were created successfully
        self.assertTrue((self.test_directory_path / "test_output_journal.txt").is_file())

        for step in range(10):
            self.assertTrue(file_path_generator(step).is_file())

        # Erase odd steps
        journal.EraseIf(lambda entry: bool(entry["step"].GetInt() % 2))

        # Check files and journal state
        for step in range(10):
            file_path = file_path_generator(step)
            file_exists = file_path.is_file()
            file_in_journal = any(file_path == pathlib.Path(entry["file_path"].GetString()) for entry in journal)
            if step % 2:
                self.assertFalse(file_exists)
                self.assertFalse(file_in_journal)
            else:
                self.assertTrue(file_exists)
                self.assertTrue(file_in_journal)

        # Update entry at step 6
        # ==> the entry at step 6 should be updated, not duplicated
        entry_base_at_step_6 = KratosMultiphysics.Parameters("""{"extra" : "stuff"}""")
        journal.SetExtractor(lambda model: extractor(model, output_base = entry_base_at_step_6))
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 6
        journal.Push(model)

        self.assertEqual(sum(entry["step"].GetInt() == 6 for entry in journal), 1) # exactly one entry at step 6
        entry_at_step_6 = next((entry for entry in journal if entry["step"].GetInt() == 6), None)

        self.assertNotEqual(entry_at_step_6, None)
        self.assertTrue(entry_at_step_6.Has("extra"))
        self.assertEqual(entry_at_step_6["extra"].GetString(), "stuff")

        # Clear the journal ==> no files should be left
        journal.Clear()
        self.assertFalse((self.test_directory_path / "test_output_journal.txt").is_file())
        for step in range(10):
            self.assertFalse(file_path_generator(step).is_file())


    def test_InvalidEntry(self) -> None:
        model = KratosMultiphysics.Model()
        extractor = lambda model: KratosMultiphysics.Parameters("""{"not" : "what", "you'd" : "expect"}""")
        journal = OutputJournal(self.test_directory_path / "invalid_journal.jrn", extractor)
        self.assertRaises(Exception, lambda: journal.Push(model))


if __name__ == "__main__":
    KratosUnittest.main()
