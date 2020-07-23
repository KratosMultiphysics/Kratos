import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestFileLogger(KratosUnittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """The log files are written here ONCE as otherwise every test would append
        its output to all previously added loggers."""

        #logs everything
        KM.Logger.AddOutput(KM.FileLoggerOutput("Kratos.log"))

        #logs warnings only
        warning_logger = KM.FileLoggerOutput("KratosWarning.log")
        warning_logger.SetSeverity(KM.Logger.Severity.WARNING)
        KM.Logger.AddOutput(warning_logger)

        KM.Logger.Print("TestFileLogger", "This is printed by 'Print'.")
        KM.Logger.PrintInfo("TestFileLogger", "This is printed by 'PrintInfo'.")
        KM.Logger.PrintWarning("TestFileLogger", "This is printed by 'PrintWarning'.")
        KM.Logger.PrintOnAllRanks("TestFileLogger", "This is printed by 'PrintOnAllRanks'.")
        KM.Logger.PrintInfoOnAllRanks("TestFileLogger", "This is printed by 'PrintInfoOnAllRanks'.")
        KM.Logger.PrintWarningOnAllRanks("TestFileLogger", "This is printed by 'PrintWarningOnAllRanks'.")
        KM.Logger.Flush()

    def test_file_logging(self):
        with open("Kratos.log", "r") as f:
            file_output = f.read()

        expected_output = "TestFileLogger This is printed by 'Print'. \n"
        expected_output += "TestFileLogger: This is printed by 'PrintInfo'. \n"
        expected_output += "[WARNING] TestFileLogger: This is printed by 'PrintWarning'. \n"
        expected_output += "Rank 0: TestFileLogger This is printed by 'PrintOnAllRanks'. \n"
        expected_output += "Rank 0: TestFileLogger: This is printed by 'PrintInfoOnAllRanks'. \n"
        expected_output += "[WARNING] Rank 0: TestFileLogger: This is printed by 'PrintWarningOnAllRanks'. \n"

        self.assertMultiLineEqual(file_output, expected_output)

    def test_warning_only_file_logging(self):
        with open("KratosWarning.log", "r") as f:
            file_output = f.read()

        expected_output = "[WARNING] TestFileLogger: This is printed by 'PrintWarning'. \n"
        expected_output += "[WARNING] Rank 0: TestFileLogger: This is printed by 'PrintWarningOnAllRanks'. \n"

        self.assertMultiLineEqual(file_output, expected_output)

    @classmethod
    def tearDownClass(self):
        kratos_utils.DeleteFileIfExisting("Kratos.log")
        kratos_utils.DeleteFileIfExisting("KratosWarning.log")

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
