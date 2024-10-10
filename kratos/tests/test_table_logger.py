
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTableLogger(KratosUnittest.TestCase):

    def setUp(self):
        self.comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        self.work_folder = "test_files"
        self.size = self.comm.Size()
        self.rank = self.comm.Rank()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.comm.Barrier()
            KratosUtilities.DeleteFileIfExisting("table_logger_test.out")
            self.comm.Barrier()

    def test_TableLogger(self):
        self.file_name = "table_logger_test.out"
        # OUT file pre process
        table_logger_settings = KratosMultiphysics.Parameters("""
        {
                "file_header" : "My custom header",
                "file_name"   : "table_logger_test.out",
                "label"       : "CUSTOM",
                "columns" : [
                    {
                        "column_label" : "h1",
                        "column_header": "header 1"
                    },
                    {
                        "column_label" : "h2",
                        "column_header": "header 2"
                    },
                    {
                        "column_label" : "h3",
                        "column_header": "header 3"
                    }
                ]
        }
        """)
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            table_logger = KratosMultiphysics.LoggerTableOutput(table_logger_settings)
            KratosMultiphysics.Logger.AddOutput(table_logger)
            KratosMultiphysics.Logger.PrintInfo("h1", "C1V")
            KratosMultiphysics.Logger.PrintInfo("h2", "C2V")
            KratosMultiphysics.Logger.PrintInfo("h3", "C3V")
            KratosMultiphysics.Logger.PrintInfo("my_example.h3", "C3V")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h1", "C1V")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h2", "C2V")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h3", "C3V")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h1", "C1V2")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h2", "C2V2")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h3", "C3V2")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h1", "header 1")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h2", "header 2")
            KratosMultiphysics.Logger.PrintInfo("CUSTOM.h3", "header 3")
            if self.rank == 0:
                with open("table_logger_test.out", 'r') as f:
                    lines = f.readlines()
                    self.assertTrue(lines[0]=="My custom header\n")
                    self.assertTrue(lines[1]=="\n")
                    self.assertTrue(lines[2]==" header 1  header 2  header 3 \n")
                    self.assertTrue(lines[3]==" --------  --------  -------- \n")
                    self.assertTrue(lines[4]=="    C1V       C2V       C3V   \n")
                    self.assertTrue(lines[5]=="   C1V2      C2V2      C3V2   \n")
                    self.assertTrue(lines[6]==" header 1  header 2  header 3 \n")

if __name__ == '__main__':
    KratosUnittest.main()
