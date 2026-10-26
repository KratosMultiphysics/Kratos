import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestCSVDatabaseIO(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.n = 5
        cls.database = {
            "bool_1": [ True, False, True, True, False],
            "bool_2": [False,  True, True, True, False],
            "int_1": [ 1,  -2, 5, 11, -5],
            "int_2": [ 10,  -200, 50, 110, 500],
            "float_1": [ 10.0,  -2.00, 5.0, 11.0, 50.0],
            "float_2": [ 10.0,  -2.00, -5.0, -11.0, 50.0],
            "str_1": [ "test1",  "test2", "test3", "test4", "test5"],
            "str_2": [ "test6",  "test7", "test8", "test9", "test10"],
            "float_3": [ 10.0,  -2.00, -5.0, -11.0, 50.0],
            "float_4": [ 10.0,  None, -5.0, None, 50.0],
        }

    @classmethod
    def tearDownClass(cls):
        with KratosUnittest.WorkFolderScope(".", __file__):
            kratos_utils.DeleteFileIfExisting("csv_database.csv")
            kratos_utils.DeleteFileIfExisting("csv_database_test.csv")

    def test_Read1(self):
        csv_database_io = Kratos.CSVDatabaseIO(
                "csv_database_ref_0.csv",
                Kratos.Testing.GetDefaultDataCommunicator(),
                row_id_name="CUSTOM_STEP",
                boolean_true_value="yes",
                echo_level=0)

        csv_database_io.Initialize()

        self.assertEqual(csv_database_io.GetRowIds(), [0,1,2,3,4])

        for step in range(self.n):
            for k, values in self.database.items():
                if k.startswith("bool"):
                    value = csv_database_io.ReadBool(step, k)
                elif k.startswith("int"):
                    value = csv_database_io.ReadInt(step, k)
                elif k.startswith("float"):
                    value = csv_database_io.ReadFloat(step, k)
                elif k.startswith("str"):
                    value = csv_database_io.ReadString(step, k)
                if (values[step] is not None):
                    self.assertEqual(value, values[step])
                else:
                    self.assertEqual(value, 0)

        csv_database_io.Finalize()

    def test_Read2(self):
        csv_database_io = Kratos.CSVDatabaseIO(
                "csv_database_no_column_info_ref.csv",
                Kratos.Testing.GetDefaultDataCommunicator(),
                row_id_name="CUSTOM_STEP",
                boolean_true_value="yes",
                echo_level=0)

        csv_database_io.Initialize()

        self.assertEqual(csv_database_io.GetRowIds(), [0,1,2,3,4])

        for step in range(self.n):
            for k, values in self.database.items():
                if k.startswith("int"):
                    value = csv_database_io.ReadFloat(step, k)
                    if (values[step] is not None):
                        self.assertEqual(value, values[step])
                    else:
                        self.assertEqual(value, 0)
                elif k.startswith("float"):
                    value = csv_database_io.ReadFloat(step, k)
                    if (values[step] is not None):
                        self.assertEqual(value, values[step])
                    else:
                        self.assertEqual(value, 0)

        csv_database_io.Finalize()

    def test_Write1(self):
        csv_database_io = Kratos.CSVDatabaseIO(
                                "csv_database_test.csv",
                                Kratos.Testing.GetDefaultDataCommunicator(),
                                "Testing output",
                                "# Testing header information \n# New line header with id = <TABLE_ID>\n",
                                row_id_name="CUSTOM_STEP",
                                write_kratos_version=False,
                                write_time_stamp=False,
                                boolean_false_value="no",
                                boolean_true_value="yes",
                                echo_level=0)

        csv_database_io.Initialize()

        for step in range(self.n):
            for k, values in self.database.items():
                if (values[step] is not None):
                    csv_database_io.Write(values[step], step, k)

        csv_database_io.Finalize()

        CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "csv_database_ref_0.csv",
                "output_file_name"      : "csv_database_test.csv",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

    def test_ReadWrite1(self):
        csv_database_io = Kratos.CSVDatabaseIO("csv_database_ref_0.csv", Kratos.Testing.GetDefaultDataCommunicator(), row_id_name="CUSTOM_STEP")

        csv_database_io.Initialize()

        with self.assertRaises(RuntimeError):
            csv_database_io.Write(1, 1, "test")

        with self.assertRaises(RuntimeError):
            csv_database_io.ReadBool(0, "CUSTOM_STEP")

        with self.assertRaises(RuntimeError):
            csv_database_io.ReadInt(10, "CUSTOM_STEP")

        csv_database_io.Finalize()

    def test_WriteRead(self):
        csv_database_io = Kratos.CSVDatabaseIO(
                                "csv_database.csv",
                                Kratos.Testing.GetDefaultDataCommunicator(),
                                "Testing output",
                                "# Testing header information \n# New line header with id = <TABLE_ID>\n",
                                row_id_name="CUSTOM_STEP",
                                write_kratos_version=False,
                                write_time_stamp=False,
                                boolean_false_value="no",
                                boolean_true_value="yes",
                                echo_level=0)

        csv_database_io.Initialize()

        csv_database_io.Write(3.0, 1, "hello")
        with self.assertRaises(RuntimeError):
            csv_database_io.ReadInt(1, "test")

        with self.assertRaises(RuntimeError):
            csv_database_io.Write(2.0, 1, "hello")

        with self.assertRaises(RuntimeError):
            csv_database_io.Write(True, 2, "hello")

        with self.assertRaises(RuntimeError):
            csv_database_io.Write(True, 2, "NEW_KEY")

        csv_database_io.Finalize()

    def test_WriteTables(self):
        csv_database_io = Kratos.CSVDatabaseIO(
                                "csv_database_test_<ID>.csv",
                                Kratos.Testing.GetDefaultDataCommunicator(),
                                "Testing output",
                                "# Testing header information \n# New line header with id = <ID>\n",
                                row_id_name="CUSTOM_STEP",
                                table_id_tag="<ID>",
                                write_kratos_version=False,
                                write_time_stamp=False,
                                boolean_false_value="no",
                                boolean_true_value="yes",
                                echo_level=0)

        for i in range(4):
            csv_database_io.Initialize(i)

            for step in range(self.n):
                for k, values in self.database.items():
                    if (values[step] is not None):
                        if isinstance(values[step], bool):
                            csv_database_io.Write(bool((values[step] * (i+1))) , step, k)
                        else:
                            csv_database_io.Write(values[step] * (i + 1), step, k)

            csv_database_io.Finalize(i)

            CompareTwoFilesCheckProcess(Kratos.Parameters("""
                {
                    "reference_file_name"   : "csv_database_ref_""" + str(i) + """.csv",
                    "output_file_name"      : "csv_database_test_""" + str(i) + """.csv",
                    "remove_output_file"    : true,
                    "comparison_type"       : "deterministic"
                }""")).Execute()

    def test_ReadTables(self):

        csv_database_io = Kratos.CSVDatabaseIO(
                "csv_database_ref_<TABLE_ID>.csv",
                Kratos.Testing.GetDefaultDataCommunicator(),
                row_id_name="CUSTOM_STEP",
                boolean_true_value="yes",
                echo_level=0)

        for i in range(4):
            csv_database_io.Initialize(i)

            self.assertEqual(csv_database_io.GetRowIds(), [0,1,2,3,4])

            for step in range(self.n):
                for k, values in self.database.items():
                    if k.startswith("bool"):
                        value = csv_database_io.ReadBool(step, k)
                    elif k.startswith("int"):
                        value = csv_database_io.ReadInt(step, k)
                    elif k.startswith("float"):
                        value = csv_database_io.ReadFloat(step, k)
                    elif k.startswith("str"):
                        value = csv_database_io.ReadString(step, k)
                    if (values[step] is not None):
                        if isinstance(value, bool):
                            self.assertEqual(value, bool((values[step] * (i+1))))
                        else:
                            self.assertEqual(value, values[step] * (i + 1))
                    else:
                        self.assertEqual(value, 0)

            csv_database_io.Finalize(i)

if __name__ == '__main__':
    KratosUnittest.main()
