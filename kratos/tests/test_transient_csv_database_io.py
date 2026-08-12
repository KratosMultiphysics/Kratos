import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestTransientCSVDatabaseIO(KratosUnittest.TestCase):
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

    def test_Read1(self):
        params = Kratos.Parameters("""{
            "file_name"           : "transient_csv_database_ref.csv",
            "write_kratos_version": true,
            "write_time_stamp"    : true,
            "echo_level"          : 0,
            "format_settings"     : {
                "int_length"     : 7,
                "float_precision": 9,
                "bool_values"    : ["no", "yes"],
                "string_length"  : 30
            }
        }""")
        transient_csv_database_io = Kratos.TransientCSVDatabaseIO(params, Kratos.Testing.GetDefaultDataCommunicator())
        transient_csv_database_io.Initialize()

        for step in range(self.n):
            for k, values in self.database.items():
                if k.startswith("bool"):
                    value = transient_csv_database_io.ReadBool(step, k)
                elif k.startswith("int"):
                    value = transient_csv_database_io.ReadInt(step, k)
                elif k.startswith("float"):
                    value = transient_csv_database_io.ReadFloat(step, k)
                elif k.startswith("str"):
                    value = transient_csv_database_io.ReadString(step, k)
                if (values[step] is not None):
                    self.assertEqual(value, values[step])
                else:
                    self.assertEqual(value, 0)

        transient_csv_database_io.Finalize()

    def test_Read2(self):
        params = Kratos.Parameters("""{
            "file_name"           : "transient_csv_database_no_column_info_ref.csv",
            "write_kratos_version": true,
            "write_time_stamp"    : true,
            "echo_level"          : 0,
            "format_settings"     : {
                "int_length"     : 7,
                "float_precision": 9,
                "bool_values"    : ["no", "yes"],
                "string_length"  : 30
            }
        }""")
        transient_csv_database_io = Kratos.TransientCSVDatabaseIO(params, Kratos.Testing.GetDefaultDataCommunicator())
        transient_csv_database_io.Initialize()

        for step in range(self.n):
            for k, values in self.database.items():
                if k.startswith("int"):
                    value = transient_csv_database_io.ReadFloat(step, k)
                    if (values[step] is not None):
                        self.assertEqual(value, values[step])
                    else:
                        self.assertEqual(value, 0)
                elif k.startswith("float"):
                    value = transient_csv_database_io.ReadFloat(step, k)
                    if (values[step] is not None):
                        self.assertEqual(value, values[step])
                    else:
                        self.assertEqual(value, 0)

        transient_csv_database_io.Finalize()

    def test_Write1(self):
        params = Kratos.Parameters("""{
            "file_name"           : "transient_csv_database.csv",
            "write_kratos_version": false,
            "write_time_stamp"    : false,
            "format_settings"     : {
                "int_length"     : 7,
                "float_precision": 9,
                "bool_values"    : ["no", "yes"],
                "string_length"  : 30
            }
        }""")
        transient_csv_database_io = Kratos.TransientCSVDatabaseIO(params, Kratos.Testing.GetDefaultDataCommunicator())
        transient_csv_database_io.Initialize()
        transient_csv_database_io.SetHeaderInformation("Testing header information\nSecond line of testing information")

        for step in range(self.n):
            for k, values in self.database.items():
                if (values[step] is not None):
                    transient_csv_database_io.Write(values[step], step, k)

        transient_csv_database_io.Finalize()

        CompareTwoFilesCheckProcess(Kratos.Parameters("""
            {
                "reference_file_name"   : "transient_csv_database_ref.csv",
                "output_file_name"      : "transient_csv_database.csv",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic",
                "tolerance"             : 1e-6,
                "relative_tolerance"    : 1e-6
            }""")).Execute()

    def test_ReadWrite1(self):
        params = Kratos.Parameters("""{
            "file_name"           : "transient_csv_database.csv",
            "write_kratos_version": false,
            "write_time_stamp"    : false,
            "format_settings"     : {
                "int_length"     : 7,
                "float_precision": 9,
                "bool_values"    : ["no", "yes"],
                "string_length"  : 30
            }
        }""")
        transient_csv_database_io = Kratos.TransientCSVDatabaseIO(params, Kratos.Testing.GetDefaultDataCommunicator())
        transient_csv_database_io.Initialize()
        transient_csv_database_io.SetHeaderInformation("Testing header information\nSecond line of testing information")

        transient_csv_database_io.Write(1, 1, "test")
        transient_csv_database_io.Write(3.0, 1, "hello")
        with self.assertRaises(RuntimeError):
            transient_csv_database_io.ReadInt(1, "test")

        with self.assertRaises(RuntimeError):
            transient_csv_database_io.Write(2.0, 1, "test")

        with self.assertRaises(RuntimeError):
            transient_csv_database_io.Write(True, 2, "NEW_KEY")

if __name__ == '__main__':
    KratosUnittest.main()
