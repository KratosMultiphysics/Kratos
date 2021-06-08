import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting

class TestFileNameDataCollector(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = KratosMultiphysics.Model()
        cls.model_part = cls.model.CreateModelPart("test_model_part")
        cls.model_part.ProcessInfo[KratosMultiphysics.TIME] = 1e-3
        cls.model_part.ProcessInfo[KratosMultiphysics.STEP] = 3
        cls.sub_model_part = cls.model_part.CreateSubModelPart("test_sub_model_part")

    def test_GetFileName(self):
        self.assertEqual(
            KratosMultiphysics.FileNameDataCollector(
                self.model_part,
                "<model_part_name>-<rank>-<step>-<time>.h5",
                {"<time>":"%0.4f"}
            ).GetFileName(), "test_model_part-0-3-0.0010.h5")

        self.assertEqual(
            KratosMultiphysics.FileNameDataCollector(
                self.model_part,
                "<model_part_name>-<rank>-<time>.h5",
                {"<time>":"0.4e"}
            ).GetFileName(), "test_model_part-0-1.0000e-03.h5")

        self.assertEqual(
            KratosMultiphysics.FileNameDataCollector(
                self.sub_model_part,
                "<model_part_full_name>-<time>.<rank>--<time>.h5",
                {
                    "<time>":"0.4E",
                    "<rank>": "%4d"
                }
            ).GetFileName(), "test_model_part.test_sub_model_part-1.0000E-03.   0--1.0000E-03.h5")

    def test_RetrieveFileNameData(self):
        file_name_data_type = KratosMultiphysics.FileNameDataCollector.FileNameData

        check_list_file_name_data = []
        def add_file_name(file_name, rank, step, time_step):
            check_list_file_name_data.append([file_name, file_name_data_type(file_name, rank, step, time_step)])

        add_file_name("test_model_part-10.0.h5", -1, -1, 10.0)           # 0
        add_file_name("10.0-test_model_part.h5", -1, -1, 10.0)           # 1
        add_file_name("10.0-1-test_model_part.h5", 1, -1, 10.0)          # 2
        add_file_name("test_model_part-1-3-10.0.h5", 1, 3, 10.0)         # 3
        add_file_name("test_model_part-1-10.0.h5", 1, -1, 10.0)          # 4
        add_file_name("test_model_part-1-1e-4.h5", 1, -1, 1e-4)          # 5
        add_file_name("test_model_part-1-0.05.h5", 1, -1, 0.05)          # 6
        add_file_name("test_model_part_10.0.h5", -1, -1, 10.0)           # 7
        add_file_name("test_model_part_1_3_10.0.h5", 1, 3, 10.0)         # 8
        add_file_name("test_model_part_1_10.0.h5", 1, -1, 10.0)          # 9
        add_file_name("test_model_part_1_1e-4.h5", 1, -1, 1e-4)          # 10
        add_file_name("test_model_part_1_0.05.h5", 1, -1, 0.05)          # 11
        add_file_name("test_model_part_10.0.h5", -1, -1, 10.0)           # 12
        add_file_name("test_model_part_1_3-10.0.h5", 1, 3, 10.0)         # 13
        add_file_name("test_model_part_1-10.0.h5", 1, -1, 10.0)          # 14
        add_file_name("test_model_part_1-1e-4.h5", 1, -1, 1e-4)          # 15
        add_file_name("test_model_part_1-0.05.h5", 1, -1, 0.05)          # 16
        add_file_name("test_model_part_stat_1_0.05.h5", 1, -1, 0.05)     # 17
        add_file_name("test_model_part_stat    1_0.05.h5", 1, -1, 0.05)  # 18

        def check_file_names(list_of_matching_indices, pattern):
            file_name_data_collector = KratosMultiphysics.FileNameDataCollector(self.model_part, pattern, {})
            f_data = file_name_data_type()
            for index, check_data in enumerate(check_list_file_name_data):
                file_name_without_path = check_data[0]
                if (file_name_data_collector.RetrieveFileNameData(f_data, file_name_without_path)):
                    self.assertTrue(index in list_of_matching_indices)
                    check_data[1].SetFileName(file_name_data_collector.GetPath() + "/" + file_name_without_path)
                    self.assertEqual(f_data, check_data[1])
                else:
                    self.assertFalse(index in list_of_matching_indices)

        check_file_names([0], "/test_folder/<model_part_full_name>/<model_part_name>-<time>.h5")
        check_file_names([1], "/test_folder/<model_part_full_name>/<time>-<model_part_name>.h5")
        check_file_names([2], "/test_folder/<model_part_full_name>/<time>-<rank>-<model_part_name>.h5")
        check_file_names([3], "/test_folder/<model_part_full_name>/<model_part_name>-<rank>-<step>-<time>.h5")
        check_file_names([4, 5, 6], "/test_folder/<model_part_full_name>/<model_part_name>-<rank>-<time>.h5")
        check_file_names([7, 12], "/test_folder/<model_part_full_name>/<model_part_name>_<time>.h5")
        check_file_names([8], "/test_folder/<model_part_full_name>/<model_part_name>_<rank>_<step>_<time>.h5")
        check_file_names([9, 10, 11], "/test_folder/<model_part_full_name>/<model_part_name>_<rank>_<time>.h5")
        check_file_names([13], "/test_folder/<model_part_full_name>/<model_part_name>_<rank>_<step>-<time>.h5")
        check_file_names([14, 15, 16], "/test_folder/<model_part_full_name>/<model_part_name>_<rank>-<time>.h5")
        check_file_names([17], "/test_folder/<model_part_full_name>/test_model_part_stat_<rank>_<time>.h5")
        check_file_names([18], "/test_folder/<model_part_full_name>/test_model_part_stat <rank>_<time>.h5")

    def test_SortListOfFileNameData(self):
        list_of_file_names = [
            "test_model_part-1-5-20.0.h5",  # 0
            "test_model_part-1-3-11.0.h5",  # 1
            "test_model_part-1-4-12.0.h5",  # 2
            "test_model_part-2-3-10.0.h5",  # 3
            "test_model_part-2-4-15.0.h5",  # 4
            "test_model_part-2-5-12.0.h5",  # 5
            "test_model_part-3-3-15.0.h5",  # 6
            "test_model_part-3-4-11.0.h5",  # 7
            "test_model_part-3-5-01.0.h5"   # 8
        ]

        file_name_data_collector = KratosMultiphysics.FileNameDataCollector(self.model_part, "<model_part_name>-<rank>-<step>-<time>.h5", {})

        list_of_file_name_data = []
        for file_name in list_of_file_names:
            f_data = KratosMultiphysics.FileNameDataCollector.FileNameData()
            self.assertTrue(file_name_data_collector.RetrieveFileNameData(f_data, file_name))
            list_of_file_name_data.append(f_data)

        def check_sorted_list(sorting_order, sorted_indices_list):
            for index, sorted_file_name_data in enumerate(KratosMultiphysics.FileNameDataCollector.GetSortedListOfFileNameData(list_of_file_name_data, sorting_order)):
                f_data = KratosMultiphysics.FileNameDataCollector.FileNameData()
                self.assertTrue(file_name_data_collector.RetrieveFileNameData(f_data, list_of_file_names[sorted_indices_list[index]]))
                self.assertEqual(sorted_file_name_data, f_data)

        check_sorted_list(["<rank>", "<step>", "<time>"], [1, 2, 0, 3, 4, 5, 6, 7, 8])
        check_sorted_list(["<step>", "<time>"],           [3, 1, 6, 7, 2, 4, 8, 5, 0])
        check_sorted_list(["<time>", "<rank>"],           [8, 3, 1, 7, 2, 5, 4, 6, 0])

    def test_GetSortedFileNamesList(self):
        dir_name = "GSFNL_test_model_part"
        self.addCleanup(lambda: DeleteDirectoryIfExisting(dir_name))
        list_of_file_names = [
            os.path.join(dir_name, "test_model_part-1-4-1e-3.vtk"),     # 0
            os.path.join(dir_name, "test_model_part-1-3-1e-2.h5"),      # 1
            os.path.join(dir_name, "test_model_part-1-2-1e-1.vtk"),     # 2
            os.path.join(dir_name, "test_model_part-1-5-1e+1.vtk"),     # 3
            os.path.join(dir_name, "test_model_part-2-2-1e-3.vtk"),     # 4
            os.path.join(dir_name, "test_model_part-2-3-1e-2.vtk"),     # 5
            os.path.join(dir_name, "test_model_part-2-4-1e-1.vtk"),     # 6
            os.path.join(dir_name, "test_model_part-3-4-1e-3.vtk"),     # 7
            os.path.join(dir_name, "test_model_part-3-3-1e-2.vtk"),     # 8
            os.path.join(dir_name, "test_model_part-3-2-1e-1.vtk"),     # 9
            os.path.join(dir_name, "test_model_part-3-1e-1.vtk"),       # 10
            os.path.join(dir_name, "test_model_part-1e-1.vtk"),         # 11
            os.path.join(dir_name, "test_model_part-1e-1-3.vtk")        # 12
        ]

        # check and create the dir
        self.assertFalse(os.path.isdir(dir_name))
        os.mkdir(dir_name)
        self.assertTrue(os.path.isdir(dir_name))

        for file_name in list_of_file_names:
            file_out = open(file_name, "w")
            file_out.close()
            self.assertTrue(os.path.isfile(file_name))

        file_name_data_collector = KratosMultiphysics.FileNameDataCollector(self.model_part, "GSFNL_<model_part_full_name>/<model_part_name>-<rank>-<step>-<time>.vtk", {})

        def check_sorted_list(sorting_order, sorted_indices_list):
            for index, sorted_file_name in enumerate(file_name_data_collector.GetSortedFileNamesList(sorting_order)):
                self.assertEqual(os.path.split(sorted_file_name), os.path.split(list_of_file_names[sorted_indices_list[index]]))

        check_sorted_list(["<rank>", "<step>"], [2, 0, 3, 4, 5, 6, 9, 8, 7])
        check_sorted_list(["<step>", "<rank>"], [2, 4, 9, 5, 8, 0, 6, 7, 3])
        check_sorted_list(["<time>", "<rank>"], [0, 4, 7, 5, 8, 2, 6, 9, 3])

    def test_ExtractFileNamePattern(self):
        extract_pattern = KratosMultiphysics.FileNameDataCollector.ExtractFileNamePattern

        self.assertEqual((True, "/test/TestModelParte-<time>.h<step>"), extract_pattern("/test/TestModelParte-0.00.h5", ["<time>", "<step>"]))
        self.assertEqual((True, "TestModelParte-<time>.h5"), extract_pattern("TestModelParte-0.00.h5", ["<time>"]))
        self.assertEqual((True, "TestModelParte-<step>.00.h5"), extract_pattern("TestModelParte-0.00.h5", ["<step>"]))
        self.assertEqual((True, "TestModelPart_10-<time>.h5"), extract_pattern("TestModelPart_10-0.00.h5", ["<skip_float>","<time>"]))
        self.assertEqual((True, "TestModelPart_10.0-<time>.h5"), extract_pattern("TestModelPart_10.0-0.00.h5", ["<skip_float>","<time>"]))
        self.assertEqual((True, "TestModelPart_10.<time>-0.00.h5"), extract_pattern("TestModelPart_10.0-0.00.h5", ["<skip_int>","<time>"]))
        self.assertEqual((False, "TestModelPart_-0.00.vtk"), extract_pattern("TestModelPart_-0.00.vtk", ["<skip_float>", "<time>"]))
        self.assertEqual((True, "TestModelParte-0-<time>.vtk"), extract_pattern("TestModelParte-0-1e-2.vtk", ["<skip_float>", "<time>"]))

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

