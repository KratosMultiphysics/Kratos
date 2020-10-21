//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Vicente Mataix Ferrándiz
//

// System includes
#include <fstream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "testing/testing.h"
#include "includes/kratos_filesystem.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(FileSystemExists, KratosCoreFastSuite)
{
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists("StupidNameShouldNotExist"));

    const std::string file_name("dummy_file.txt");
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
    std::ofstream output_file;
    output_file.open(file_name);
    KRATOS_CHECK(output_file.good());
    output_file.close();

    KRATOS_CHECK(Kratos::filesystem::exists(file_name));

    KRATOS_CHECK(Kratos::filesystem::remove(file_name));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
}

KRATOS_TEST_CASE_IN_SUITE(FileSystemJoinPaths, KratosCoreFastSuite)
{
    std::vector<std::string> paths_1 {"eee", "ccc", "gt"};
    KRATOS_CHECK_STRING_EQUAL(Kratos::FilesystemExtensions::JoinPaths(paths_1), "eee/ccc/gt");

    std::vector<std::string> paths_2;
    KRATOS_CHECK_STRING_EQUAL(Kratos::FilesystemExtensions::JoinPaths(paths_2), "");
}

KRATOS_TEST_CASE_IN_SUITE(FileSystemJoinEmptyPaths, KratosCoreFastSuite)
{
    std::vector<std::string> paths_1 {"eee", "", "gt"};
    KRATOS_CHECK_STRING_EQUAL(Kratos::FilesystemExtensions::JoinPaths(paths_1), "eee/gt");

    std::vector<std::string> paths_2 {"", "", "gt"};
    KRATOS_CHECK_STRING_EQUAL(Kratos::FilesystemExtensions::JoinPaths(paths_2), "gt");
}

KRATOS_TEST_CASE_IN_SUITE(FileSystemIsRegularFile, KratosCoreFastSuite)
{
    const std::string file_name("dummy_file.txt");
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));

    std::ofstream output_file;
    output_file.open(file_name);
    KRATOS_CHECK(output_file.good());
    output_file.close();

    KRATOS_CHECK(Kratos::filesystem::exists(file_name));

    KRATOS_CHECK(Kratos::filesystem::is_regular_file(file_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::is_directory(file_name));

    KRATOS_CHECK(Kratos::filesystem::remove(file_name));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
}

KRATOS_TEST_CASE_IN_SUITE(FileSystemIsDirectory, KratosCoreFastSuite)
{
    const std::string dir_name("MyCustomDir");
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name));

    KRATOS_CHECK(Kratos::filesystem::create_directory(dir_name));

    KRATOS_CHECK(Kratos::filesystem::exists(dir_name));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::is_regular_file(dir_name));
    KRATOS_CHECK(Kratos::filesystem::is_directory(dir_name));

    KRATOS_CHECK(Kratos::filesystem::remove(dir_name));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name));
}

KRATOS_TEST_CASE_IN_SUITE(FileSystemSubDirectories, KratosCoreFastSuite)
{
    const std::string base_dir_name("MyCustomDir2");
    const std::string sub_dir_name("TheSubDir");
    const std::string sub_sub_dir_name("TheSubSubDir");

    const std::string full_dir_name = Kratos::FilesystemExtensions::JoinPaths({base_dir_name, sub_dir_name, sub_sub_dir_name});

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(base_dir_name));

    KRATOS_CHECK(Kratos::filesystem::create_directory(base_dir_name)); // first create only the base-dir

    KRATOS_CHECK(Kratos::filesystem::exists(base_dir_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(full_dir_name));

    KRATOS_CHECK(Kratos::filesystem::create_directories(full_dir_name)); // now create the entire hierarchie of directories

    KRATOS_CHECK(Kratos::filesystem::exists(base_dir_name));
    KRATOS_CHECK(Kratos::filesystem::exists(full_dir_name));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::is_regular_file(base_dir_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::is_regular_file(full_dir_name));
    KRATOS_CHECK(Kratos::filesystem::is_directory(base_dir_name));
    KRATOS_CHECK(Kratos::filesystem::is_directory(full_dir_name));

    // creating a file in a subdir to make sure the deletion also works
    const std::string file_name(Kratos::FilesystemExtensions::JoinPaths({base_dir_name, sub_dir_name, "dummy_file.txt"}));
    std::ofstream output_file;
    output_file.open(file_name);
    KRATOS_CHECK(output_file.good());
    output_file.close();

    Kratos::filesystem::remove_all(base_dir_name);

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(base_dir_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(full_dir_name));
}

KRATOS_TEST_CASE_IN_SUITE(FileSystemRenameFile, KratosCoreFastSuite)
{
    const std::string file_name("dummy_file.txt");
    const std::string file_name_new("dummy_file_renamed.txt");
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name_new));

    std::ofstream output_file;
    output_file.open(file_name);
    KRATOS_CHECK(output_file.good());
    output_file.close();

    KRATOS_CHECK(Kratos::filesystem::exists(file_name));

    Kratos::filesystem::rename(file_name, file_name_new);

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
    KRATOS_CHECK(Kratos::filesystem::exists(file_name_new));

    KRATOS_CHECK(Kratos::filesystem::remove(file_name_new));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name_new));
}

KRATOS_TEST_CASE_IN_SUITE(FileSystemRenameDirectory, KratosCoreFastSuite)
{
    const std::string dir_name("DirNameOld");
    const std::string dir_name_new("DirNameRenamed");

    const std::string raw_file_name("dummy_file.txt");
    const std::string file_name(Kratos::FilesystemExtensions::JoinPaths({dir_name, raw_file_name}));
    const std::string file_name_new(Kratos::FilesystemExtensions::JoinPaths({dir_name_new, raw_file_name}));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name_new));

    KRATOS_CHECK(Kratos::filesystem::create_directory(dir_name));

    std::ofstream output_file;
    output_file.open(file_name);
    KRATOS_CHECK(output_file.good());
    output_file.close();

    KRATOS_CHECK(Kratos::filesystem::exists(file_name));

    Kratos::filesystem::rename(dir_name, dir_name_new);

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
    KRATOS_CHECK(Kratos::filesystem::exists(dir_name_new));
    KRATOS_CHECK(Kratos::filesystem::exists(file_name_new));

    Kratos::filesystem::remove_all(dir_name_new);

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name_new));
}

KRATOS_TEST_CASE_IN_SUITE(FileNameInformationCollector_RetrieveFileNameInformation, KratosCoreFastSuite)
{
    Model model;
    const auto& model_part = model.CreateModelPart("test_model_part");

    using fic =  filesystem::FileNameInformationCollector;
    using file_name_data = filesystem::FileNameData;
    file_name_data f_data;

    const std::vector<std::pair<std::string, file_name_data>>& list_of_checking_file_names = {
        std::make_pair("test_model_part-10.0.h5", file_name_data{-1, -1, 10.0}),        // 0
        std::make_pair("10.0-test_model_part.h5", file_name_data{-1, -1, 10.0}),        // 1
        std::make_pair("10.0-1-test_model_part.h5", file_name_data{1, -1, 10.0}),       // 2
        std::make_pair("test_model_part-1-3-10.0.h5", file_name_data{1, 3, 10.0}),      // 3
        std::make_pair("test_model_part-1-10.0.h5", file_name_data{1, -1, 10.0}),       // 4
        std::make_pair("test_model_part-1-1e-4.h5", file_name_data{1, -1, 1e-4}),       // 5
        std::make_pair("test_model_part-1-0.05.h5", file_name_data{1, -1, 0.05}),       // 6
        std::make_pair("test_model_part_10.0.h5", file_name_data{-1, -1, 10.0}),        // 7
        std::make_pair("test_model_part_1_3_10.0.h5", file_name_data{1, 3, 10.0}),      // 8
        std::make_pair("test_model_part_1_10.0.h5", file_name_data{1, -1, 10.0}),       // 9
        std::make_pair("test_model_part_1_1e-4.h5", file_name_data{1, -1, 1e-4}),       // 10
        std::make_pair("test_model_part_1_0.05.h5", file_name_data{1, -1, 0.05}),       // 11
        std::make_pair("test_model_part_10.0.h5", file_name_data{-1, -1, 10.0}),        // 12
        std::make_pair("test_model_part_1_3-10.0.h5", file_name_data{1, 3, 10.0}),      // 13
        std::make_pair("test_model_part_1-10.0.h5", file_name_data{1, -1, 10.0}),       // 14
        std::make_pair("test_model_part_1-1e-4.h5", file_name_data{1, -1, 1e-4}),       // 15
        std::make_pair("test_model_part_1-0.05.h5", file_name_data{1, -1, 0.05}),       // 16
        std::make_pair("test_model_part_stat_1_0.05.h5", file_name_data{1, -1, 0.05})   // 17
    };

    const auto& check_file_names = [&](
        const std::vector<int>& rIsMatchingIndexList,
        const std::string& rPattern)
    {
        KRATOS_TRY

        fic pattern(model_part, rPattern);

        for (std::size_t i = 0; i < list_of_checking_file_names.size(); ++i) {
            const auto& file_name_data = list_of_checking_file_names[i].second;
            const auto& file_name = list_of_checking_file_names[i].first;
            if (pattern.RetrieveFileNameInformation(f_data, file_name)) {
                KRATOS_CHECK_NOT_EQUAL(std::find(rIsMatchingIndexList.begin(), rIsMatchingIndexList.end(), i), rIsMatchingIndexList.end());
                KRATOS_CHECK_EQUAL(f_data.Rank, file_name_data.Rank);
                KRATOS_CHECK_EQUAL(f_data.Step, file_name_data.Step);
                KRATOS_CHECK_EQUAL(f_data.TimeStep, file_name_data.TimeStep);
            } else {
                KRATOS_CHECK_EQUAL(std::find(rIsMatchingIndexList.begin(), rIsMatchingIndexList.end(), i), rIsMatchingIndexList.end());
            }
        }

        KRATOS_CATCH("");
    };

    check_file_names({0}, "/test_folder/<model_part_full_name>/<model_part_name>-<time>.h5");
    check_file_names({1}, "/test_folder/<model_part_full_name>/<time>-<model_part_name>.h5");
    check_file_names({2}, "/test_folder/<model_part_full_name>/<time>-<rank>-<model_part_name>.h5");
    check_file_names({3}, "/test_folder/<model_part_full_name>/<model_part_name>-<rank>-<step>-<time>.h5");
    check_file_names({4, 5, 6}, "/test_folder/<model_part_full_name>/<model_part_name>-<rank>-<time>.h5");
    check_file_names({7, 12}, "/test_folder/<model_part_full_name>/<model_part_name>_<time>.h5");
    check_file_names({8}, "/test_folder/<model_part_full_name>/<model_part_name>_<rank>_<step>_<time>.h5");
    check_file_names({9, 10, 11}, "/test_folder/<model_part_full_name>/<model_part_name>_<rank>_<time>.h5");
    check_file_names({13}, "/test_folder/<model_part_full_name>/<model_part_name>_<rank>_<step>-<time>.h5");
    check_file_names({14, 15, 16}, "/test_folder/<model_part_full_name>/<model_part_name>_<rank>-<time>.h5");
    check_file_names({17}, "/test_folder/<model_part_full_name>/test_model_part_stat_<rank>_<time>.h5");
}

KRATOS_TEST_CASE_IN_SUITE(FileNameInformationCollector_SortListOfFileNameData, KratosCoreFastSuite)
{
    Model model;
    const auto& model_part = model.CreateModelPart("test_model_part");

    using fic = filesystem::FileNameInformationCollector;
    using file_name_data = filesystem::FileNameData;

    const std::vector<std::string>& list_of_file_names = {
        "test_model_part-1-5-20.0.h5",  // 0
        "test_model_part-1-3-11.0.h5",  // 1
        "test_model_part-1-4-12.0.h5",  // 2
        "test_model_part-2-3-10.0.h5",  // 3
        "test_model_part-2-4-15.0.h5",  // 4
        "test_model_part-2-5-12.0.h5",  // 5
        "test_model_part-3-3-15.0.h5",  // 6
        "test_model_part-3-4-11.0.h5",  // 7
        "test_model_part-3-5-1.0.h5",   // 8
    };

    std::vector<file_name_data> f_data;
    const auto& get_file_name_data = [&](const std::string& rPattern) {
        f_data.clear();
        fic pattern(model_part, rPattern);
        for (const auto& r_file_name : list_of_file_names) {
            file_name_data current_f_data;
            pattern.RetrieveFileNameInformation(current_f_data, r_file_name);
            current_f_data.FileName = r_file_name;
            f_data.push_back(current_f_data);
        }
    };

    const auto& check_sorted_list = [&](const std::string& rPattern, const std::vector<std::string>& rSortingOrder, const std::vector<int>& rCheckIndices) {
        KRATOS_TRY

        get_file_name_data(rPattern);
        fic::SortListOfFileNameData(f_data, rSortingOrder);

        for (std::size_t i = 0; i < f_data.size(); ++i) {
            KRATOS_CHECK_EQUAL(f_data[i].FileName, list_of_file_names[rCheckIndices[i]]);
        }

        KRATOS_CATCH("");
    };

    check_sorted_list("<model_part_name>-<rank>-<step>-<time>.h5", {"<rank>", "<step>", "<time>"}, {1, 2, 0, 3, 4, 5, 6, 7, 8});
    check_sorted_list("<model_part_name>-<rank>-<step>-<time>.h5", {"<step>", "<time>"},           {3, 1, 6, 7, 2, 4, 8, 5, 0});
    check_sorted_list("<model_part_name>-<rank>-<step>-<time>.h5", {"<time>", "<rank>"},           {8, 3, 1, 7, 2, 5, 4, 6, 0});

}

KRATOS_TEST_CASE_IN_SUITE(FileNameInformationCollector_GetFileName, KratosCoreFastSuite)
{
    KRATOS_TRY

    Model model;
    auto& model_part = model.CreateModelPart("test_model_part");
    model_part.GetProcessInfo().SetValue(TIME, 1e-3);
    model_part.GetProcessInfo().SetValue(STEP, 3);
    using fic = filesystem::FileNameInformationCollector;

    KRATOS_CHECK_EQUAL(fic(model_part, "<model_part_name>-<rank>-<step>-<time>.h5", "%0.4f").GetFileName(), "test_model_part-0-3-0.0010.h5");
    KRATOS_CHECK_EQUAL(fic(model_part, "<model_part_name>-<rank>-<time>.h5", "0.4e").GetFileName(), "test_model_part-0-1.0000e-03.h5");

    KRATOS_CATCH("");
}

} // namespace Testing
} // namespace Kratos
