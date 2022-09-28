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
#include "testing/testing.h"
#include "includes/kratos_filesystem.h"
#include "utilities/parallel_utilities.h"

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

KRATOS_TEST_CASE_IN_SUITE(FileSystemParentPathFilename, KratosCoreFastSuite)
{
    std::vector<std::string> paths_1 {"sl", "", "uom", "dssc"};
    const std::string& path = Kratos::FilesystemExtensions::JoinPaths(paths_1);
    const auto& parent_path = Kratos::filesystem::parent_path(path);
    const auto& parent_parent_path = Kratos::filesystem::parent_path(parent_path);
    KRATOS_CHECK_STRING_EQUAL(parent_parent_path, "sl");
    KRATOS_CHECK_STRING_EQUAL(Kratos::filesystem::filename(parent_path), "uom");
    KRATOS_CHECK_STRING_EQUAL(Kratos::filesystem::filename(path), "dssc");
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

KRATOS_TEST_CASE_IN_SUITE(ListDirectory, KratosCoreFastSuite)
{
    const std::string dir_name("ListDirectoryTest");
    const std::string dir_name_2("ListDirectoryTestSub");

    const std::string raw_file_name_1("dummy_file_1.txt");
    const std::string raw_file_name_2("dummy_file_2.txt");
    const std::string sub_dir(Kratos::FilesystemExtensions::JoinPaths({dir_name, dir_name_2}));
    const std::string file_name_1(Kratos::FilesystemExtensions::JoinPaths({dir_name, raw_file_name_1}));
    const std::string file_name_2(Kratos::FilesystemExtensions::JoinPaths({dir_name, raw_file_name_2}));

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name));
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(sub_dir));
    KRATOS_CHECK(Kratos::filesystem::create_directories(sub_dir));

    std::ofstream output_file;
    output_file.open(file_name_1);
    output_file.close();
    output_file.open(file_name_2);
    output_file.close();

    KRATOS_CHECK(Kratos::filesystem::exists(file_name_1));
    KRATOS_CHECK(Kratos::filesystem::exists(file_name_2));
    KRATOS_CHECK(Kratos::filesystem::exists(sub_dir));

    const auto& list_of_dirs = Kratos::FilesystemExtensions::ListDirectory(dir_name);
    const std::vector<std::string> check_list = {sub_dir, file_name_1, file_name_2};
    KRATOS_CHECK_EQUAL(check_list.size(), list_of_dirs.size());

    for (const auto& r_dir : list_of_dirs) {
        bool found_check_dir = false;
        for (const auto& check_dir : check_list) {
            if (filesystem::parent_path(r_dir) == filesystem::parent_path(check_dir) &&
                filesystem::filename(r_dir) == filesystem::filename(check_dir)) {
                found_check_dir = true;
                break;
            }
        }
        KRATOS_CHECK(found_check_dir);
    }

    Kratos::filesystem::remove_all(dir_name);

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(dir_name));
}

KRATOS_TEST_CASE_IN_SUITE(MPISafeCreateDirectories, KratosCoreFastSuite)
{
    auto create_dir_test_fct = [](const std::string& rDirName){
        // make sure the dir does not exist already
        KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(rDirName));

        IndexPartition<> index_part(100);
        index_part.for_each([&rDirName](std::size_t i){
            FilesystemExtensions::MPISafeCreateDirectories(rDirName);
        });

        KRATOS_CHECK(Kratos::filesystem::exists(rDirName));

        // cleanup afterwards
        Kratos::filesystem::remove_all(rDirName);
        KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(rDirName));
    };

    const std::string base_dir_name("MyCustomDir2");
    const std::string sub_dir_name("TheSubDir");
    const std::string sub_sub_dir_name("TheSubSubDir");

    const std::string full_dir_name_1 = Kratos::FilesystemExtensions::JoinPaths({base_dir_name, sub_dir_name});
    const std::string full_dir_name_2 = Kratos::FilesystemExtensions::JoinPaths({base_dir_name, sub_dir_name, sub_sub_dir_name});

    create_dir_test_fct(base_dir_name);
    create_dir_test_fct(full_dir_name_1);
    create_dir_test_fct(full_dir_name_2);

    // final cleanup after test
    Kratos::filesystem::remove_all(base_dir_name);
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(base_dir_name));
}

} // namespace Testing
} // namespace Kratos
