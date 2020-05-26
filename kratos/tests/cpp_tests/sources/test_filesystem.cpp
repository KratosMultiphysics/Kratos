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

} // namespace Testing
} // namespace Kratos
