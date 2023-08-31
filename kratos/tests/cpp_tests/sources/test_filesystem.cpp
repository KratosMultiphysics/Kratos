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
#include <tuple>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_filesystem.h"
#include "utilities/parallel_utilities.h"
#include "testing/scoped_file.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ListDirectory, KratosCoreFastSuite)
{
    const std::filesystem::path dir_name("ListDirectoryTest");
    const std::filesystem::path dir_name_2("ListDirectoryTestSub");

    const std::filesystem::path raw_file_name_1("dummy_file_1.txt");
    const std::filesystem::path raw_file_name_2("dummy_file_2.txt");
    const std::filesystem::path sub_dir(dir_name / dir_name_2);
    const std::filesystem::path file_name_1(dir_name / raw_file_name_1);
    const std::filesystem::path file_name_2(dir_name / raw_file_name_2);

    KRATOS_EXPECT_FALSE(std::filesystem::exists(dir_name));
    KRATOS_EXPECT_FALSE(std::filesystem::exists(sub_dir));
    KRATOS_EXPECT_TRUE(std::filesystem::create_directories(sub_dir));

    std::ofstream output_file;
    output_file.open(file_name_1);
    output_file.close();
    output_file.open(file_name_2);
    output_file.close();

    KRATOS_EXPECT_TRUE(std::filesystem::exists(file_name_1));
    KRATOS_EXPECT_TRUE(std::filesystem::exists(file_name_2));
    KRATOS_EXPECT_TRUE(std::filesystem::exists(sub_dir));

    const std::vector<std::filesystem::path>& list_of_dirs = Kratos::FilesystemExtensions::ListDirectory(dir_name);
    const std::vector<std::filesystem::path> check_list{sub_dir, file_name_1, file_name_2};

    KRATOS_EXPECT_EQ(check_list.size(), list_of_dirs.size());

    for (const auto& r_dir : list_of_dirs) {
        bool found_check_dir = false;
        for (const auto& check_dir : check_list) {
            if (r_dir.parent_path() == check_dir.parent_path() &&
                r_dir.filename() == check_dir.filename()) {
                found_check_dir = true;
                break;
            }
        }
        KRATOS_EXPECT_TRUE(found_check_dir);
    }

    std::filesystem::remove_all(dir_name);

    KRATOS_EXPECT_FALSE(std::filesystem::exists(dir_name));
}

KRATOS_TEST_CASE_IN_SUITE(MPISafeCreateDirectories, KratosCoreFastSuite)
{
    auto create_dir_test_fct = [](const std::filesystem::path& rDirName){
        // make sure the dir does not exist already
        KRATOS_EXPECT_FALSE(std::filesystem::exists(rDirName));

        IndexPartition(100).for_each([&rDirName](std::size_t i){
            FilesystemExtensions::MPISafeCreateDirectories(rDirName);
        });

        KRATOS_EXPECT_TRUE(std::filesystem::exists(rDirName));

        // cleanup afterwards
        std::filesystem::remove_all(rDirName);
        KRATOS_EXPECT_FALSE(std::filesystem::exists(rDirName));
    };

    const std::filesystem::path base_dir_name("MyCustomDir2");
    const std::filesystem::path sub_dir_name("TheSubDir");
    const std::filesystem::path sub_sub_dir_name("TheSubSubDir");

    const std::filesystem::path full_dir_name_1 = base_dir_name / sub_dir_name;
    const std::filesystem::path full_dir_name_2 = base_dir_name / sub_dir_name / sub_sub_dir_name;

    create_dir_test_fct(base_dir_name);
    create_dir_test_fct(full_dir_name_1);
    create_dir_test_fct(full_dir_name_2);

    // final cleanup after test
    std::filesystem::remove_all(base_dir_name);
    KRATOS_EXPECT_FALSE(std::filesystem::exists(base_dir_name));
}

KRATOS_TEST_CASE_IN_SUITE(ResolveSymlinksToFile, KratosCoreFastSuite)
{
    using Path = std::filesystem::path;

    const ScopedDirectory working_dir("test_resolve_symlinks");
    const Path file_path = working_dir / Path("file.test");

    std::vector<ScopedSymlink> symlinks;
    symlinks.reserve(3);

    // Create symlinks
    for (unsigned short level=0; level<=2; ++level) {
            symlinks.emplace_back(working_dir / Path("symlink_" + std::to_string(level) + ".test"), level ? symlinks.back() : file_path);
    }

    {
        // Level 0-2 indirections to a non-existent path
        for (const auto& r_symlink : symlinks) {
            KRATOS_EXPECT_EQ(FilesystemExtensions::ResolveSymlinks(r_symlink), file_path);
        }
    }

    {
        // Level 0-2 indirections to an existing path
        const ScopedFile file(file_path);
        for (const auto& r_symlink : symlinks) {
            KRATOS_EXPECT_EQ(FilesystemExtensions::ResolveSymlinks(r_symlink), file_path);
        }

        // Input is a file
        KRATOS_EXPECT_EQ(FilesystemExtensions::ResolveSymlinks(file), file_path);
    }

    {
        // Cyclic indirections 2-4 cycles
        for (const auto& r_symlink : symlinks) {
            ScopedSymlink loopback(file_path, r_symlink);
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(std::ignore=FilesystemExtensions::ResolveSymlinks(r_symlink), "cyclic");
        }

        // 1-cycle
        ScopedSymlink loopback(file_path, file_path);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(std::ignore=FilesystemExtensions::ResolveSymlinks(loopback), "cyclic");
    }
}

KRATOS_TEST_CASE_IN_SUITE(ResolveSymlinksToDirectory, KratosCoreFastSuite)
{
    using Path = std::filesystem::path;

    const ScopedDirectory working_dir("test_resolve_symlinks");
    const Path directory_path = working_dir / Path("directory");

    std::vector<ScopedSymlink> symlinks;
    symlinks.reserve(3);

    // Create symlinks
    for (unsigned short level=0; level<=2; ++level) {
        symlinks.emplace_back(working_dir / Path("symlink_" + std::to_string(level) + ".test"), level ? symlinks.back() : directory_path);
    }

    {
        // Level 0-2 indirections to a non-existent path
        for (const auto& r_symlink : symlinks) {
            KRATOS_EXPECT_EQ(FilesystemExtensions::ResolveSymlinks(r_symlink), directory_path);
        }
    }

    {
        // Level 0-2 indirections to an existing directory
        const ScopedDirectory directory(directory_path);
        for (const auto& r_symlink : symlinks) {
            KRATOS_EXPECT_EQ(FilesystemExtensions::ResolveSymlinks(r_symlink), directory_path);
        }

        // Input is a directory
        KRATOS_EXPECT_EQ(FilesystemExtensions::ResolveSymlinks(directory), directory_path);
    }

    {
        // Cyclic indirections 2-4 cycles
        for (const auto& r_symlink : symlinks) {
            ScopedSymlink loopback(directory_path, r_symlink);
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(std::ignore=FilesystemExtensions::ResolveSymlinks(r_symlink), "cyclic");
        }

        // 1-cycle
        ScopedSymlink loopback(directory_path, directory_path);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(std::ignore=FilesystemExtensions::ResolveSymlinks(loopback), "cyclic");
    }
}

} // namespace Testing
} // namespace Kratos
