//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Vicente Mataix Ferrándiz
//

// System includes
#include <fstream>
#include <string>
#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_filesystem.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(FilesystemExists, KratosCoreFastSuite)
{
    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists("StupidNameShouldNotExist"));

    const std::string file_name("dummy_file.txt");
    std::ofstream output_file;
    output_file.open(file_name);
    output_file.close();

    KRATOS_CHECK(Kratos::filesystem::exists(file_name));

    Kratos::filesystem::remove(file_name);

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
}

KRATOS_TEST_CASE_IN_SUITE(FilesystemCurrentPath, KratosCoreFastSuite)
{
    // This test check that the current path method works
    const std::string current_path = (Kratos::filesystem::current_path()).string();
    
    char buff[FILENAME_MAX];
    GetCurrentDir( buff, FILENAME_MAX );
    const std::string current_working_dir(buff);
    
    KRATOS_CHECK_STRING_EQUAL(current_working_dir, current_path);
}

} // namespace Testing
} // namespace Kratos
