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
//

// System includes
#include <fstream>
#include <string>

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

    KRATOS_CHECK(Kratos::filesystem::exists(file_name));

    Kratos::filesystem::remove(file_name);

    KRATOS_CHECK_IS_FALSE(Kratos::filesystem::exists(file_name));
}

} // namespace Testing
} // namespace Kratos
