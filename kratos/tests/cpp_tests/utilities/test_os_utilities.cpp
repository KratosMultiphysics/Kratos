//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "utilities/os_utilities.h"
#include "../auxiliar_files_for_cpp_unnitest/aux_os_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(GetCurrentWorkingDir, KratosCoreFastSuite)
{
    // Using directly the utilities
    const std::string current_dir = OSUtilities::GetCurrentWorkingDir();

    // We call an auxilir utility, the folder should be the same
    const std::string external_current_dir = AuxiliarGetCurrentWorkingDir();

    // We check
    KRATOS_CHECK_STRING_EQUAL(current_dir, external_current_dir);
}

}   // namespace Testing
}  // namespace Kratos.
