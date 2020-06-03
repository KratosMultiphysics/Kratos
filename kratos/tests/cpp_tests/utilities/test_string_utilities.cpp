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
#include "utilities/string_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ConvertCammelCaseToSnakeCase, KratosCoreFastSuite)
{
    const std::string CammelCase = "TestInCammelCase";
    const std::string snake_case = StringUtilities::ConvertCammelCaseToSnakeCase(CammelCase);
    KRATOS_CHECK_STRING_EQUAL(snake_case, "test_in_cammel_case");
}

}   // namespace Testing
}  // namespace Kratos.
