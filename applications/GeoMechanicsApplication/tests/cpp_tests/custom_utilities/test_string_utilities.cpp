// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#include "custom_utilities/string_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoStringUtilities_ConvertStringToLower, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto str = "ConvertMeToLowerCase";

    // Act
    const auto lower_case_string = GeoStringUtilities::ToLower(str);

    // Assert
    const auto expected_string = "convertmetolowercase";
    KRATOS_CHECK_STRING_EQUAL(lower_case_string, expected_string);
}

} // namespace Kratos::Testing
