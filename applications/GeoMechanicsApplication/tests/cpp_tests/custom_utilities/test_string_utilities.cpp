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

#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoStringUtilities_ConvertStringToLower, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto str = "ConvertMeToLowerCase"s;

    // Act
    const auto lower_case_string = GeoStringUtilities::ToLower(str);

    // Assert
    const auto expected_string = "convertmetolowercase"s;
    EXPECT_EQ(lower_case_string, expected_string);
}

KRATOS_TEST_CASE_IN_SUITE(GeoStringUtilities_JoiningAnEmptyVectorOfStringsYieldsAnEmptyString,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto no_strings = std::vector<std::string>{};

    // Act
    const auto joined_string = GeoStringUtilities::Join(no_strings, ", "s);

    // Assert
    const auto expected_string = ""s;
    EXPECT_EQ(joined_string, expected_string);
}

KRATOS_TEST_CASE_IN_SUITE(GeoStringUtilities_JoiningASingleStringYieldsItself, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto one_string = std::vector{"Foo"s};

    // Act
    const auto joined_string = GeoStringUtilities::Join(one_string, ", "s);

    // Assert
    const auto expected_string = "Foo"s;
    EXPECT_EQ(joined_string, expected_string);
}

KRATOS_TEST_CASE_IN_SUITE(GeoStringUtilities_JoiningMultipleStringsYieldsASingleStringWithEachElementSeparatedByTheGivenSeparator,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto strings = std::vector{"Foo"s, "Bar"s, "Baz"s};

    // Act
    const auto joined_string = GeoStringUtilities::Join(strings, ", "s);

    // Assert
    const auto expected_string = "Foo, Bar, Baz"s;
    EXPECT_EQ(joined_string, expected_string);
}

} // namespace Kratos::Testing
