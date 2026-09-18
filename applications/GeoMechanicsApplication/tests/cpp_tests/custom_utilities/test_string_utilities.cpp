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

struct JoinStringsTestData {
    std::vector<std::string> InputStrings;
    std::string              Separator;
    std::string              ExpectedResult;
};

class ParameterizedJoinStringsTest : public ::testing::TestWithParam<JoinStringsTestData>
{
};

TEST_P(ParameterizedJoinStringsTest, JoiningMultipleStringsYieldsASingleStringWithEachElementSeparatedByTheGivenSeparator)
{
    // Arrange
    const auto& [input_strings, separator, expected_result] = GetParam();

    // Act
    const auto joined_string = GeoStringUtilities::Join(input_strings, separator);

    // Assert
    EXPECT_EQ(joined_string, expected_result);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParameterizedJoinStringsTest,
                         testing::Values(JoinStringsTestData{{}, ", "s, {}},
                                         JoinStringsTestData{{"Foo"s}, ", "s, "Foo"s},
                                         JoinStringsTestData{{"Foo"s, "Bar"s, "Baz"s}, ", "s, "Foo, Bar, Baz"s},
                                         JoinStringsTestData{{"Foo"s, "Bar"s, "Baz"s}, {}, "FooBarBaz"s}));

} // namespace Kratos::Testing
