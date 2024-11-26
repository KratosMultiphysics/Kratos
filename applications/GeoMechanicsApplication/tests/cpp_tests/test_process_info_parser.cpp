// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_utilities/json_process_info_parser.h"
#include "geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

const std::string parameterString = R"(
{
    "key" : "value"
})";

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_ReturnsExpectedProcesses_ForAllListTypes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Given
    const std::string process_list = R"(
    {
        "constraints_process_list":
        [
            {
                "process_name": "ConstraintProcess1",
                "Parameters":)" + parameterString +
                                     R"(
            },
            {
                "process_name": "ConstraintProcess2",
                "Parameters":)" + parameterString +
                                     R"(
            }
        ],
        "loads_process_list":
        [
            {
                "process_name": "LoadProcess1",
                "Parameters":)" + parameterString +
                                     R"(
            },
            {
                "process_name": "LoadProcess2",
                "Parameters":)" + parameterString +
                                     R"(
            }
        ],
        "auxiliar_process_list":
        [
            {
                "process_name": "AuxiliarProcess1",
                "Parameters":)" + parameterString +
                                     R"(
            },
            {
                "process_name": "AuxiliarProcess2",
                "Parameters":)" + parameterString +
                                     R"(
            }
        ]
    }
    )";

    // When
    JsonProcessInfoParser parser;
    const auto            actual_result = parser.GetProcessList(process_list);

    // Then
    const std::vector<ProcessParameters> expected_result = {
        ProcessParameters{"ConstraintProcess1", Parameters{parameterString}},
        ProcessParameters{"ConstraintProcess2", Parameters{parameterString}},
        ProcessParameters{"LoadProcess1", Parameters{parameterString}},
        ProcessParameters{"LoadProcess2", Parameters{parameterString}},
        ProcessParameters{"AuxiliarProcess1", Parameters{parameterString}},
        ProcessParameters{"AuxiliarProcess2", Parameters{parameterString}}};

    KRATOS_EXPECT_EQ(expected_result, actual_result);
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_GivesDuplicates_ForProcessesWithIdenticalNames,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Given
    const std::string process_list_with_duplicate_names = R"(
    {
      "constraints_process_list":
      [
          {
            "process_name": "ApplyVectorConstraintTableProcess",
            "Parameters":)" + parameterString + R"(
          },
          {
            "process_name": "ApplyVectorConstraintTableProcess",
            "Parameters":)" + parameterString + R"(
          }
      ]
    }
    )";

    // When
    JsonProcessInfoParser parser;
    const auto            actual_result = parser.GetProcessList(process_list_with_duplicate_names);

    // Then
    const std::vector<ProcessParameters> expected_result = {
        ProcessParameters{"ApplyVectorConstraintTableProcess", Parameters{parameterString}},
        ProcessParameters{"ApplyVectorConstraintTableProcess", Parameters{parameterString}}};

    KRATOS_EXPECT_EQ(expected_result, actual_result);
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_Throws_WhenProcessDoesNotHaveParameters, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Given
    const std::string process_without_parameters = R"(
    {
      "constraints_process_list":
      [
          {
            "process_name": "ApplyVectorConstraintTableProcess"
          }
      ]
    }
    )";

    // When
    JsonProcessInfoParser parser;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        parser.GetProcessList(process_without_parameters),
        "Getting a value that does not exist. entry string : Parameters")
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_ReturnsEmptyList_WhenNoProcessesAreDefined, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    JsonProcessInfoParser parser;
    KRATOS_EXPECT_TRUE(parser.GetProcessList(Parameters{}).empty())
}

} // namespace Kratos::Testing