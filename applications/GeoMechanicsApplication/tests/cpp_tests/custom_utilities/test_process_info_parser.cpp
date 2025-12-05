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
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

using namespace Kratos;

namespace Kratos::Testing
{

const std::string parameterString = R"(
{
    "key" : "value"
})";

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetProcessList_ReturnsExpectedProcesses_ForAllListTypes)
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
        "auxiliary_process_list":
        [
            {
                "process_name": "AuxiliaryProcess1",
                "Parameters":)" + parameterString +
                                     R"(
            },
            {
                "process_name": "AuxiliaryProcess2",
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
        ProcessParameters{"AuxiliaryProcess1", Parameters{parameterString}},
        ProcessParameters{"AuxiliaryProcess2", Parameters{parameterString}}};

    EXPECT_EQ(expected_result, actual_result);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetProcessList_GivesDuplicates_ForProcessesWithIdenticalNames)
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

    EXPECT_EQ(expected_result, actual_result);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetProcessList_Throws_WhenProcessDoesNotHaveParameters)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetProcessList_ReturnsEmptyList_WhenNoProcessesAreDefined)
{
    JsonProcessInfoParser parser;
    KRATOS_EXPECT_TRUE(parser.GetProcessList(Parameters{}).empty())
}

} // namespace Kratos::Testing
