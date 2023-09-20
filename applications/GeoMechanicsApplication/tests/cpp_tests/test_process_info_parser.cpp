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

#include "custom_utilities/process_info_parser.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing {

const std::string parameterString = R"(
{
    "key" : "value"
})";

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_ReturnsExpectedProcesses_ForAllListTypes, WorkInProgress)
{
    // Given
    const std::string process_list = R"(
    {
        "constraints_process_list":
        [
            {
                "process_name": "ConstraintProcess1",
                "Parameters":)" + parameterString + R"(
            },
            {
                "process_name": "ConstraintProcess2",
                "Parameters":)" + parameterString + R"(
            }
        ],
        "loads_process_list":
        [
            {
                "process_name": "LoadProcess1",
                "Parameters":)" + parameterString + R"(
            },
            {
                "process_name": "LoadProcess2",
                "Parameters":)" + parameterString + R"(
            }
        ],
        "auxiliar_process_list":
        [
            {
                "process_name": "AuxiliarProcess1",
                "Parameters":)" + parameterString + R"(
            },
            {
                "process_name": "AuxiliarProcess2",
                "Parameters":)" + parameterString + R"(
            }
        ]
    }
    )";

    // When
    ProcessInfoParser parser;
    const auto actual_result = parser.GetProcessList(process_list);

    // Then
    const std::vector<ProcessInfo> expected_result = {ProcessInfo{Parameters{parameterString}, "ConstraintProcess1"},
                                                      ProcessInfo{Parameters{parameterString}, "ConstraintProcess2"},
                                                      ProcessInfo{Parameters{parameterString}, "LoadProcess1"},
                                                      ProcessInfo{Parameters{parameterString}, "LoadProcess2"},
                                                      ProcessInfo{Parameters{parameterString}, "AuxiliarProcess1"},
                                                      ProcessInfo{Parameters{parameterString}, "AuxiliarProcess2"}};

    KRATOS_EXPECT_EQ(expected_result, actual_result);
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_GivesDuplicates_ForProcessesWithIdenticalNames, WorkInProgress)
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
    ProcessInfoParser parser;
    const auto actual_result = parser.GetProcessList(process_list_with_duplicate_names);

    // Then
    const std::vector<ProcessInfo> expected_result = {ProcessInfo{Parameters{parameterString}, "ApplyVectorConstraintTableProcess"},
                                                      ProcessInfo{Parameters{parameterString}, "ApplyVectorConstraintTableProcess"}};

    KRATOS_EXPECT_EQ(expected_result, actual_result);
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_Throws_WhenProcessDoesNotHaveParameters, WorkInProgress)
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
    ProcessInfoParser parser;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parser.GetProcessList(process_without_parameters), "Getting a value that does not exist. entry string : Parameters")
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_ReturnsEmptyString_WhenNoProcessesAreDefined, WorkInProgress)
{
    ProcessInfoParser parser;
    KRATOS_EXPECT_EQ(parser.GetProcessList(Parameters{}).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessList_ReturnsCorrectList_ForOutputProcess, WorkInProgress)
{
    const std::string output_process_list = R"(
    {
        "gid_output":
        [{
            "python_module": "gid_output_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "GiDOutputProcess",
            "Parameters":)" + parameterString + R"(
        }]
    }
    )";
    ProcessInfoParser parser;
    const auto actual_result = parser.GetProcessList(output_process_list);

    const std::vector<ProcessInfo> expected_result = {ProcessInfo{Parameters{parameterString}, "GiDOutputProcess"}};

    KRATOS_EXPECT_EQ(expected_result, actual_result);
}



}