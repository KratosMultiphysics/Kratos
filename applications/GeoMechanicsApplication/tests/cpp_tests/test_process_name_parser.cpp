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

#include "custom_utilities/process_name_parser.h"
#include "testing/testing.h"
#include "testing_utilities.h"

using namespace Kratos;

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(GetProcessNamesForMultipleConstraintProcesses, WorkInProgress)
{
    const std::string constraint_process_list = R"(
    {
      "constraints_process_list":
      [
          {
            "process_name": "ApplyVectorConstraintTableProcess"
          },
          {
            "process_name": "ApplyScalarConstraintTableProcess"
          }
      ]
    }
    )";

    ProcessNameParser parser;
    const auto actual_result = parser.GetProcessNames(constraint_process_list);

    const std::vector<std::string> expected_result = {"ApplyVectorConstraintTableProcess", "ApplyScalarConstraintTableProcess"};
    TestingUtilities::Expect_Equal(actual_result, expected_result);
}

KRATOS_TEST_CASE_IN_SUITE(GetProcessNamesForMultipleProcessesWithIdenticalNamesGivesDuplicates, WorkInProgress)
{
    const std::string constraint_process_list = R"(
    {
      "constraints_process_list":
      [
          {
            "process_name": "ApplyVectorConstraintTableProcess"
          },
          {
            "process_name": "ApplyVectorConstraintTableProcess"
          }
      ]
    }
    )";

    ProcessNameParser parser;
    const auto actual_result = parser.GetProcessNames(constraint_process_list);

    const std::vector<std::string> expected_result = {"ApplyVectorConstraintTableProcess", "ApplyVectorConstraintTableProcess"};
    TestingUtilities::Expect_Equal(actual_result, expected_result);
}


}