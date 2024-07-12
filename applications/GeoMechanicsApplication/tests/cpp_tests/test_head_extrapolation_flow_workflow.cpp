// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#include <iostream>
#include <string>

// Project includes
#include "custom_workflows/dgeoflow.h"
#include "flow_stubs.h"
#include "geo_mechanics_fast_suite.h"
#include "test_utilities.h"

namespace
{

using namespace Kratos;
using namespace Kratos::Testing;

const auto workingDirectory = std::filesystem::path{"."} / "applications" / "GeoMechanicsApplication" /
                              "tests" / "test_head_extrapolation_custom_workflow_flow";

int RunTestCase(int test_case_number)
{
    const auto projectFile = std::string{"ProjectParameters_"} + std::to_string(test_case_number) + ".json";

    auto execute = KratosExecute();

    const KratosExecute::CriticalHeadInfo  critical_head_info(0, 0, 0);
    const KratosExecute::CallBackFunctions call_back_functions(
        &flow_stubs::emptyLog, &flow_stubs::emptyProgress, &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

    return execute.ExecuteFlowAnalysis(workingDirectory.generic_string(), projectFile,
                                       critical_head_info, "", call_back_functions);
}

void CompareResults(int test_case_number)
{
    const auto original =
        (workingDirectory / ("test_head_extrapolate_" + std::to_string(test_case_number) + ".orig.res"))
            .generic_string();
    const auto result =
        (workingDirectory / ("test_head_extrapolate_" + std::to_string(test_case_number) + ".post.res"))
            .generic_string();

    KRATOS_EXPECT_TRUE(TestUtilities::CompareFiles(original, result))
}
}


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_1, KratosGeoMechanicsIntegrationSuite)
{
    constexpr auto test_case_number = 1;
    KRATOS_EXPECT_EQ(RunTestCase(test_case_number), 0);
    CompareResults(test_case_number);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_2, KratosGeoMechanicsIntegrationSuite)
{
    constexpr auto test_case_number = 2;
    KRATOS_EXPECT_EQ(RunTestCase(test_case_number), 0);
    CompareResults(test_case_number);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_3, KratosGeoMechanicsIntegrationSuite)
{
    constexpr auto test_case_number = 3;
    KRATOS_EXPECT_EQ(RunTestCase(test_case_number), 0);
    CompareResults(test_case_number);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_4, KratosGeoMechanicsIntegrationSuite)
{
    constexpr auto test_case_number = 4;
    KRATOS_EXPECT_EQ(RunTestCase(test_case_number), 0);
    CompareResults(test_case_number);
}

} // namespace Kratos::Testing
