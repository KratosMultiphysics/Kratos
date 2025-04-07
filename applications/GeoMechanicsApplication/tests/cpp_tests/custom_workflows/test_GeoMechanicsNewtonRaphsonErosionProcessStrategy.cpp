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

// System includes
#include <map>

/* Project includes */
#include "custom_workflows/dgeoflow.h"
#include "tests/cpp_tests/flow_stubs.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using namespace Kratos;

class ParametrizedDGeoFlowTests : public ::testing::TestWithParam<std::string>
{
};

TEST_P(ParametrizedDGeoFlowTests, ErosionProcessStrategy)
{
    // Arrange
    auto workingDirectory =
        "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/" + GetParam();
    auto projectFile = "ProjectParameters.json";

    auto execute = KratosExecute();

    std::vector<double> progress_report_values;
    auto                reportProgress = [&progress_report_values](double progress) {
        progress_report_values.push_back(progress);
    };

    std::vector<std::string> reported_textual_progress;
    auto reportTextualProgress = [&reported_textual_progress](const char* message) {
        reported_textual_progress.emplace_back(message);
    };

    const KratosExecute::CriticalHeadInfo  critical_head_info(3, 4, 0.1);
    const KratosExecute::CallBackFunctions call_back_functions(
        &flow_stubs::emptyLog, reportProgress, reportTextualProgress, &flow_stubs::emptyCancel);

    // Act
    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info,
                                                   "PorousDomain.Left_head", call_back_functions);

    // Assert
    KRATOS_EXPECT_EQ(status, 0);
    KRATOS_EXPECT_EQ(reported_textual_progress.front(), "Calculating head level 3m (1/12)");
    KRATOS_EXPECT_EQ(reported_textual_progress.back(), "Calculating head level 3.8m (9/12)");
    KRATOS_EXPECT_EQ(reported_textual_progress.size(), 9);
    KRATOS_EXPECT_EQ(progress_report_values.front(), 0.0);
    KRATOS_EXPECT_EQ(progress_report_values.back(), 0.75);
    KRATOS_EXPECT_EQ(progress_report_values.size(), 10);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsIntegrationSuite,
                         ParametrizedDGeoFlowTests,
                         ::testing::Values("HeightAquiferD10L30.gid", "HeightAquiferD10L30line"));

} // namespace Kratos::Testing
