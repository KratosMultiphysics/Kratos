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

#include <filesystem>
#include <string>

// Project includes
#include "containers/model.h"
#include "custom_workflows/custom_workflow_factory.h"
#include "custom_workflows/dgeosettlement.h"
#include "flow_stubs.h"
#include "test_utilities.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing
{

// KRATOS_TEST_CASE_IN_SUITE(SettlementWorkflow, KratosGeoMechanicsFastSuite)
[[maybe_unused]] void TestSettlement()
{
    const auto working_directory = std::filesystem::path{
        "./applications/GeoMechanicsApplication/tests/"
        "test_settlement_workflow"};

    auto settlement = CustomWorkflowFactory::CreateKratosGeoSettlement();
    for (int i = 0; i < 4; ++i)
    {
        const auto project_file =
            "ProjectParameters_stage" + std::to_string(i + 1) + ".json";
        const int status = settlement->RunStage(
            working_directory, project_file, &flow_stubs::emptyLog,
            &flow_stubs::emptyProgress, &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

        const std::string original_file =
            "test_model_stage" + std::to_string(i + 1) + ".post.orig.res";
        const std::string result_file =
            "test_model_stage" + std::to_string(i + 1) + ".post.res";

        KRATOS_EXPECT_EQ(status, 0);
        KRATOS_EXPECT_TRUE(TestUtilities::CompareFiles(
            working_directory / original_file, working_directory / result_file))
    }
}

} // namespace Kratos::Testing