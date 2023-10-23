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

#include <string>
#include <filesystem>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_workflows/dgeosettlement.h"
#include "flow_stubs.h"
#include "custom_workflows/custom_workflow_factory.h"
#include "test_utilities.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SettlementWorkflow, KratosGeoMechanicsIntegrationSuite)
{
    const auto working_directory = std::filesystem::path{"./applications/GeoMechanicsApplication/tests/test_settlement_workflow"};

    auto settlement = CustomWorkflowFactory::CreateKratosGeoSettlement();
    for (int i = 0; i < 3; ++i) {
        auto projectFile = "ProjectParameters_stage"+ std::to_string(i + 1) + ".json";
        int status = settlement->RunStage(working_directory, projectFile,
                                          &flow_stubs::emptyLog, &flow_stubs::emptyProgress,
                                          &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

        std::string original = working_directory.generic_string() + std::string("/test_model_stage") + std::to_string(i+1) + ".post.orig.res";
        std::string result = working_directory.generic_string() + std::string("/test_model_stage") + std::to_string(i+1) + ".post.res";

        KRATOS_EXPECT_EQ(status, 0);
        KRATOS_EXPECT_TRUE(TestUtilities::CompareFiles(original, result))
    }
}

}