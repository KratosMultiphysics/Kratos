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
#include "testing/testing.h"
#include "custom_workflows/dgeoflow.h"
#include "flow_stubs.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategy, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30.gid";
    auto projectFile = "ProjectParameters.json";

    auto execute = KratosExecute();
    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(3, 4, 0.1);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(&flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyProgress,
                                                                 &flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyCancel);

    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info, "PorousDomain.Left_head", call_back_functions);


    KRATOS_EXPECT_EQ(status, 0);
}

KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategyTextualProgressReport, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30.gid";
    auto projectFile = "ProjectParameters.json";

    auto execute = KratosExecute();

    bool firstMessageFound = false;
    bool finalMessageFound = false;
    int messageCount = 0;

    std::function<void(const char*)> reportTextualProgress = [&firstMessageFound, &finalMessageFound, &messageCount](const char* message)
    {
        messageCount++;

        if(strcmp(message, "Calculating head level 3m (1/12)") == 0) {
            firstMessageFound = true;
        }

        if(strcmp(message, "Calculating head level 3.8m (9/12)") == 0) {
            finalMessageFound = true;
        }
    };

    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(3, 4, 0.1);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(&flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyProgress,
                                                                 reportTextualProgress,
                                                                 &flow_stubs::emptyCancel);

    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info, "PorousDomain.Left_head", call_back_functions);

    KRATOS_EXPECT_EQ(status, 0);
    KRATOS_EXPECT_EQ(firstMessageFound, true);
    KRATOS_EXPECT_EQ(finalMessageFound, true);
    KRATOS_EXPECT_EQ(messageCount, 9);
}

KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategyProgressReport, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30.gid";
    auto projectFile = "ProjectParameters.json";

    auto execute = KratosExecute();

    bool startProgressFound = false;
    bool endProgressFound = false;
    int progressUpdates = 0;

    std::function<void(double)> reportProgress = [&startProgressFound, &endProgressFound, &progressUpdates](double progress)
    {
        progressUpdates++;

        if(progress == 0.0) {
            startProgressFound = true;
        }

        if(progress == 0.75) {
            endProgressFound = true;
        }
    };

    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(3, 4, 0.1);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(&flow_stubs::emptyLog,
                                                                 reportProgress,
                                                                 &flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyCancel);

    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info, "PorousDomain.Left_head", call_back_functions);


    KRATOS_EXPECT_EQ(status, 0);
    KRATOS_EXPECT_EQ(startProgressFound, true);
    KRATOS_EXPECT_EQ(endProgressFound, true);
    KRATOS_EXPECT_EQ(progressUpdates, 10);
}

}

