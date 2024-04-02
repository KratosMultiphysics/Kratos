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

#include <string>
#include <iostream>

// Project includes
#include "testing/testing.h"
#include "custom_workflows/dgeoflow.h"
#include "flow_stubs.h"
#include "test_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_1, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_1.json";

    auto execute = Kratos::KratosExecute();

    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(0, 0, 0);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(&flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyProgress,
                                                                 &flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyCancel);

    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info, "", call_back_functions);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_1.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_1.post.res";

    KRATOS_EXPECT_TRUE(TestUtilities::CompareFiles(original, result))
}

 KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_2, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_2.json";

    auto execute = Kratos::KratosExecute();

    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(0, 0, 0);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(&flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyProgress,
                                                                 &flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyCancel);

    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info, "", call_back_functions);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_2.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_2.post.res";

    KRATOS_EXPECT_TRUE(TestUtilities::CompareFiles(original, result))
}

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_3, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_3.json";

    auto execute = Kratos::KratosExecute();

    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(0, 0, 0);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(&flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyProgress,
                                                                 &flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyCancel);

    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info, "", call_back_functions);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_3.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_3.post.res";

    KRATOS_EXPECT_TRUE(TestUtilities::CompareFiles(original, result))
}

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_4, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_4.json";

    auto execute = Kratos::KratosExecute();

    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(0, 0, 0);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(&flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyProgress,
                                                                 &flow_stubs::emptyLog,
                                                                 &flow_stubs::emptyCancel);

    const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info, "", call_back_functions);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_4.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_4.post.res";

    KRATOS_EXPECT_TRUE(TestUtilities::CompareFiles(original, result))
}

}
