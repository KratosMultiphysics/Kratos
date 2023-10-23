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

#pragma once
#include <string>
#include <iostream>

#include <fstream>
#include <iterator>
#include <algorithm>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_workflows/dgeoflow.h"
#include "flow_stubs.h"

bool compareFiles(const std::string &p1, const std::string &p2)
{
    std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

    if (f1.fail() || f2.fail())
    {
        return false; // file problem
    }

    if (f1.tellg() != f2.tellg())
    {
        return false; // size mismatch
    }

    // seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(f2.rdbuf()));
}

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_1, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_1.json";

    auto execute = Kratos::KratosExecute();
    int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile,
                                             0, 0, 0,
                                             "", &flow_stubs::emptyLog, &flow_stubs::emptyProgress,
                                             &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_1.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_1.post.res";

    KRATOS_EXPECT_TRUE(compareFiles(original, result))
}

 KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_2, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_2.json";

    auto execute = Kratos::KratosExecute();
    int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile,
                                             0, 0, 0,
                                             "", &flow_stubs::emptyLog, &flow_stubs::emptyProgress,
                                             &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_2.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_2.post.res";

    KRATOS_EXPECT_TRUE(compareFiles(original, result))
}

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_3, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_3.json";

    auto execute = Kratos::KratosExecute();
    int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile,
                                             0, 0, 0,
                                             "", &flow_stubs::emptyLog, &flow_stubs::emptyProgress,
                                             &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_3.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_3.post.res";

    KRATOS_EXPECT_TRUE(compareFiles(original, result))
}

KRATOS_TEST_CASE_IN_SUITE(CalculateExtrapolatedHeadFlow_4, KratosGeoMechanicsIntegrationSuite)
{
    auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_head_extrapolation_custom_workflow_flow";
    auto projectFile = "ProjectParameters_4.json";

    auto execute = Kratos::KratosExecute();
    int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile,
                                             0, 0, 0,
                                             "", &flow_stubs::emptyLog, &flow_stubs::emptyProgress,
                                             &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

    KRATOS_EXPECT_EQ(status, 0);

    // output_files
    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_4.orig.res";
    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_4.post.res";

    KRATOS_EXPECT_TRUE(compareFiles(original, result))
}

}
