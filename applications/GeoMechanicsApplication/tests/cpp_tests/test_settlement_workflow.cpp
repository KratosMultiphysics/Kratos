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
#include <iostream>

#include <fstream>
#include <iterator>
#include <filesystem>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_workflows/dgeosettlement.h"
#include "flow_stubs.h"
#include "custom_workflows/custom_workflow_factory.h"
namespace
{


bool compareFiles(const std::string &p1, const std::string &p2)
{
    std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

    if (f1.fail() || f2.fail()) {
        return false; // file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; // size mismatch
    }

    // seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(f2.rdbuf()));
}
}

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SettlementWorkflow, KratosGeoMechanicsFastSuite)
{
    const auto working_directory = std::filesystem::path{"./applications/GeoMechanicsApplication/tests/test_settlement_workflow"};

    auto settlement = CustomWorkflowFactory::CreateKratosGeoSettlement();
    std::ofstream stream;
    stream.open("./applications/GeoMechanicsApplication/tests/test_settlement_workflow/test_output.txt", std::ios_base::out);
    auto log_callback = [&stream](const char* output){stream << output;};
    for (int i = 0; i < 1; ++i) {
        auto projectFile = "ProjectParameters_stage"+ std::to_string(i + 1) + ".json";
        int status = settlement->RunStage(working_directory, projectFile,
                                          log_callback, &flow_stubs::emptyProgress,
                                          &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

        KRATOS_EXPECT_EQ(status, 0);
    }


    // output_files
//    std::string original = (std::string) workingDirectory + "/test_head_extrapolate_1.orig.res";
//    std::string result = (std::string) workingDirectory + "/test_head_extrapolate_1.post.res";

//    KRATOS_EXPECT_TRUE(compareFiles(original, result))
}

}