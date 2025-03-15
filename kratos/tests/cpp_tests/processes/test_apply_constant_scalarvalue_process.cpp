//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Richard Faasse
//

#include "processes/apply_constant_scalarvalue_process.h"
#include "tests/test_utilities/cpp_tests_utilities.h"
#include "tests/test_utilities/test_suite.h"

#include <algorithm>

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantScalarValueProcess_FreesDoFAfterFinalize_ForDoubleVariable,
    KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main",2);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
    }

    Parameters parameters( R"(
      {
          "model_part_name" : "Main",
          "mesh_id"         : 0,
          "variable_name"   : "DISPLACEMENT_X",
          "is_fixed"        : true,
          "value"           : 1.0
      }  )" );

    ApplyConstantScalarValueProcess process(r_model_part, parameters);
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(std::all_of(
        r_model_part.NodesBegin(), r_model_part.NodesEnd(),
        [](const auto& rNode) {
          return rNode.IsFixed(DISPLACEMENT_X) &&
                 rNode.FastGetSolutionStepValue(DISPLACEMENT_X) == 1.0;
        }))

    process.ExecuteFinalize();
    KRATOS_EXPECT_TRUE(std::all_of(
        r_model_part.NodesBegin(), r_model_part.NodesEnd(),
        [](const auto& rNode) {
          return !rNode.IsFixed(DISPLACEMENT_X) &&
                 rNode.FastGetSolutionStepValue(DISPLACEMENT_X) == 1.0;
        }))
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantScalarValueProcess_FinalizeDoesNothing_ForIntVariable,
    KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main",2);
    r_model_part.AddNodalSolutionStepVariable(TIME_STEPS);

    CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

    Parameters parameters( R"(
      {
          "model_part_name" : "Main",
          "mesh_id"         : 0,
          "variable_name"   : "TIME_STEPS",
          "is_fixed"        : false,
          "value"           : 1
      }  )" );

    ApplyConstantScalarValueProcess process(r_model_part, parameters);
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(std::all_of(
        r_model_part.NodesBegin(), r_model_part.NodesEnd(),
        [](const auto& rNode) {
          return rNode.FastGetSolutionStepValue(TIME_STEPS) == 1;
        }))

    process.ExecuteFinalize();
    KRATOS_EXPECT_TRUE(std::all_of(
        r_model_part.NodesBegin(), r_model_part.NodesEnd(),
        [](const auto& rNode) {
          return rNode.FastGetSolutionStepValue(TIME_STEPS) == 1;
        }))
}


}