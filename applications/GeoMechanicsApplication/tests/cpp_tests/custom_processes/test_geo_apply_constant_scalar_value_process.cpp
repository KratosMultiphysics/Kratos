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

#include "custom_processes/geo_apply_constant_scalar_value_process.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include <algorithm>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoApplyConstantScalarValueProcess_FreesDoFAfterFinalize_ForDoubleVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      current_model;
    const auto nodal_variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X)};
    ModelPart& r_model_part =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(current_model, nodal_variables);

    Parameters parameters(R"(
      {
          "model_part_name" : "Main",
          "variable_name"   : "DISPLACEMENT_X",
          "is_fixed"        : true,
          "value"           : 1.0
      }  )");

    GeoApplyConstantScalarValueProcess process(r_model_part, parameters);
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(std::all_of(r_model_part.NodesBegin(), r_model_part.NodesEnd(), [](const auto& rNode) {
        return rNode.IsFixed(DISPLACEMENT_X) && rNode.FastGetSolutionStepValue(DISPLACEMENT_X) == 1.0;
    }))

    process.ExecuteFinalize();
    KRATOS_EXPECT_TRUE(std::all_of(r_model_part.NodesBegin(), r_model_part.NodesEnd(), [](const auto& rNode) {
        return !rNode.IsFixed(DISPLACEMENT_X) && rNode.FastGetSolutionStepValue(DISPLACEMENT_X) == 1.0;
    }))
}

KRATOS_TEST_CASE_IN_SUITE(GeoApplyConstantScalarValueProcess_FinalizeDoesNothing_ForIntVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main", 2);
    r_model_part.AddNodalSolutionStepVariable(TIME_STEPS);

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    Parameters parameters(R"(
      {
          "model_part_name" : "Main",
          "variable_name"   : "TIME_STEPS",
          "is_fixed"        : false,
          "value"           : 1
      }  )");

    GeoApplyConstantScalarValueProcess process(r_model_part, parameters);
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(std::all_of(r_model_part.NodesBegin(), r_model_part.NodesEnd(), [](const auto& rNode) {
        return rNode.FastGetSolutionStepValue(TIME_STEPS) == 1;
    }))

    process.ExecuteFinalize();
    KRATOS_EXPECT_TRUE(std::all_of(r_model_part.NodesBegin(), r_model_part.NodesEnd(), [](const auto& rNode) {
        return rNode.FastGetSolutionStepValue(TIME_STEPS) == 1;
    }))
}

} // namespace Kratos::Testing