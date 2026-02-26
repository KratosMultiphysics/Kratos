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
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <algorithm>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoApplyConstantScalarValueProcess_FreesDoFAfterFinalize_ForDoubleVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      current_model;
    const auto nodal_variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X)};
    auto&      r_model_part =
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
        return rNode.IsFixed(DISPLACEMENT_X) && rNode.FastGetSolutionStepValue(DISPLACEMENT_X) == 0.0;
    }))

    process.ExecuteInitializeSolutionStep();
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
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("Main", 2);
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
        return rNode.FastGetSolutionStepValue(TIME_STEPS) == 0;
    }))

    process.ExecuteInitializeSolutionStep();
    KRATOS_EXPECT_TRUE(std::all_of(r_model_part.NodesBegin(), r_model_part.NodesEnd(), [](const auto& rNode) {
        return rNode.FastGetSolutionStepValue(TIME_STEPS) == 1;
    }))

    process.ExecuteFinalize();
    KRATOS_EXPECT_TRUE(std::all_of(r_model_part.NodesBegin(), r_model_part.NodesEnd(), [](const auto& rNode) {
        return rNode.FastGetSolutionStepValue(TIME_STEPS) == 1;
    }))
}

KRATOS_TEST_CASE_IN_SUITE(GeoApplyConstantScalarValueProcess_Throws_WhenTryingToFixIntVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(TIME_STEPS);

    Parameters parameters(R"(
      {
          "model_part_name" : "Main",
          "variable_name"   : "TIME_STEPS",
          "is_fixed"        : true,
          "value"           : 1
      }  )");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoApplyConstantScalarValueProcess process(r_model_part, parameters),
        "It is not possible to fix the variable 'TIME_STEPS' which is not of type Variable<double>")
}

KRATOS_TEST_CASE_IN_SUITE(GeoApplyConstantScalarValueProcess_ThrowsWhenNodalVariableNotInModelPart,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("Main");

    Parameters parameters(R"(
      {
          "model_part_name" : "Main",
          "variable_name"   : "DISPLACEMENT_X",
          "is_fixed"        : true,
          "value"           : 1.0
      }  )");

    const std::string expected_error_message = "Trying to fix a variable that is not in ModelPart "
                                               "'Main' - variable name is DISPLACEMENT_X";
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoApplyConstantScalarValueProcess process(r_model_part, parameters), expected_error_message)
}

KRATOS_TEST_CASE_IN_SUITE(GeoApplyConstantScalarValueProcess_ThrowsWhenVariableNameIsMissing,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("Main");

    Parameters parameters_without_variable_name(R"(
      {
          "model_part_name" : "Main",
          "is_fixed"        : true,
          "value"           : 1.0
      }  )");

    const std::string expected_error_message = "Missing 'variable_name' parameter in the "
                                               "parameters of 'GeoApplyConstantScalarValueProcess'";
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoApplyConstantScalarValueProcess process(r_model_part, parameters_without_variable_name),
        expected_error_message);
}

KRATOS_TEST_CASE_IN_SUITE(GeoApplyConstantScalarValueProcess_ThrowsWhenValueIsMissing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("Main");

    Parameters parameters_without_value(R"(
      {
          "model_part_name" : "Main",
          "variable_name"   : "DISPLACEMENT_X",
          "is_fixed"        : true
      }  )");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoApplyConstantScalarValueProcess process(r_model_part, parameters_without_value),
        "Missing 'value' parameter in the parameters of 'GeoApplyConstantScalarValueProcess'");
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoGeoApplyConstantScalarValueProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model      model;
    const auto nodal_variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X)};
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model, nodal_variables);
    Parameters                               parameters(R"(
      {
          "model_part_name" : "Main",
          "variable_name"   : "DISPLACEMENT_X",
          "is_fixed"        : true,
          "value"           : 1.0
      }  )");
    const GeoApplyConstantScalarValueProcess process(r_model_part, parameters);
    // Act & assert
    KRATOS_EXPECT_EQ(process.Info(), "GeoApplyConstantScalarValueProcess");
}

} // namespace Kratos::Testing