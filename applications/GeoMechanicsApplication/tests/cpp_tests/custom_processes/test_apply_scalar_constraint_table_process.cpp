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

#include "custom_processes/apply_scalar_constraint_table_process.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include "geo_mechanics_application_variables.h"

namespace
{

using namespace Kratos;
using namespace Kratos::Testing;

void AssertNodesHaveCorrectValueAndFixity(double                               ExpectedValue,
                                          bool                                 ExpectedFixity,
                                          const ModelPart::NodesContainerType& rNodes)
{
    for (const auto& r_node : rNodes) {
        KRATOS_EXPECT_EQ(r_node.IsFixed(DISPLACEMENT_X), ExpectedFixity);
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(DISPLACEMENT_X), ExpectedValue);
    }
}

ModelPart& SetupModelPart(const Table<double>::Pointer& rTable, Model& model)
{
    const auto nodal_variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X)};
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model, nodal_variables);
    r_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER] = 1.0;
    r_model_part.AddTable(1, rTable);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
    }

    return r_model_part;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_FreesDoFAfterFinalize_ForDoubleVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  table = std::make_shared<Table<double>>();
    table->SetNameOfX("TIME"); // Table can be minimal, since we only do Initialize and Finalize
    auto& r_model_part = SetupModelPart(table, model);

    Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "DISPLACEMENT_X",
          "is_fixed":        true,
          "table":           1,
          "value":           0.5
      }  )");

    ApplyScalarConstraintTableProcess process(r_model_part, parameters);

    // Act
    process.ExecuteInitialize();
    process.ExecuteFinalize();

    // Assert
    constexpr double expected_value =
        0.5; // Same as the initial value, since we have not initialized any solution step
    constexpr bool expected_fixity = false;
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part.Nodes());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_AppliesCorrectValuesThroughTime_ForDoubleVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  table = std::make_shared<Table<double>>();
    table->insert(0.0, 0.5);
    table->insert(1.0, 1.5);
    table->SetNameOfX("TIME");
    table->SetNameOfY("DISPLACEMENT_X");
    auto& r_model_part = SetupModelPart(table, model);

    Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "DISPLACEMENT_X",
          "is_fixed":        true,
          "table":           1,
          "value":           0.3
      }  )");

    ApplyScalarConstraintTableProcess process(r_model_part, parameters);

    // Act & Assert
    process.ExecuteInitialize();
    double expected_value = 0.3; // Initial value, since we haven't initialized a solution step
    constexpr bool expected_fixity = true;
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part.Nodes());

    r_model_part.GetProcessInfo()[TIME] = 0.5;
    process.ExecuteInitializeSolutionStep();
    expected_value = 1.0;
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part.Nodes());

    r_model_part.GetProcessInfo()[TIME] = 0.8;
    process.ExecuteInitializeSolutionStep();
    expected_value = 1.3;
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part.Nodes());

    r_model_part.GetProcessInfo()[TIME] = 2.0;
    process.ExecuteInitializeSolutionStep();
    expected_value = 2.5; // Extrapolated value
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part.Nodes());
}

} // namespace Kratos::Testing