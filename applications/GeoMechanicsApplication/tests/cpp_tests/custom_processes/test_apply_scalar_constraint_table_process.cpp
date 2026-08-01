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
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include "geo_mechanics_application_variables.h"

#include <string>

namespace
{

using namespace Kratos;
using namespace Kratos::Testing;
using namespace std::string_literals;

template <typename T>
void AssertNodesHaveCorrectValueAndFixity(const Variable<T>&                   rVariable,
                                          double                               ExpectedValue,
                                          bool                                 ExpectedFixity,
                                          const ModelPart::NodesContainerType& rNodes)
{
    for (const auto& r_node : rNodes) {
        KRATOS_EXPECT_EQ(r_node.IsFixed(rVariable), ExpectedFixity);
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(rVariable), ExpectedValue);
    }
}

ModelPart& SetupModelPartWith2D3NElement(const Table<double>::Pointer& rTable, Model& rModel, const std::string& ModelPartName)
{
    const auto nodal_variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(WATER_PRESSURE)};
    auto& r_model_part =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(rModel, nodal_variables, ModelPartName);
    r_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER] = 1.0;
    if (rTable) r_model_part.AddTable(1, rTable);

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
    auto  p_table = std::make_shared<Table<double>>();
    p_table->SetNameOfX("TIME"); // Table can be minimal, since we only do Initialize and Finalize
    auto& r_first_model_part  = SetupModelPartWith2D3NElement(p_table, model, "First"s);
    auto& r_second_model_part = SetupModelPartWith2D3NElement(p_table, model, "Second"s);

    Parameters parameters(R"(
      {
          "model_part_name_list": ["First", "Second"],
          "variable_name":   "DISPLACEMENT_X",
          "is_fixed":        true,
          "table":           1,
          "value":           0.5
      }  )");

    ApplyScalarConstraintTableProcess process(model, parameters);

    // Act
    process.ExecuteInitialize();
    process.ExecuteFinalize();

    // Assert
    constexpr double expected_value =
        0.5; // Same as the initial value, since we have not initialized any solution step
    constexpr bool expected_fixity = false;
    AssertNodesHaveCorrectValueAndFixity(DISPLACEMENT_X, expected_value, expected_fixity,
                                         r_first_model_part.Nodes());
    AssertNodesHaveCorrectValueAndFixity(DISPLACEMENT_X, expected_value, expected_fixity,
                                         r_second_model_part.Nodes());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_AppliesCorrectValuesThroughTime_ForDoubleVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  p_table = std::make_shared<Table<double>>();
    p_table->insert(0.0, 0.5);
    p_table->insert(1.0, 1.5);
    p_table->SetNameOfX("TIME");
    p_table->SetNameOfY("DISPLACEMENT_X");
    auto& r_model_part = SetupModelPartWith2D3NElement(p_table, model, "Main"s);

    Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "DISPLACEMENT_X",
          "is_fixed":        true,
          "table":           1,
          "value":           0.3
      }  )");

    ApplyScalarConstraintTableProcess process(model, parameters);

    // Act & Assert
    process.ExecuteInitialize();
    double expected_value = 0.3; // Initial value, since we haven't initialized a solution step
    constexpr bool expected_fixity = true;
    AssertNodesHaveCorrectValueAndFixity(DISPLACEMENT_X, expected_value, expected_fixity,
                                         r_model_part.Nodes());

    r_model_part.GetProcessInfo()[TIME] = 0.5;
    process.ExecuteInitializeSolutionStep();
    expected_value = 1.0;
    AssertNodesHaveCorrectValueAndFixity(DISPLACEMENT_X, expected_value, expected_fixity,
                                         r_model_part.Nodes());

    r_model_part.GetProcessInfo()[TIME] = 0.8;
    process.ExecuteInitializeSolutionStep();
    expected_value = 1.3;
    AssertNodesHaveCorrectValueAndFixity(DISPLACEMENT_X, expected_value, expected_fixity,
                                         r_model_part.Nodes());

    r_model_part.GetProcessInfo()[TIME] = 2.0;
    process.ExecuteInitializeSolutionStep();
    expected_value = 2.5; // Extrapolated value
    AssertNodesHaveCorrectValueAndFixity(DISPLACEMENT_X, expected_value, expected_fixity,
                                         r_model_part.Nodes());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_CheckInfoApplyScalarConstraintTableProcess,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  p_table = std::make_shared<Table<double>>();
    SetupModelPartWith2D3NElement(p_table, model, "Main"s);

    const Parameters                        parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "DISPLACEMENT_X",
          "is_fixed":        true,
          "table":           1,
          "value":           0.3
      }  )");
    const ApplyScalarConstraintTableProcess process(model, parameters);
    // Act & assert
    KRATOS_EXPECT_EQ(process.Info(), "ApplyScalarConstraintTableProcess");
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_UniformFluidPressure_NoTable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = SetupModelPartWith2D3NElement(nullptr, model, "Main"s);

    const Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "WATER_PRESSURE",
          "is_fixed":        false,
          "table":           0,
          "value":           10500.0,
          "fluid_pressure_type": "Uniform"
      }  )");

    ApplyScalarConstraintTableProcess process(model, parameters);

    // Act
    process.ExecuteInitialize();
    r_model_part.GetProcessInfo()[TIME] = 5.0;
    process.ExecuteInitializeSolutionStep();

    // Assert
    constexpr auto expected_value  = 10500.0;
    constexpr auto expected_fixity = false;
    AssertNodesHaveCorrectValueAndFixity(WATER_PRESSURE, expected_value, expected_fixity,
                                         r_model_part.Nodes());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_GenericVariableWithTable, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  p_table = std::make_shared<Table<double>>();
    p_table->insert(0.0, 1.0);
    p_table->insert(10.0, 5.0);
    p_table->SetNameOfX("TIME");
    p_table->SetNameOfY("WATER_PRESSURE");
    auto& r_model_part = SetupModelPartWith2D3NElement(p_table, model, "Main"s);

    const Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "WATER_PRESSURE",
          "table":           1,
          "value":           0.0
      }  )");

    ApplyScalarConstraintTableProcess process(model, parameters);

    // Act
    process.ExecuteInitialize();
    r_model_part.GetProcessInfo()[TIME] = 5.0;
    process.ExecuteInitializeSolutionStep();

    // Assert
    const auto     expected_value  = 3.0; // Linear interpolation between 1.0 and 5.0 at TIME=5.0
    constexpr auto expected_fixity = false;
    AssertNodesHaveCorrectValueAndFixity(WATER_PRESSURE, expected_value, expected_fixity,
                                         r_model_part.Nodes());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_ErrorOnUnknownFluidType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    SetupModelPartWith2D3NElement(nullptr, model, "Main"s);

    const Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "WATER_PRESSURE",
          "fluid_pressure_type": "NonExistent_Type"
      }  )");

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ApplyScalarConstraintTableProcess process(model, parameters),
                                      "Unknown fluid_pressure_type: NonExistent_Type");
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_HydrostaticBranch_Dispatch,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    SetupModelPartWith2D3NElement(nullptr, model, "Main"s);

    const Parameters parameters(R"(
      {
       "model_part_name": "Main",
        "fluid_pressure_type": "Hydrostatic",
        "gravity_direction": 1,
        "is_fixed": false,
        "is_seepage": false,
        "pressure_tension_cut_off": 0.0,
        "reference_coordinate": 0.0,
        "specific_weight": 10000.0,
        "table": 0,
        "variable_name": "WATER_PRESSURE",
        "value":           10500.0
      }  )");

    // Act
    ApplyScalarConstraintTableProcess process(model, parameters);

    // Assert
    EXPECT_NO_THROW(process.ExecuteInitialize());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_InterpolateLine_NoTable_FixedValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    SetupModelPartWith2D3NElement(nullptr, model, "Main"s);

    Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "WATER_PRESSURE",
          "table": 0,
          "fluid_pressure_type": "Interpolate_Line",
          "gravity_direction": 1,
          "out_of_plane_direction": 2,
          "value": 2500.0
      }  )");

    ApplyScalarConstraintTableProcess process(model, parameters);

    // Act and Assert
    EXPECT_NO_THROW(process.ExecuteInitialize());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_PhreaticMultiLine_CoordinatesParsing,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    SetupModelPartWith2D3NElement(nullptr, model, "Main"s);

    Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "WATER_PRESSURE",
          "table": 0,
          "fluid_pressure_type": "Phreatic_Multi_Line",
          "gravity_direction": 1,
          "out_of_plane_direction": 2,
          "x_coordinates": [0.0, 5.0, 10.0],
          "y_coordinates": [10.0, 8.0, 5.0],
          "z_coordinates": [0.0, 0.0, 0.0],
          "specific_weight": 9810.0
      }  )");

    ApplyScalarConstraintTableProcess process(model, parameters);

    // Act & Assert
    EXPECT_NO_THROW(process.ExecuteInitialize());
}

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_PhreaticLine_BranchCreation,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    SetupModelPartWith2D3NElement(nullptr, model, "Main"s);

    Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "WATER_PRESSURE",
          "table": 0,
          "fluid_pressure_type": "Phreatic_Line",
          "gravity_direction": 1,
          "out_of_plane_direction": 2,
          "first_reference_coordinate": [0.0, 0.0, 0.0],
          "second_reference_coordinate": [1.0, 1.0, 5.0],
          "specific_weight": 10000.0,
          "is_fixed": false
      }  )");

    ApplyScalarConstraintTableProcess process(model, parameters);

    // Act & Assert
    EXPECT_NO_THROW(process.ExecuteInitialize());
}

} // namespace Kratos::Testing