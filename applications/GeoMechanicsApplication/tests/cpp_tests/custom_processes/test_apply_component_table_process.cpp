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

#include "custom_processes/apply_component_table_process.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include <algorithm>

namespace
{

using namespace Kratos;

void AssertNodesHaveCorrectValueAndFixity(double expected_value, bool expected_fixity, const ModelPart& rModelPart)
{
    for (const auto& r_node : rModelPart.Nodes()) {
        KRATOS_CHECK_EQUAL(r_node.IsFixed(DISPLACEMENT_X), expected_fixity);
        KRATOS_CHECK_EQUAL(r_node.FastGetSolutionStepValue(DISPLACEMENT_X), expected_value);
    }
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyScalarConstraintTableProcess_FreesDoFAfterFinalize_ForDoubleVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      model;
    const auto nodal_variables = Geo::ConstVariableRefs{
        std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y), std::cref(WATER_PRESSURE)};
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model, nodal_variables);
    r_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER] = 1.0;

    auto table = std::make_shared<Table<double>>();
    table->insert(0.0, 0.5);
    table->insert(1.0, 1.5);
    table->SetNameOfX("TIME");
    table->SetNameOfY("DISPLACEMENT_X");
    r_model_part.AddTable(1, table);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
    }

    Parameters parameters(R"(
      {
          "model_part_name": "Main",
          "variable_name":   "DISPLACEMENT_X",
          "is_fixed":        true,
          "table":           1,
          "value":           0.5
      }  )");

    ApplyComponentTableProcess process(r_model_part, parameters);
    process.ExecuteInitialize();

    double expected_value  = 0.5; // Initial value, since we haven't initialized a solution step
    bool   expected_fixity = true;
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part);

    r_model_part.GetProcessInfo()[TIME] = 0.5;
    process.ExecuteInitializeSolutionStep();
    expected_value  = 1.0;
    expected_fixity = true;
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part);

    process.ExecuteFinalize();
    expected_value  = 1.0;
    expected_fixity = false;
    AssertNodesHaveCorrectValueAndFixity(expected_value, expected_fixity, r_model_part);
}

} // namespace Kratos::Testing