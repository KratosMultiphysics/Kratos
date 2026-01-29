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

#include "includes/smart_pointers.h"

#include "containers/model.h"
#include "custom_processes/apply_constant_interpolate_line_pressure_process.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/point.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "testing/testing.h"

namespace
{

using namespace Kratos;
using namespace Kratos::Testing;

ModelPart& CreateTestModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("TestPart");
    result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_properties = result.CreateNewProperties(0);

    // Create nodes directly in the model part
    result.CreateNewNode(1, 0.0, 0.0, 0.0);
    result.CreateNewNode(2, 1.0, 0.0, 0.0);
    result.CreateNewNode(3, 1.0, 1.0, 0.0);
    result.CreateNewNode(4, 0.0, 1.0, 0.0);

    // Create elements using the node IDs
    std::vector<ModelPart::IndexType> elem1_nodes = {1, 2, 3};
    std::vector<ModelPart::IndexType> elem2_nodes = {1, 3, 4};
    result.CreateNewElement("Element2D3N", 1, elem1_nodes, p_properties);
    result.CreateNewElement("Element2D3N", 2, elem2_nodes, p_properties);

    return result;
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_Construction, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    Parameters params(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 0.0,
        "table": 1
    })");

    EXPECT_NO_THROW(ApplyConstantInterpolateLinePressureProcess process(r_model_part, params));
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_ThrowsOnInvalidDirections, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    Parameters params(R"({

        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 1,
        "pressure_tension_cut_off": 0.0,
        "table": 1
    })");

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        ApplyConstantInterpolateLinePressureProcess process(r_model_part, params),
        "Gravity direction cannot be the same as Out-of-Plane directions");
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_ExecuteInitializeSolutionStep,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    // Set initial pressure values for boundary nodes
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 10.0 * static_cast<double>(r_node.Id());
    }

    Parameters params(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 100.0,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);
    process.ExecuteInitializeSolutionStep();

    // Check that pressure is set and fixed
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_CHECK(r_node.IsFixed(WATER_PRESSURE))
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 10.0 * r_node.Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_SeepageBranch, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    // Set up so that CalculatePressure returns a value less than the cut-off for node 1,
    // and a value greater than or equal to the cut-off for node 2.
    // We'll do this by setting the initial WATER_PRESSURE and using a high cut-off.

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    }

    // Set a high cut-off so all nodes will be below the cut-off
    Parameters params_below(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": true,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 100.0,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process_below(r_model_part, params_below);
    process_below.ExecuteInitializeSolutionStep();

    // All nodes should be set and fixed (since pressure < cut-off)
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_CHECK(r_node.IsFixed(WATER_PRESSURE));
        // The value is set by CalculatePressure, which in this test setup is 0.0
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    }

    // Now set a low cut-off so all nodes will be above the cut-off
    Parameters params_above(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": true,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": -100.0,
        "table": 1
    })");

    // Reset nodes to unfixed
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Free(WATER_PRESSURE);
    }

    ApplyConstantInterpolateLinePressureProcess process_above(r_model_part, params_above);
    process_above.ExecuteInitializeSolutionStep();

    // All nodes should be free (since pressure >= cut-off)
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_CHECK_IS_FALSE(r_node.IsFixed(WATER_PRESSURE));
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_Info, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    Parameters params(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": false,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 0.0,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);
    KRATOS_EXPECT_EQ(process.Info(), "ApplyConstantInterpolateLinePressureProcess");
}

} // namespace