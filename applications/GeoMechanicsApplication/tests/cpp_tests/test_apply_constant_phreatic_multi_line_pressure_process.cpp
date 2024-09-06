// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//
#include "containers/model.h"
#include "custom_processes/apply_constant_phreatic_multi_line_pressure_process.h"
#include "geo_mechanics_fast_suite.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/checks.h"
#include "processes/structured_mesh_generator_process.h"

using namespace Kratos;

namespace
{

bool CanCreateInstanceOfApplyConstantPhreaticMultiLinePressureProcessWithoutFailure(ModelPart& rModelPart,
                                                                                    const Parameters& rProcessParameters)
{
    try {
        ApplyConstantPhreaticMultiLinePressureProcess{rModelPart, rProcessParameters};
    } catch (const Exception&) {
        return false;
    }

    return true;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultiLinePressureProcessThrowsWhenLessThanTwoCoordinatesInGravityDirectionAreGiven,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model           = Model{};
    auto& r_model_part    = model.CreateModelPart("foo");
    auto  test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "variable_name": "WATER_PRESSURE",
                "x_coordinates": [0.0, 1.0],
                "y_coordinates": [],
                "z_coordinates": [0.0, 0.0],
                "gravity_direction": 1
            }  )"};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyConstantPhreaticMultiLinePressureProcess{r_model_part, test_parameters}),
        "At least two coordinates in gravity direction must be given, but got 0")

    test_parameters.RemoveValue("y_coordinates");
    test_parameters.AddVector("y_coordinates", ScalarVector{1, 0.0});
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyConstantPhreaticMultiLinePressureProcess{r_model_part, test_parameters}),
        "At least two coordinates in gravity direction must be given, but got 1")
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultiLinePressureProcessThrowsWhenLessThanTwoCoordinatesInHorizontalDirectionAreGiven,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model           = Model{};
    auto& r_model_part    = model.CreateModelPart("foo");
    auto  test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "variable_name": "WATER_PRESSURE",
                "x_coordinates": [],
                "y_coordinates": [0.0, 1.0],
                "z_coordinates": [0.0, 0.0],
                "gravity_direction": 1,
                "out_of_plane_direction": 2
            }  )"};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyConstantPhreaticMultiLinePressureProcess{r_model_part, test_parameters}),
        "At least two coordinates in horizontal direction must be given, but got 0")

    test_parameters.RemoveValue("x_coordinates");
    test_parameters.AddVector("x_coordinates", ScalarVector{1, 0.0});
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyConstantPhreaticMultiLinePressureProcess{r_model_part, test_parameters}),
        "At least two coordinates in horizontal direction must be given, but got 1")
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultiLinePressureProcessDoesNotThrowWhenAtLeastTwoCoordinatesInGravityDirectionAreGiven,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model           = Model{};
    auto& r_model_part    = model.CreateModelPart("foo");
    auto  test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "variable_name": "WATER_PRESSURE",
                "x_coordinates": [0.0, 1.0],
                "y_coordinates": [1.0, 2.0],
                "z_coordinates": [0.0, 0.0],
                "gravity_direction": 1
            }  )"};

    KRATOS_EXPECT_TRUE(CanCreateInstanceOfApplyConstantPhreaticMultiLinePressureProcessWithoutFailure(
        r_model_part, test_parameters))

    test_parameters.RemoveValue("y_coordinates");
    test_parameters.AddVector("y_coordinates", ScalarVector{5, 1.0});

    KRATOS_EXPECT_TRUE(CanCreateInstanceOfApplyConstantPhreaticMultiLinePressureProcessWithoutFailure(
        r_model_part, test_parameters))
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultiLinePressureProcessDoesNotThrowWhenAtLeastTwoCoordinatesInHorizontalDirectionAreGiven,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model           = Model{};
    auto& r_model_part    = model.CreateModelPart("foo");
    auto  test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "variable_name": "WATER_PRESSURE",
                "x_coordinates": [1.0, 2.0],
                "y_coordinates": [0.0, 1.0],
                "z_coordinates": [0.0, 0.0],
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "table": [0, 0]
            }  )"};

    KRATOS_EXPECT_TRUE(CanCreateInstanceOfApplyConstantPhreaticMultiLinePressureProcessWithoutFailure(
        r_model_part, test_parameters))

    test_parameters.RemoveValue("x_coordinates");
    test_parameters.AddVector("x_coordinates", ScalarVector{2, 1.0});

    KRATOS_EXPECT_TRUE(CanCreateInstanceOfApplyConstantPhreaticMultiLinePressureProcessWithoutFailure(
        r_model_part, test_parameters))
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultilinePressureProcess_AppliesCorrectPressure_WhenPhreaticLineDoesNotSpanEntireDomain,
                          KratosGeoMechanicsFastSuite)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Set up the test model part mesh
    auto                   p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto                   p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    auto                   p_point_3 = Kratos::make_intrusive<Node>(3, 5.0, 1.0, 0.0);
    auto                   p_point_4 = Kratos::make_intrusive<Node>(4, 5.0, 0.0, 0.0);
    Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters             mesher_parameters(R"({
        "number_of_divisions": 2,
        "element_name": "Element2D3N",
        "condition_name": "LineCondition",
        "create_skin_sub_model_part": true
    })");
    StructuredMeshGeneratorProcess(domain_geometry, r_model_part, mesher_parameters).Execute();

    auto test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "variable_name": "WATER_PRESSURE",
                "specific_weight": 10000.0,
                "x_coordinates": [0.0, 1.0, 2.0],
                "y_coordinates": [1.0, 1.0, 1.0],
                "z_coordinates": [0.0, 0.0, 0.0],
                "table": [0, 0, 0],
                "gravity_direction": 1
            }  )"};

    ApplyConstantPhreaticMultiLinePressureProcess process{r_model_part, test_parameters};
    process.ExecuteInitialize();

    const auto phreatic_line_position = 1.0;
    const auto specific_weight        = 10000.0;
    block_for_each(r_model_part.Nodes(), [&phreatic_line_position, &specific_weight](auto& node) {
        const auto water_pressure = (node.Y() - phreatic_line_position) * specific_weight;
        KRATOS_EXPECT_DOUBLE_EQ(node.FastGetSolutionStepValue(WATER_PRESSURE, 0), water_pressure);
    });
}

} // namespace Kratos::Testing