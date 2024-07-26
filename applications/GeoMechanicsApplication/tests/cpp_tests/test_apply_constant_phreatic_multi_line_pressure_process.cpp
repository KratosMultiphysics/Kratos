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
#include "includes/checks.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultiLinePressureProcessThrowsWhenLessThanTwoCoordinatesInGravityDirectionAreGiven,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");

    auto test_parameters = Parameters{R"(
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
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");

    auto test_parameters = Parameters{R"(
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

} // namespace Kratos::Testing