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

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultiLinePressureProcessThrowsWhenCoordinatesInGravityDirectionAreEmpty,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");

    const auto test_parameters = Parameters{R"(
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
        "Coordinates in gravity direction must not be empty")
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantPhreaticMultiLinePressureProcessThrowsWhenHorizontalCoordinatesAreEmpty, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");

    const auto test_parameters = Parameters{R"(
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
        "Coordinates in horizontal direction must not be empty")
}

} // namespace Kratos::Testing