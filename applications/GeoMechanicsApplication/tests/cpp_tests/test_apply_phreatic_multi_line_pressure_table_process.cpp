// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//
#include "containers/model.h"
#include "custom_processes/apply_phreatic_multi_line_pressure_table_process.h"
#include "geo_mechanics_fast_suite.h"
#include "includes/checks.h"

using namespace Kratos;

namespace
{

bool CanCreateInstanceOfApplyPhreaticMultiLinePressureTableProcessWithoutFailure(ModelPart& rModelPart,
                                                                                 const Parameters& rProcessParameters)
{
    try {
        ApplyPhreaticMultiLinePressureTableProcess{rModelPart, rProcessParameters};
    } catch (const Exception&) {
        return false;
    }

    return true;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyPhreaticMultiLinePressureTableProcessThrowsWhenTableLengthDoesNotMatchCoordinates,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model           = Model{};
    auto& r_model_part    = model.CreateModelPart("foo");
    auto  test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "variable_name": "WATER_PRESSURE",
                "x_coordinates": [0.0, 1.0],
                "y_coordinates": [0.0, 1.0],
                "z_coordinates": [0.0, 0.0],
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "table": [1, 2, 3, 4]
            }  )"};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyPhreaticMultiLinePressureTableProcess{r_model_part, test_parameters}), "Got 2 coordinates and 4 table references. The number of coordinates and table references should be equal.")
}

KRATOS_TEST_CASE_IN_SUITE(ApplyPhreaticMultiLinePressureTableProcessDoesNotThrowWhenTableLengthMatchesCoordinates,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model           = Model{};
    auto& r_model_part    = model.CreateModelPart("foo");
    auto  test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "variable_name": "WATER_PRESSURE",
                "x_coordinates": [0.0, 1.0, 2.0],
                "y_coordinates": [0.0, 1.0, 2.0],
                "z_coordinates": [0.0, 0.0, 0.0],
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "table": [0, 0, 3]
            }  )"};

    KRATOS_EXPECT_TRUE(CanCreateInstanceOfApplyPhreaticMultiLinePressureTableProcessWithoutFailure(
        r_model_part, test_parameters))
}

} // namespace Kratos::Testing