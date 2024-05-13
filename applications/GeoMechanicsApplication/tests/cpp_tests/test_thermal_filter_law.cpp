// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#include "custom_constitutive/thermal_filter_law.h"
#include "geo_mechanics_application.h"
#include "includes/ublas_interface.h"
#include "testing/testing.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(CalculateThermalFilterLawMatrix, KratosGeoMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");

    auto cond_prop = r_model_part.CreateNewProperties(0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_WATER, 1000.0);

    const SizeType dimension = 1;
    GeoThermalDispersionLaw geo_thermal_filter_law(dimension);
    ProcessInfo info;

    const Matrix thermal_filter_matrix =
        geo_thermal_filter_law.CalculateThermalDispersionMatrix(*cond_prop, info);

    Matrix expected_solution = ZeroMatrix(1, 1);
    expected_solution(0, 0)  = 1000.0;

    constexpr double tolerance{1.0e-6};

     KRATOS_EXPECT_NEAR(thermal_filter_matrix(0, 0), expected_solution(0, 0), tolerance);

}

KRATOS_TEST_CASE_IN_SUITE(GetWorkingSpaceDimension_ReturnsFilterDimensionValue, KratosGeoMechanicsFastSuite)
{
    constexpr SizeType dimension = 1;
    GeoThermalDispersionLaw geo_thermal_filter_law(dimension);

    KRATOS_EXPECT_EQ(geo_thermal_filter_law.WorkingSpaceDimension(), dimension);
}

} // namespace Kratos::Testing
