// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_constitutive/thermal_dispersion_law.h"
#include "geo_mechanics_application.h"
#include "includes/ublas_interface.h"
#include "testing/testing.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(CalculateThermalDispersionMatrix2D, KratosGeoMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");

    auto cond_prop = r_model_part.CreateNewProperties(0);
    cond_prop->SetValue(POROSITY, 0.5);
    cond_prop->SetValue(SATURATED_SATURATION, 0.75);
    cond_prop->SetValue(RETENTION_LAW, "SaturatedLaw");
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_WATER, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XX, 1500.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XY, 2000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_YY, 1200.0);

    const SizeType dimension = 2;
    GeoThermalDispersionLaw geo_thermal_dispersion_2D_law(dimension);
    ProcessInfo info;

    const Matrix thermal_dispersion_matrix =
        geo_thermal_dispersion_2D_law.CalculateThermalDispersionMatrix(
            *cond_prop, info);

    Matrix expected_solution = ZeroMatrix(2, 2);
    expected_solution(0, 0) = 1125.0;
    expected_solution(0, 1) = 1000.0;
    expected_solution(1, 0) = expected_solution(0, 1);
    expected_solution(1, 1) = 975.0;

    constexpr double tolerance{1.0e-6};

    for (unsigned int i = 0; i < thermal_dispersion_matrix.size1(); ++i) {
        for (unsigned int j = 0; j < thermal_dispersion_matrix.size2(); ++j) {
            KRATOS_EXPECT_NEAR(thermal_dispersion_matrix(i, j),
                               expected_solution(i, j), tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(CalculateThermalDispersionMatrix3D, KratosGeoMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");

    auto cond_prop = r_model_part.CreateNewProperties(0);
    cond_prop->SetValue(POROSITY, 0.5);
    cond_prop->SetValue(SATURATED_SATURATION, 0.75);
    cond_prop->SetValue(RETENTION_LAW, "SaturatedLaw");
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_WATER, 800.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XX, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XY, 2000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_YY, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XZ, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_YZ, 2000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_ZZ, 500.0);

    const SizeType dimension = 3;
    GeoThermalDispersionLaw geo_thermal_dispersion_3D_law(dimension);
    ProcessInfo info;

    const Matrix thermal_dispersion_matrix =
        geo_thermal_dispersion_3D_law.CalculateThermalDispersionMatrix(
            *cond_prop, info);

    Matrix expected_solution = ZeroMatrix(3, 3);
    expected_solution(0, 0) = 800.0;
    expected_solution(0, 1) = 1000.0;
    expected_solution(0, 2) = 500.0;
    expected_solution(1, 0) = expected_solution(0, 1);
    expected_solution(1, 1) = expected_solution(0, 0);
    expected_solution(1, 2) = 1000.0;
    expected_solution(2, 0) = expected_solution(0, 2);
    expected_solution(2, 1) = expected_solution(1, 2);
    expected_solution(2, 2) = 550.0;

    constexpr double tolerance{1.0e-6};

    for (unsigned int i = 0; i < thermal_dispersion_matrix.size1(); ++i) {
        for (unsigned int j = 0; j < thermal_dispersion_matrix.size2(); ++j) {
            KRATOS_EXPECT_NEAR(thermal_dispersion_matrix(i, j),
                               expected_solution(i, j), tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(GetWorkingSpaceDimension_ReturnsCorrectValue, KratosGeoMechanicsFastSuite)
{
    constexpr SizeType dimension = 3;
    GeoThermalDispersionLaw geo_thermal_dispersion_3D_law(dimension);

    KRATOS_EXPECT_EQ(geo_thermal_dispersion_3D_law.WorkingSpaceDimension(), dimension);
}

} // namespace Kratos::Testing
