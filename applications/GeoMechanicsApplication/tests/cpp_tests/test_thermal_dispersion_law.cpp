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
//                   Gennady Markelov
//

#include "custom_constitutive/thermal_dispersion_2D_law.hpp"
#include "geo_mechanics_application.h"
#include "includes/ublas_interface.h"
#include "testing/testing.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(CalculateThermalDispersionMatrix2D, KratosGeoMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    auto cond_prop = r_model_part.CreateNewProperties(0);
    cond_prop->SetValue(POROSITY, 0.5);
    cond_prop->SetValue(SATURATION, 0.75);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_WATER, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XX, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XY, 2000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_YY, 1000.0);

    Matrix ThermalDispersionMatrix = ZeroMatrix(2, 2);
    GeoThermalDispersion2DLaw mGeoThermalDispersion2DLaw;

    mGeoThermalDispersion2DLaw.CalculateThermalDispersionMatrix(
        ThermalDispersionMatrix, *cond_prop);

    Matrix expected_solution = ZeroMatrix(2, 2);
    expected_solution(0, 0) = 875.0;
    expected_solution(0, 1) = 1000.0;
    expected_solution(1, 0) = expected_solution(0, 1);
    expected_solution(1, 1) = expected_solution(0, 0);

    constexpr double tolerance{1.0e-6};

    for (unsigned int i = 0; i < ThermalDispersionMatrix.size1(); ++i) {
        for (unsigned int j = 0; j < ThermalDispersionMatrix.size2(); ++j) {
            KRATOS_EXPECT_NEAR(ThermalDispersionMatrix(i, j),
                               expected_solution(i, j), tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(CalculateThermalDispersionMatrix3D, KratosGeoMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    auto cond_prop = r_model_part.CreateNewProperties(0);
    cond_prop->SetValue(POROSITY, 0.5);
    cond_prop->SetValue(SATURATION, 0.75);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_WATER, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XX, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XY, 2000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_YY, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_XZ, 1000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_YZ, 2000.0);
    cond_prop->SetValue(THERMAL_CONDUCTIVITY_SOLID_ZZ, 500.0);

    Matrix ThermalDispersionMatrix = ZeroMatrix(3, 3);
    GeoThermalDispersion2DLaw mGeoThermalDispersion2DLaw;

    mGeoThermalDispersion2DLaw.CalculateThermalDispersionMatrix(
        ThermalDispersionMatrix, *cond_prop);

    Matrix expected_solution = ZeroMatrix(3, 3);
    expected_solution(0, 0) = 875.0;
    expected_solution(0, 1) = 1000.0;
    expected_solution(0, 2) = 500.0;
    expected_solution(1, 0) = expected_solution(0, 1);
    expected_solution(1, 1) = expected_solution(0, 0);
    expected_solution(1, 2) = 1000.0;
    expected_solution(2, 0) = expected_solution(0, 2);
    expected_solution(2, 1) = expected_solution(1, 2);
    expected_solution(2, 2) = 625.0;

    constexpr double tolerance{1.0e-6};

    for (unsigned int i = 0; i < ThermalDispersionMatrix.size1(); ++i) {
        for (unsigned int j = 0; j < ThermalDispersionMatrix.size2(); ++j) {
            KRATOS_EXPECT_NEAR(ThermalDispersionMatrix(i, j),
                               expected_solution(i, j), tolerance);
        }
    }
}
} // namespace Kratos::Testing
