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

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateThermalFilterLawMatrix, KratosGeoMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");

    auto p_cond_prop = r_model_part.CreateNewProperties(0);
    p_cond_prop->SetValue(THERMAL_CONDUCTIVITY_WATER, 1000.0);

    GeoThermalFilterLaw geo_thermal_filter_law;

    const Matrix thermal_filter_matrix = geo_thermal_filter_law.CalculateThermalDispersionMatrix(*p_cond_prop);

    Matrix expected_solution = ScalarMatrix(1, 1, 1000.0);

    constexpr double tolerance{1.0e-6};

    KRATOS_EXPECT_MATRIX_NEAR(thermal_filter_matrix, expected_solution, tolerance)
}

} // namespace Kratos::Testing
