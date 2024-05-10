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

#include "custom_utilities/equation_of_motion_utilities.h"
#include "includes/checks.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateDampingMatrixGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    constexpr std::size_t n = 10;

    constexpr double mass_matrix_value = 10;
    const auto       mass_matrix       = scalar_matrix(n, n, mass_matrix_value);

    constexpr double stiffness_matrix_value = 20;
    const auto       stiffness_matrix       = scalar_matrix(n, n, stiffness_matrix_value);

    double rayleigh_alpha = 0.0;
    double rayleigh_beta  = 1.0;
    auto   damping_matrix = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        rayleigh_alpha, rayleigh_beta, mass_matrix, stiffness_matrix);

    KRATOS_CHECK_MATRIX_NEAR(damping_matrix, stiffness_matrix, 1e-4);

    rayleigh_alpha = 1.0;
    rayleigh_beta  = 0.0;
    damping_matrix = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        rayleigh_alpha, rayleigh_beta, mass_matrix, stiffness_matrix);

    KRATOS_CHECK_MATRIX_NEAR(damping_matrix, mass_matrix, 1e-4);

    rayleigh_alpha = 0.5;
    rayleigh_beta  = 0.5;
    damping_matrix = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        rayleigh_alpha, rayleigh_beta, mass_matrix, stiffness_matrix);

    const double expected_matrix_value = rayleigh_alpha * mass_matrix_value + rayleigh_beta * stiffness_matrix_value;
    const auto expected_damping_matrix = scalar_matrix(n, n, expected_matrix_value);

    KRATOS_CHECK_MATRIX_NEAR(damping_matrix, expected_damping_matrix, 1e-4)
}

} // namespace Kratos::Testing