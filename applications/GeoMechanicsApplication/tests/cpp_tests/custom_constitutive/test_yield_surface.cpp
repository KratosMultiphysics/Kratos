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

#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/tension_cutoff.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include "utilities/math_utils.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestCoulombYieldSurface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    double           friction_angle = MathUtils<>::DegreesToRadians(45.0);
    constexpr double cohesion       = 2.0;

    CoulombYieldSurface coulombYieldSurface(friction_angle, cohesion, 0.0);

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;
    double yield_value       = coulombYieldSurface.YieldFunctionValue(principal_stress);
    double expected_solution = 1.0;
    KRATOS_EXPECT_NEAR(yield_value, expected_solution, Defaults::absolute_tolerance);

    principal_stress <<= 1.7071067811865475, 1.0, 0.2928932188134525;
    yield_value       = coulombYieldSurface.YieldFunctionValue(principal_stress);
    expected_solution = 0.0;
    KRATOS_EXPECT_NEAR(yield_value, expected_solution, Defaults::absolute_tolerance);

    principal_stress <<= 0.1715728752538099, -1.0, -1.8284271247461901;
    yield_value       = coulombYieldSurface.YieldFunctionValue(principal_stress);
    expected_solution = -1.0;
    KRATOS_EXPECT_NEAR(yield_value, expected_solution, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(TestTensionCutoff, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    double tension_cutoff = 2.0;

    TensionCutoff tensionCutoff(tension_cutoff);

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;
    double tension_value     = tensionCutoff.YieldFunctionValue(principal_stress);
    double expected_solution = 1.0;
    KRATOS_EXPECT_NEAR(tension_value, expected_solution, Defaults::absolute_tolerance);

    principal_stress <<= 2.0, 1.5, 1.0;
    tension_value     = tensionCutoff.YieldFunctionValue(principal_stress);
    expected_solution = 0.0;
    KRATOS_EXPECT_NEAR(tension_value, expected_solution, Defaults::absolute_tolerance);

    principal_stress <<= 1.0, 0.5, 0.1;
    tension_value     = tensionCutoff.YieldFunctionValue(principal_stress);
    expected_solution = -1.0;
    KRATOS_EXPECT_NEAR(tension_value, expected_solution, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
