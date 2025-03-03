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

#include "custom_constitutive/coulomb_yield_surface.hpp"
#include "custom_constitutive/tension_cutoff.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestCoulombYieldSurface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    double friction_angle = 45.0 * Globals::Pi / 180.0;
    double cohesion       = 2.0;

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;

    constexpr double tolerance{1.0e-12};

    CoulombYieldSurface coulombYieldSurface(friction_angle, cohesion, 0.0);
    double              yield_value = coulombYieldSurface.YieldFunctionValue(principal_stress);
    double              expected_solution = 1.0;
    KRATOS_EXPECT_NEAR(yield_value, expected_solution, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(TestTensionCutoff, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    double tension_cutoff = 2.0;

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;

    constexpr double tolerance{1.0e-12};

    TensionCutoff tensionCutoff(tension_cutoff);
    double        tension_value     = tensionCutoff.YieldFunctionValue(principal_stress);
    double        expected_solution = -1.0;
    KRATOS_EXPECT_NEAR(tension_value, expected_solution, tolerance);
}

} // namespace Kratos::Testing
