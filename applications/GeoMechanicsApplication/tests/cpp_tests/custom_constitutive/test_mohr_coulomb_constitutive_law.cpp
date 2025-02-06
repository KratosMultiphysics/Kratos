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

#include "custom_constitutive/coulomb_yield_function.hpp"
#include "custom_constitutive/tension_cutoff_function.hpp"
#include "geo_mechanics_application.h"
#include "includes/ublas_interface.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateYieldFunctions, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    double friction_angle = 35.0;
    double cohesion = 0.75;
    double tension_cutoff = 1.0;

    Vector principal_stress = ZeroVector(3);
    principal_stress(0) = 135.736961146391;
    principal_stress(1) = 22.5224297324582;
    principal_stress(2) = -8.25939087884923;

    constexpr double tolerance{1.0e-12};

    CoulombYieldFunction coulombYieldFunction(friction_angle, cohesion);
    double yield_value = coulombYieldFunction(principal_stress);
    double expected_solution = 34.824746745030032;
    KRATOS_EXPECT_NEAR(yield_value, expected_solution, tolerance);

    TensionCutoffFunction tensionCutoffFunction(tension_cutoff);
    double tension_value = tensionCutoffFunction(principal_stress);
    expected_solution = 9.25939087884923;
    KRATOS_EXPECT_NEAR(tension_value, expected_solution, tolerance);
}

} // namespace Kratos::Testing
