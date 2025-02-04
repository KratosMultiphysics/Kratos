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

#include "custom_constitutive/mohr_coulomb_constitutive_law.hpp"
#include "geo_mechanics_application.h"
#include "includes/ublas_interface.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateYieldFunction, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");

    auto cond_prop = r_model_part.CreateNewProperties(0);
    cond_prop->SetValue(GEO_FRICTION_ANGLE, 35.0);
    cond_prop->SetValue(GEO_COHESION, 0.75);
    cond_prop->SetValue(GEO_TENSION_CUTOFF, 1.0);

    Vector principal_stress = ZeroVector(3);
    principal_stress(0) = 135.736961146391;
    principal_stress(1) = 22.5224297324582;
    principal_stress(2) = -8.25939087884923;

    constexpr double tolerance{1.0e-6};

    MohrCoulombConstitutiveLaw mohr_coulomb_constitutive_law;
    double yield_value = mohr_coulomb_constitutive_law.CalculateCoulombYieldFunction(principal_stress, *cond_prop);
    double expected_solution = 34.824746745030032;
    KRATOS_EXPECT_NEAR(yield_value, expected_solution, tolerance);

    double tension_value = mohr_coulomb_constitutive_law.CalculateTensionYieldFunction(principal_stress, *cond_prop);
    expected_solution = 9.25939087884923;
    KRATOS_EXPECT_NEAR(tension_value, expected_solution, tolerance);
}


KRATOS_TEST_CASE_IN_SUITE(CalculatePrinsipalStresses3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector cauchy_stresses = ZeroVector(6);
    cauchy_stresses(0) = 80.0;
    cauchy_stresses(1) = 50.0;
    cauchy_stresses(2) = 20.0;
    cauchy_stresses(3) = 40.0;
    cauchy_stresses(4) = 35.0;
    cauchy_stresses(5) = 45.0;

    MohrCoulombConstitutiveLaw mohr_coulomb_constitutive_law;
    Vector principal_stresses = mohr_coulomb_constitutive_law.CalculatePrincipalStresses(cauchy_stresses);

    Vector expected_solution = ZeroVector(3);
    expected_solution(0) = 135.736961146391;
    expected_solution(1) = 22.5224297324582;
    expected_solution(2) = -8.25939087884923;

    constexpr double tolerance{1.0e-6};

    for (unsigned int i = 0; i < principal_stresses.size(); ++i) {
        KRATOS_EXPECT_NEAR(principal_stresses(i), expected_solution(i), tolerance);
    }
}

} // namespace Kratos::Testing
