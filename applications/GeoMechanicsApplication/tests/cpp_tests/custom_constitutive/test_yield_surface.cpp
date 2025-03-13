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
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include "utilities/math_utils.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestCoulombYieldSurface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto       friction_angle  = MathUtils<>::DegreesToRadians(45.0);
    constexpr double cohesion        = 2.0;
    const auto       dilatancy_angle = MathUtils<>::DegreesToRadians(0.0);

    CoulombYieldSurface coulomb_yield_surface(friction_angle, cohesion, dilatancy_angle);

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(principal_stress), 1.0, Defaults::absolute_tolerance);

    principal_stress <<= 1.7071067811865475, 1.0, 0.2928932188134525;
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(principal_stress), 0.0, Defaults::absolute_tolerance);

    principal_stress <<= 0.1715728752538099, -1.0, -1.8284271247461901;
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(principal_stress), -1.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CoulombYieldSurface_CanBeSavedAndLoadedThroughInterface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     scoped_registration     = ScopedSerializerRegistrationOfAllYieldSurfaces{};
    const auto     friction_angle          = MathUtils<>::DegreesToRadians(60.0);
    constexpr auto cohesion                = 2.0;
    const auto     dilatancy_angle         = MathUtils<>::DegreesToRadians(30.0);
    auto           p_coulomb_yield_surface = std::unique_ptr<YieldSurface>(
        std::make_unique<CoulombYieldSurface>(friction_angle, cohesion, dilatancy_angle));
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_coulomb_yield_surface);
    auto p_loaded_coulomb_yield_surface = std::unique_ptr<YieldSurface>{};
    serializer.load("test_tag"s, p_loaded_coulomb_yield_surface);

    // Assert
    ASSERT_NE(p_loaded_coulomb_yield_surface, nullptr);
    auto principal_stresses = Vector(3);
    principal_stresses <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_NEAR(p_loaded_coulomb_yield_surface->YieldFunctionValue(principal_stresses),
                       0.5 * std::sqrt(3.0) - 1, Defaults::absolute_tolerance);
    auto expected_derivative = Vector(3);
    expected_derivative <<= 0.75, 0.0, -0.25;
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_coulomb_yield_surface->DerivativeOfFlowFunction(principal_stresses),
                              expected_derivative, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(TestTensionCutoff, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    double tension_cutoff = 2.0;

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;

    TensionCutoff tensionCutoff(tension_cutoff);
    KRATOS_EXPECT_NEAR(tensionCutoff.YieldFunctionValue(principal_stress), 1.0, Defaults::absolute_tolerance);

    principal_stress <<= 2.0, 1.5, 1.0;
    KRATOS_EXPECT_NEAR(tensionCutoff.YieldFunctionValue(principal_stress), 0.0, Defaults::absolute_tolerance);

    principal_stress <<= 1.0, 0.5, 0.1;
    KRATOS_EXPECT_NEAR(tensionCutoff.YieldFunctionValue(principal_stress), -1.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(TensionCutOff_CanBeSavedAndLoadedThroughInterface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto     scoped_registration = ScopedSerializerRegistrationOfAllYieldSurfaces{};
    constexpr auto tensile_strength    = 2.0;
    auto p_tension_cut_off = std::unique_ptr<YieldSurface>(std::make_unique<TensionCutoff>(tensile_strength));
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_tension_cut_off);
    auto p_loaded_tension_cut_off = std::unique_ptr<YieldSurface>{};
    serializer.load("test_tag"s, p_loaded_tension_cut_off);

    // Assert
    ASSERT_NE(p_loaded_tension_cut_off, nullptr);
    auto principal_stresses = Vector(3);
    principal_stresses <<= tensile_strength, 0.0, 0.0;
    KRATOS_EXPECT_NEAR(p_loaded_tension_cut_off->YieldFunctionValue(principal_stresses), 0.0,
                       Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
