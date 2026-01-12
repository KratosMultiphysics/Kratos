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

#include "custom_constitutive/compression_cap_yield_surface.h"
#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/tension_cutoff.h"
#include "custom_utilities/registration_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
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
    auto material_properties                 = Properties{};
    material_properties[GEO_FRICTION_ANGLE]  = 45.0;
    material_properties[GEO_COHESION]        = 2.0;
    material_properties[GEO_DILATANCY_ANGLE] = 0.0;

    const auto coulomb_yield_surface = CoulombYieldSurface{material_properties};

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;
    auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress);
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(sigma_tau), 1.0, Defaults::absolute_tolerance);

    principal_stress <<= 1.7071067811865475, 1.0, 0.2928932188134525;
    sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress);
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(sigma_tau), 0.0, Defaults::absolute_tolerance);

    principal_stress <<= 0.1715728752538099, -1.0, -1.8284271247461901;
    sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress);
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(sigma_tau), -1.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CoulombYieldSurface_CanBeSavedAndLoadedThroughInterface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("CoulombYieldSurface"s, CoulombYieldSurface{})};

    auto material_properties                 = Properties{};
    material_properties[GEO_FRICTION_ANGLE]  = 60.0;
    material_properties[GEO_COHESION]        = 2.0;
    material_properties[GEO_DILATANCY_ANGLE] = 30.0;
    auto p_coulomb_yield_surface =
        std::unique_ptr<YieldSurface>{std::make_unique<CoulombYieldSurface>(material_properties)};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_coulomb_yield_surface);
    auto p_loaded_coulomb_yield_surface = std::unique_ptr<YieldSurface>{};
    serializer.load("test_tag"s, p_loaded_coulomb_yield_surface);

    // Assert
    ASSERT_NE(p_loaded_coulomb_yield_surface, nullptr);
    auto principal_stresses = Vector(3);
    principal_stresses <<= 1.0, 1.0, 1.0;
    const auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stresses);
    KRATOS_EXPECT_NEAR(p_loaded_coulomb_yield_surface->YieldFunctionValue(sigma_tau),
                       0.5 * std::sqrt(3.0) - 1, Defaults::absolute_tolerance);
    auto expected_derivative = Vector(2);
    expected_derivative <<= 0.5, 1.0;
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_coulomb_yield_surface->DerivativeOfFlowFunction(sigma_tau),
                              expected_derivative, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(TestTensionCutoff, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    double tension_cutoff = 2.0;

    Vector principal_stress(3);
    principal_stress <<= 3.0, 2.0, 1.0;
    auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress);
    TensionCutoff tensionCutoff(tension_cutoff);
    KRATOS_EXPECT_NEAR(tensionCutoff.YieldFunctionValue(sigma_tau), 1.0, Defaults::absolute_tolerance);

    principal_stress <<= 2.0, 1.5, 1.0;
    sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress);
    KRATOS_EXPECT_NEAR(tensionCutoff.YieldFunctionValue(sigma_tau), 0.0, Defaults::absolute_tolerance);

    principal_stress <<= 1.0, 0.5, 0.1;
    sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress);
    KRATOS_EXPECT_NEAR(tensionCutoff.YieldFunctionValue(sigma_tau), -1.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(TensionCutOff_CanBeSavedAndLoadedThroughInterface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("TensionCutoff"s, TensionCutoff{})};
    constexpr auto tensile_strength = 2.0;
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
    const auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stresses);
    KRATOS_EXPECT_NEAR(p_loaded_tension_cut_off->YieldFunctionValue(sigma_tau), 0.0, Defaults::absolute_tolerance);
    auto expected_derivative = Vector(2);
    expected_derivative <<= 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_tension_cut_off->DerivativeOfFlowFunction(sigma_tau),
                              expected_derivative, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CoulombYieldSurface_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties(3);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CoulombYieldSurface{properties},
                                      "GEO_COHESION does not exist in the property with Id 3.")
    properties.SetValue(GEO_COHESION, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CoulombYieldSurface{properties},
                                      "GEO_COHESION in the property with Id 3 has an invalid "
                                      "value: -1 is out of the range [0, -).")
    properties.SetValue(GEO_COHESION, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "GEO_FRICTION_ANGLE does not exist in the property with Id 3.")
    properties.SetValue(GEO_FRICTION_ANGLE, -30.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CoulombYieldSurface{properties},
                                      "GEO_FRICTION_ANGLE in the property with Id 3 has an invalid "
                                      "value: -30 is out of the range [0, -).")
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "GEO_DILATANCY_ANGLE does not exist in the property with Id 3.")
    properties.SetValue(GEO_DILATANCY_ANGLE, -30.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CoulombYieldSurface{properties},
                                      "GEO_DILATANCY_ANGLE in the property with Id 3 has an "
                                      "invalid value: -30 is out of the range [0, 30].")
    properties.SetValue(GEO_DILATANCY_ANGLE, 40.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CoulombYieldSurface{properties},
                                      " GEO_DILATANCY_ANGLE in the property with Id 3 has an "
                                      "invalid value: 40 is out of the range [0, 30].")
    properties.SetValue(GEO_DILATANCY_ANGLE, 30.0);

    properties.SetValue(GEO_COULOMB_HARDENING_TYPE, "Linear");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "GEO_COHESION_FUNCTION_COEFFICIENTS does not exist in the property with Id 3.")
    properties.SetValue(GEO_COHESION_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({-1.0}));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "Error: Entry 0 in GEO_COHESION_FUNCTION_COEFFICIENTS out of range. Value: -1")
    properties.SetValue(GEO_COHESION_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({1.0}));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS does not exist in the property with Id 3.")
    properties.SetValue(GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({-1.0}));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "Error: Entry 0 in GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS out of range. Value: -1")
    properties.SetValue(GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({1.0}));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "GEO_DILATANCY_ANGLE_FUNCTION_COEFFICIENTS does not exist in the property with Id 3.")
    properties.SetValue(GEO_DILATANCY_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({-1.0}));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CoulombYieldSurface{properties},
        "Error: Entry 0 in GEO_DILATANCY_ANGLE_FUNCTION_COEFFICIENTS out of range. Value: -1")
}

KRATOS_TEST_CASE_IN_SUITE(TestCompressionCapYieldSurface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto material_properties                          = Properties{};
    material_properties[GEO_CAP_HARDENING_TYPE]       = "None";
    material_properties[GEO_COMPRESSION_CAP_SIZE]     = 4.0;
    material_properties[GEO_COMPRESSION_CAP_LOCATION] = 20.0;

    const auto cap_yield_surface = CompressionCapYieldSurface{material_properties};

    Vector principal_stress(3);
    principal_stress <<= 30.0, 20.0, 10.0;
    auto p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stress);
    KRATOS_EXPECT_NEAR(cap_yield_surface.YieldFunctionValue(p_q), 18.75, Defaults::absolute_tolerance);

    principal_stress <<= 20.0, 15.0, 10.0;
    p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stress);
    KRATOS_EXPECT_NEAR(cap_yield_surface.YieldFunctionValue(p_q), -170.3125, Defaults::absolute_tolerance);

    principal_stress <<= 46.1880215351700611607, 0.0, -46.1880215351700611607;
    p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stress);
    KRATOS_EXPECT_NEAR(cap_yield_surface.YieldFunctionValue(p_q), 0.0, Defaults::absolute_tolerance);

    auto expected_derivative = Vector(2);
    expected_derivative <<= 0.0, 10.0;
    KRATOS_EXPECT_VECTOR_NEAR(cap_yield_surface.DerivativeOfFlowFunction(p_q),
                              expected_derivative, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
