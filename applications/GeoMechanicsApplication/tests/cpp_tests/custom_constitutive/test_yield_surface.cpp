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
#include "custom_constitutive/principal_stresses.hpp"
#include "custom_constitutive/sigma_tau.hpp"
#include "custom_constitutive/tension_cutoff.h"
#include "custom_utilities/registration_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include "utilities/math_utils.h"

#include <numbers>
#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

class ParametrizedYieldFunctionValuesOfCoulombYieldSurfaceFixture
    : public ::testing::TestWithParam<std::tuple<Geo::PrincipalStresses, double>>
{
};

TEST_P(ParametrizedYieldFunctionValuesOfCoulombYieldSurfaceFixture, CoulombYieldSurface_CalculateYieldFunctionValues)
{
    // Arrange
    auto material_properties                 = Properties{};
    material_properties[GEO_FRICTION_ANGLE]  = 45.0;
    material_properties[GEO_COHESION]        = 2.0;
    material_properties[GEO_DILATANCY_ANGLE] = 0.0;
    const auto coulomb_yield_surface         = CoulombYieldSurface{material_properties};

    const auto& [principal_stresses, expected_value] = GetParam();

    // Act & Assert
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(principal_stresses), expected_value,
                       Defaults::absolute_tolerance);
    const auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stresses);
    KRATOS_EXPECT_NEAR(coulomb_yield_surface.YieldFunctionValue(sigma_tau), expected_value,
                       Defaults::absolute_tolerance);
}

INSTANTIATE_TEST_SUITE_P(
    KratosGeoMechanicsFastSuiteWithoutKernel,
    ParametrizedYieldFunctionValuesOfCoulombYieldSurfaceFixture,
    ::testing::Values(std::make_tuple(Geo::PrincipalStresses{3.0, 2.0, 1.0}, 1.0),
                      std::make_tuple(Geo::PrincipalStresses{1.7071067811865475, 1.0, 0.2928932188134525}, 0.0),
                      std::make_tuple(Geo::PrincipalStresses{0.1715728752538099, -1.0, -1.8284271247461901}, -1.0)));

KRATOS_TEST_CASE_IN_SUITE(CoulombYieldSurface_CanBeSavedAndLoaded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("CoulombYieldSurface"s, CoulombYieldSurface{})};

    auto material_properties                 = Properties{};
    material_properties[GEO_FRICTION_ANGLE]  = 60.0;
    material_properties[GEO_COHESION]        = 2.0;
    material_properties[GEO_DILATANCY_ANGLE] = 30.0;
    const auto coulomb_yield_surface         = CoulombYieldSurface{material_properties};
    auto       serializer                    = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, coulomb_yield_surface);
    auto loaded_coulomb_yield_surface = CoulombYieldSurface{};
    serializer.load("test_tag"s, loaded_coulomb_yield_surface);

    // Assert
    const auto principal_stresses = Geo::PrincipalStresses{1.0, 1.0, 1.0};
    const auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stresses);
    KRATOS_EXPECT_NEAR(loaded_coulomb_yield_surface.YieldFunctionValue(sigma_tau),
                       0.5 * std::sqrt(3.0) - 1, Defaults::absolute_tolerance);
    const auto expected_derivative = UblasUtilities::CreateVector({0.5, 1.0});
    KRATOS_EXPECT_VECTOR_NEAR(loaded_coulomb_yield_surface.DerivativeOfFlowFunction(sigma_tau),
                              expected_derivative, Defaults::absolute_tolerance);
}

class ParametrizedYieldFunctionValuesOfTensionCutOffFixture
    : public ::testing::TestWithParam<std::tuple<Geo::PrincipalStresses, double>>
{
};

TEST_P(ParametrizedYieldFunctionValuesOfTensionCutOffFixture, TensionCutOff_CalculateYieldFunctionValues)
{
    // Arrange
    constexpr auto tensile_strength                  = 2.0;
    const auto     tension_cut_off                   = TensionCutoff{tensile_strength};
    const auto& [principal_stresses, expected_value] = GetParam();

    // Act & Assert
    KRATOS_EXPECT_NEAR(tension_cut_off.YieldFunctionValue(principal_stresses), expected_value,
                       Defaults::absolute_tolerance);
    const auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stresses);
    KRATOS_EXPECT_NEAR(tension_cut_off.YieldFunctionValue(sigma_tau), expected_value, Defaults::absolute_tolerance);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedYieldFunctionValuesOfTensionCutOffFixture,
                         ::testing::Values(std::make_tuple(Geo::PrincipalStresses{3.0, 2.0, 1.0}, 1.0),
                                           std::make_tuple(Geo::PrincipalStresses{2.0, 1.5, 1.0}, 0.0),
                                           std::make_tuple(Geo::PrincipalStresses{1.0, 0.5, 0.1}, -1.0)));

class ParametrizedDerivativeOfFlowFunctionOfTensionCutOff : public ::testing::TestWithParam<Geo::SigmaTau>
{
};

TEST_P(ParametrizedDerivativeOfFlowFunctionOfTensionCutOff, TensionCutOff_DerivativeOfFlowFunctionIsIndependentOfGivenSigmaTau)
{
    // Arrange
    constexpr auto tensile_strength = 2.0;
    const auto     tension_cut_off  = TensionCutoff{tensile_strength};
    const auto&    sigma_tau        = GetParam();
    const auto averaging_type = Geo::PrincipalStresses::PrincipalStressesAveragingType::NO_AVERAGING;

    // Act & Assert
    const auto expected_derivative = UblasUtilities::CreateVector({1.0, 1.0});
    KRATOS_EXPECT_VECTOR_NEAR(tension_cut_off.DerivativeOfFlowFunction(sigma_tau, averaging_type),
                              expected_derivative, Defaults::absolute_tolerance);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedDerivativeOfFlowFunctionOfTensionCutOff,
                         ::testing::Values(Geo::SigmaTau{-1.0, -1.0},
                                           Geo::SigmaTau{0.0, 2.0},
                                           Geo::SigmaTau{1.0, -3.0}));

KRATOS_TEST_CASE_IN_SUITE(TensionCutOff_CanBeSavedAndLoaded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("TensionCutoff"s, TensionCutoff{})};
    constexpr auto tensile_strength = 2.0;
    const auto     tension_cut_off  = TensionCutoff{tensile_strength};
    auto           serializer       = StreamSerializer{};
    const auto averaging_type = Geo::PrincipalStresses::PrincipalStressesAveragingType::NO_AVERAGING;

    // Act
    serializer.save("test_tag"s, tension_cut_off);
    auto loaded_tension_cut_off = TensionCutoff{};
    serializer.load("test_tag"s, loaded_tension_cut_off);

    // Assert
    const auto principal_stresses = Geo::PrincipalStresses{tensile_strength, 0.0, 0.0};
    const auto sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stresses);
    KRATOS_EXPECT_NEAR(loaded_tension_cut_off.YieldFunctionValue(sigma_tau), 0.0, Defaults::absolute_tolerance);
    const auto expected_derivative = UblasUtilities::CreateVector({1.0, 1.0});
    KRATOS_EXPECT_VECTOR_NEAR(loaded_tension_cut_off.DerivativeOfFlowFunction(sigma_tau, averaging_type),
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
    auto material_properties                    = Properties{};
    material_properties[GEO_CAP_HARDENING_TYPE] = "None";
    material_properties[GEO_FRICTION_ANGLE] = (180.0 / std::numbers::pi) * std::asin(24.0 / 25.0);
    material_properties[GEO_PRECONSOLIDATION_STRESS] = 20.0;

    auto principal_stress = UblasUtilities::CreateVector({30.0, 20.0, 10.0});
    auto p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stress);
    KRATOS_EXPECT_NEAR(CompressionCapYieldSurface{material_properties}.YieldFunctionValue(p_q),
                       18.75, Defaults::absolute_tolerance);

    material_properties[K0_NC] = 1.0 / 25.0; // overrule the compression cap size calculated from the friction angle
    principal_stress = UblasUtilities::CreateVector({20.0, 15.0, 10.0});
    p_q              = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stress);
    KRATOS_EXPECT_NEAR(CompressionCapYieldSurface{material_properties}.YieldFunctionValue(p_q),
                       -170.3125, Defaults::absolute_tolerance);

    material_properties[GEO_COMPRESSION_CAP_SIZE] = 4.0; // overrule the compression cap size calculated from K0_NC
    const auto cap_yield_surface = CompressionCapYieldSurface{material_properties};
    principal_stress = UblasUtilities::CreateVector({46.1880215351700611607, 0.0, -46.1880215351700611607});
    p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stress);
    KRATOS_EXPECT_NEAR(cap_yield_surface.YieldFunctionValue(p_q), 0.0, Defaults::absolute_tolerance);

    auto expected_derivative = UblasUtilities::CreateVector({0.0, 10.0});
    KRATOS_EXPECT_VECTOR_NEAR(cap_yield_surface.DerivativeOfFlowFunction(p_q), expected_derivative,
                              Defaults::absolute_tolerance);

    principal_stress = UblasUtilities::CreateVector({51.9615242270663, 0.0, -51.9615242270663});
    p_q              = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stress);
    expected_derivative = UblasUtilities::CreateVector({0.0, 11.25});
    KRATOS_EXPECT_VECTOR_NEAR(cap_yield_surface.DerivativeOfFlowFunction(p_q), expected_derivative,
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CompressionCapYieldSurface_CanBeSavedAndLoadedThroughInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration = ScopedSerializerRegistration{
        std::make_pair("CompressionCapYieldSurface"s, CompressionCapYieldSurface{})};

    auto material_properties                         = Properties{};
    material_properties[GEO_PRECONSOLIDATION_STRESS] = 20.0;
    material_properties[GEO_COMPRESSION_CAP_SIZE]    = 4.0;
    auto p_cap_yield_surface =
        std::unique_ptr<YieldSurface>{std::make_unique<CompressionCapYieldSurface>(material_properties)};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_cap_yield_surface);
    auto p_loaded_cap_yield_surface = std::unique_ptr<YieldSurface>{};
    serializer.load("test_tag"s, p_loaded_cap_yield_surface);

    // Assert
    ASSERT_NE(p_loaded_cap_yield_surface, nullptr);
    auto       principal_stresses = UblasUtilities::CreateVector({30.0, 20.0, 10.0});
    const auto p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stresses);
    KRATOS_EXPECT_NEAR(p_loaded_cap_yield_surface->YieldFunctionValue(p_q), 18.75, Defaults::absolute_tolerance);

    const auto expected_derivative = UblasUtilities::CreateVector({40.0, 2.1650635094610966169});
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_cap_yield_surface->DerivativeOfFlowFunction(p_q),
                              expected_derivative, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CompressionCapYieldSurface_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties material_properties(3);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CompressionCapYieldSurface{material_properties},
        "GEO_PRECONSOLIDATION_STRESS does not exist in the property with Id 3.")
    material_properties.SetValue(GEO_PRECONSOLIDATION_STRESS, -1.0);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CompressionCapYieldSurface{material_properties},
                                      "GEO_PRECONSOLIDATION_STRESS in the property with Id 3 has "
                                      "an invalid value: -1 is out of the range [0, -).")
    material_properties.SetValue(GEO_PRECONSOLIDATION_STRESS, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CompressionCapYieldSurface{material_properties},
        "Invalid material properties: cap size determination requires one of "
        "GEO_COMPRESSION_CAP_SIZE, K0_NC, or GEO_FRICTION_ANGLE to be defined.")
    material_properties.SetValue(GEO_FRICTION_ANGLE, -30.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CompressionCapYieldSurface{material_properties},
                                      "GEO_FRICTION_ANGLE in the property with Id 3 has an invalid "
                                      "value: -30 is out of the range [0, -).")
    material_properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    EXPECT_NO_THROW(CompressionCapYieldSurface{material_properties});
    material_properties.SetValue(K0_NC, -2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CompressionCapYieldSurface{material_properties},
        "K0_NC in the property with Id 3 has an invalid value: -2 is out of the range [0, -).")
    material_properties.SetValue(K0_NC, 2.0);

    EXPECT_NO_THROW(CompressionCapYieldSurface{material_properties});
    material_properties.SetValue(GEO_COMPRESSION_CAP_SIZE, -20.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CompressionCapYieldSurface{material_properties},
                                      "GEO_COMPRESSION_CAP_SIZE in the property with Id 3 has an "
                                      "invalid value: -20 is out of the range [0, -).")
    material_properties.SetValue(GEO_COMPRESSION_CAP_SIZE, 20.0);

    EXPECT_NO_THROW(CompressionCapYieldSurface{material_properties});
}

} // namespace Kratos::Testing
