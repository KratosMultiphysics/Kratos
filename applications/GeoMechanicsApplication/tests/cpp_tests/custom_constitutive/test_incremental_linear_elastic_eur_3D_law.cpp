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

#include "custom_constitutive/incremental_linear_elastic_eur_law.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/mat_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <cmath>

namespace
{

using namespace Kratos;

GeoIncrementalLinearElasticEurLaw CreateIncrementalLinearElasticEur3DLaw()
{
    return GeoIncrementalLinearElasticEurLaw{std::make_unique<ThreeDimensional>()};
}

double CalculateExpectedNormalDiagonal(double YoungsModulus, double PoissonsRatio)
{
    const auto denominator = (1.0 + PoissonsRatio) * (1.0 - 2.0 * PoissonsRatio);
    return YoungsModulus * (1.0 - PoissonsRatio) / denominator;
}

Properties CreateValidMaterialProperties()
{
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED");
    properties.SetValue(REFERENCE_HARDENING_MODULUS, 50.0);
    properties.SetValue(SWELLING_SLOPE, 1.0);
    return properties;
}

double CalculateConstitutiveNormalDiagonal(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties)
{
    ConstitutiveLaw::Parameters parameters;
    auto                        strain = Vector{ScalarVector{6, 0.0}};
    parameters.SetStrainVector(strain);
    parameters.SetMaterialProperties(rProperties);

    Matrix constitutive_matrix;
    rLaw.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    return constitutive_matrix(0, 0);
}

void InitializeLawWithFinalizedStress(GeoIncrementalLinearElasticEurLaw& rLaw, Vector& rStressVector)
{
    ConstitutiveLaw::Parameters initialize_parameters;
    auto                        initial_strain = Vector(6, 0.0);
    initialize_parameters.SetStrainVector(initial_strain);
    initialize_parameters.SetStressVector(rStressVector);
    rLaw.InitializeMaterialResponseCauchy(initialize_parameters);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLawReturnsExpectedDiagonalEntryAtReferencePressure,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    ConstitutiveLaw::Parameters parameters;
    auto                        strain = Vector{ScalarVector{6, 0.0}};
    parameters.SetStrainVector(strain);

    auto properties = CreateValidMaterialProperties();
    parameters.SetMaterialProperties(properties);

    // Act
    Matrix constitutive_matrix;
    law.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    // Assert
    constexpr auto youngs_modulus = 1.0e7;
    constexpr auto poisson_ratio  = 0.3;
    const auto     expected_value = CalculateExpectedNormalDiagonal(youngs_modulus, poisson_ratio);
    KRATOS_EXPECT_NEAR(constitutive_matrix(0, 0), expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLawScalesDiagonalEntryWithConfinement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    auto properties = CreateValidMaterialProperties();

    auto initial_stress = UblasUtilities::CreateVector({-500.0, -300.0, -100.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, initial_stress);

    // Act
    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    // Assert
    const auto eur_ref    = properties[YOUNG_MODULUS];
    const auto expected_E = eur_ref * 2.0; // (s - p)/(s + p_ref) = (0 - (-100))/50 = 2 and m = 1
    const auto expected_value = CalculateExpectedNormalDiagonal(expected_E, properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLawUsesReferencePressureAtLowConfinement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();

    // sigma3' = -20, so -sigma3' < p_ref (=50) and the lower bound should enforce p_ref.
    auto low_confinement_stress = UblasUtilities::CreateVector({-20.0, -20.0, -20.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, low_confinement_stress);

    // Act
    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    // Assert
    const auto expected_value =
        CalculateExpectedNormalDiagonal(properties[YOUNG_MODULUS], properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLawAccountsForStressShiftTerm,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();
    properties.SetValue(GEO_COHESION, 20.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 45.0);

    // sigma3' = -100 -> bounded minor principal effective stress remains -100.
    auto stress_state = UblasUtilities::CreateVector({-250.0, -150.0, -100.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, stress_state);

    // Act
    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    // Assert
    constexpr auto pi           = 3.14159265358979323846;
    const auto     phi_rad      = properties[GEO_FRICTION_ANGLE] * pi / 180.0;
    const auto     stress_shift = properties[GEO_COHESION] / std::tan(phi_rad);
    const auto     expected_E   = properties[YOUNG_MODULUS] * (stress_shift - (-100.0)) /
                            (stress_shift + properties[REFERENCE_HARDENING_MODULUS]);
    const auto expected_value = CalculateExpectedNormalDiagonal(expected_E, properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

} // namespace Kratos::Testing
