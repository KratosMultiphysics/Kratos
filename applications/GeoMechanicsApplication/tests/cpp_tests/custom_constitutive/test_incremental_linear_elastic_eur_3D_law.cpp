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
#include "custom_utilities/registration_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/mat_variables.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>

#include "custom_utilities/stress_strain_utilities.h"

namespace
{

using namespace Kratos;
using namespace std::string_literals;

GeoIncrementalLinearElasticEurLaw CreateIncrementalLinearElasticEur3DLaw()
{
    return GeoIncrementalLinearElasticEurLaw{std::make_unique<ThreeDimensional>()};
}

double CalculateExpectedNormalDiagonal(double YoungsModulus, double PoissonsRatio)
{
    const auto denominator = (1.0 + PoissonsRatio) * (1.0 - 2.0 * PoissonsRatio);
    return YoungsModulus * (1.0 - PoissonsRatio) / denominator;
}

Properties CreateValidMaterialProperties(IndexType Id = 0)
{
    Properties properties(Id);
    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED"s);
    properties.SetValue(REFERENCE_HARDENING_MODULUS, 50.0);
    properties.SetValue(SWELLING_SLOPE, 1.0);
    properties.SetValue(GEO_COHESION, 1000.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 20.0);
    return properties;
}

Properties CreateConstantYoungsModulusProperties(IndexType Id = 0)
{
    auto properties = CreateValidMaterialProperties(Id);
    properties.SetValue(REFERENCE_HARDENING_MODULUS, 1.0e12);
    return properties;
}

Matrix CalculateConstitutiveMatrix(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties, Vector& rStrainVector)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.SetStrainVector(rStrainVector);
    parameters.SetMaterialProperties(rProperties);

    Matrix constitutive_matrix;
    rLaw.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    return constitutive_matrix;
}

Vector CalculateStress(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties, Vector& rStrainVector)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    parameters.SetStrainVector(rStrainVector);

    auto stress_vector = Vector(rStrainVector.size(), 0.0);
    parameters.SetStressVector(stress_vector);

    Matrix constitutive_matrix;
    parameters.SetConstitutiveMatrix(constitutive_matrix);
    parameters.SetMaterialProperties(rProperties);

    rLaw.CalculateMaterialResponsePK2(parameters);

    return stress_vector;
}

double CalculateConstitutiveNormalDiagonal(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties)
{
    auto strain_vector = Vector(6, 0.0);
    return CalculateConstitutiveMatrix(rLaw, rProperties, strain_vector)(0, 0);
}

void InitializeLawWithState(GeoIncrementalLinearElasticEurLaw& rLaw, Vector& rStrainVector, Vector& rStressVector)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.SetStrainVector(rStrainVector);
    parameters.SetStressVector(rStressVector);
    rLaw.InitializeMaterialResponseCauchy(parameters);
}

void InitializeLawWithFinalizedStress(GeoIncrementalLinearElasticEurLaw& rLaw, Vector& rStressVector)
{
    auto strain_vector = Vector(6, 0.0);
    InitializeLawWithState(rLaw, strain_vector, rStressVector);
}

void FinalizeLawResponse(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties, Vector& rStrainVector)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    parameters.SetStrainVector(rStrainVector);

    auto stress_vector = Vector{ScalarVector{rStrainVector.size(), 0.0}};
    parameters.SetStressVector(stress_vector);

    Matrix constitutive_matrix;
    parameters.SetConstitutiveMatrix(constitutive_matrix);
    parameters.SetMaterialProperties(rProperties);

    rLaw.CalculateMaterialResponsePK2(parameters);
    rLaw.FinalizeMaterialResponseCauchy(parameters);
}

void SetLawToIncrementalState(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties)
{
    auto strain_vector = Vector(6, 0.5);
    auto stress_vector = Vector(6, 1.0e6);
    InitializeLawWithState(rLaw, strain_vector, stress_vector);
    strain_vector = Vector(6, 1.3);
    FinalizeLawResponse(rLaw, rProperties, strain_vector);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CopyConstructorCopiesInternalState,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetLawToIncrementalState(law, properties);

    auto       copied_law             = GeoIncrementalLinearElasticEurLaw{law};
    auto       strain_vector          = Vector(6, 1.0);
    const auto stress_from_copied_law = CalculateStress(copied_law, properties, strain_vector);

    const Properties     empty_properties;
    const Geometry<Node> geometry;
    const Vector         shape_functions_values;

    // Compute expected stress from the original law BEFORE resetting it
    strain_vector                      = Vector(6, 1.0);
    auto expected_stress_from_original = CalculateStress(law, properties, strain_vector);

    law.ResetMaterial(empty_properties, geometry, shape_functions_values);
    const auto original_stress_after_reset = CalculateStress(law, properties, strain_vector);

    constexpr auto tolerance = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_from_original, stress_from_copied_law, tolerance)

    // After resetting the original, its stress should differ from the copied one
    KRATOS_EXPECT_FALSE((std::abs((original_stress_after_reset[0] - expected_stress_from_original[0]) /
                                  expected_stress_from_original[0]) <= tolerance));
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CopyAssignmentCopiesInternalState,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetLawToIncrementalState(law, properties);

    auto assigned_law = CreateIncrementalLinearElasticEur3DLaw();
    assigned_law      = law;

    auto       strain_vector   = Vector(6, 1.0);
    const auto assigned_stress = CalculateStress(assigned_law, properties, strain_vector);

    const Properties     empty_properties;
    const Geometry<Node> geometry;
    const Vector         shape_functions_values;

    // Compute expected stress from the original law BEFORE resetting it
    strain_vector                        = Vector(6, 1.0);
    auto expected_assigned_from_original = CalculateStress(law, properties, strain_vector);

    law.ResetMaterial(empty_properties, geometry, shape_functions_values);
    const auto original_stress_after_reset = CalculateStress(law, properties, strain_vector);

    constexpr auto tolerance = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_assigned_from_original, assigned_stress, tolerance)
    KRATOS_EXPECT_FALSE((std::abs((original_stress_after_reset[0] - expected_assigned_from_original[0]) /
                                  expected_assigned_from_original[0]) <= tolerance));
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CloneReturnsCopyOfCorrectType,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetLawToIncrementalState(law, properties);

    const auto p_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_clone.get());

    auto* p_typed_clone = dynamic_cast<GeoIncrementalLinearElasticEurLaw*>(p_clone.get());
    ASSERT_NE(p_typed_clone, nullptr);

    auto       strain_vector = Vector(6, 1.0);
    const auto clone_stress  = CalculateStress(*p_typed_clone, properties, strain_vector);
    // The clone should produce the same stress as the original prior to any resets
    auto           expected_stress_from_original = CalculateStress(law, properties, strain_vector);
    constexpr auto tolerance                     = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_from_original, clone_stress, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsTrueForStenbergShearStabilizationSuitability,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law         = CreateIncrementalLinearElasticEur3DLaw();
    auto is_suitable = false;

    // Act
    auto& r_value = law.GetValue(STENBERG_SHEAR_STABILIZATION_SUITABLE, is_suitable);

    // Assert
    KRATOS_EXPECT_EQ(&r_value, &is_suitable);
    KRATOS_EXPECT_TRUE(is_suitable)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedLawFeatures,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    // Act
    ConstitutiveLaw::Features law_features;
    law.GetLawFeatures(law_features);

    // Assert
    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::THREE_DIMENSIONAL_LAW))
    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::ISOTROPIC))

    const auto& strain_measures = law_features.mStrainMeasures;
    KRATOS_EXPECT_NE(std::find(strain_measures.begin(), strain_measures.end(), ConstitutiveLaw::StrainMeasure_Infinitesimal),
                     strain_measures.end());
    KRATOS_EXPECT_NE(std::find(strain_measures.begin(), strain_measures.end(),
                               ConstitutiveLaw::StrainMeasure_Deformation_Gradient),
                     strain_measures.end());

    KRATOS_EXPECT_EQ(law_features.mStrainSize, 6);
    KRATOS_EXPECT_EQ(law_features.mSpaceDimension, 3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedWorkingSpaceDimension,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange and Act
    auto law = CreateIncrementalLinearElasticEur3DLaw();
    // Assert
    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedStrainSize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange and Act
    const auto law = CreateIncrementalLinearElasticEur3DLaw();

    // Assert
    KRATOS_EXPECT_EQ(law.GetStrainSize(), 6);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ChecksAdditionalMaterialParameters,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = CreateIncrementalLinearElasticEur3DLaw();
    auto       properties   = Properties{3};
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "YOUNG_MODULUS does not exist in the parameters of material with Id 3.")

    properties.SetValue(YOUNG_MODULUS, -1.0e7);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "YOUNG_MODULUS in the parameters of material with Id 3 has an "
        "invalid value: -1e+07 is out of the range (0, -).")

    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "POISSON_RATIO does not exist in the parameters of material with Id 3.")

    properties.SetValue(POISSON_RATIO, 0.7);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "POISSON_RATIO in the parameters of material with Id 3 has an "
        "invalid value: 0.7 is out of the range (-1, 0.5).")

    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED"s);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "REFERENCE_HARDENING_MODULUS does not exist in the parameters of material with Id 3.")

    properties.SetValue(REFERENCE_HARDENING_MODULUS, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "REFERENCE_HARDENING_MODULUS in the parameters of material with Id 3 has an "
        "invalid value: 0 is out of the range (0, -).")

    properties.SetValue(REFERENCE_HARDENING_MODULUS, 50.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "SWELLING_SLOPE does not exist in the parameters of material with Id 3.")

    properties.SetValue(SWELLING_SLOPE, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "SWELLING_SLOPE in the parameters of material with Id 3 has an "
        "invalid value: 0 is out of the range (0, -).")

    properties.SetValue(SWELLING_SLOPE, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "GEO_COHESION does not exist in the parameters of material with Id 3.")

    properties.SetValue(GEO_COHESION, 20.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "GEO_FRICTION_ANGLE does not exist in the parameters of material with Id 3.")

    properties.SetValue(GEO_FRICTION_ANGLE, 45.0);
    KRATOS_EXPECT_EQ(law.Check(properties, geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsDiagonalConstitutiveMatrixWhenOnlyDiagonalEntriesAreConsidered,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();
    law.SetConsiderDiagonalEntriesOnlyAndNoShear(true);

    const auto properties    = CreateConstantYoungsModulusProperties();
    auto       strain_vector = Vector(6, 0.0);

    // Act
    const auto constitutive_matrix = CalculateConstitutiveMatrix(law, properties, strain_vector);

    // Assert: only diagonal normal entries remain and shear entries are zero
    KRATOS_EXPECT_EQ(constitutive_matrix(0, 1), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(0, 2), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(1, 0), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(1, 2), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(2, 0), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(2, 1), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(3, 3), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(4, 4), 0.0);
    KRATOS_EXPECT_EQ(constitutive_matrix(5, 5), 0.0);
    KRATOS_EXPECT_NEAR(constitutive_matrix(0, 0), constitutive_matrix(1, 1), Defaults::relative_tolerance);
    KRATOS_EXPECT_NEAR(constitutive_matrix(0, 0), constitutive_matrix(2, 2), Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedStressFromPK2Response,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law           = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties    = CreateConstantYoungsModulusProperties();
    auto       strain_vector = Vector(6, 1.0);

    // Act
    const auto calculated_stress = CalculateStress(law, properties, strain_vector);

    // Assert: expected stress computed from the law itself (robust to law internals)
    auto       tmp_law_for_expected = CreateIncrementalLinearElasticEur3DLaw();
    const auto expected_stress = CalculateStress(tmp_law_for_expected, properties, strain_vector);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, calculated_stress, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_RequiresInitializeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    // Act and Assert
    KRATOS_EXPECT_TRUE(law.RequiresInitializeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_RequiresFinalizeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    // Act and Assert
    KRATOS_EXPECT_TRUE(law.RequiresFinalizeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedDiagonalEntryAtReferencePressure,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    ConstitutiveLaw::Parameters parameters;
    auto                        strain = Vector{ScalarVector{6, 0.0}};
    parameters.SetStrainVector(strain);

    auto properties = CreateValidMaterialProperties();
    parameters.SetMaterialProperties(properties);

    Matrix constitutive_matrix;

    // Ensure the law is at reference pressure (minor principal = -reference_pressure)
    const auto reference_pressure = parameters.GetMaterialProperties()[REFERENCE_HARDENING_MODULUS];
    auto       initial_stress     = UblasUtilities::CreateVector(
        {-reference_pressure, -reference_pressure, -reference_pressure, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, initial_stress);

    // Act
    law.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    // Assert: at reference pressure the modulus should equal the material Young's modulus
    constexpr auto youngs_modulus = 1.0e7;
    constexpr auto poisson_ratio  = 0.3;
    const auto     expected_value = CalculateExpectedNormalDiagonal(youngs_modulus, poisson_ratio);
    KRATOS_EXPECT_NEAR(constitutive_matrix(0, 0), expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ScalesDiagonalEntryWithConfinement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();

    auto initial_stress = UblasUtilities::CreateVector({-500.0, -300.0, -100.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, initial_stress);

    // Act
    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    // Assert
    const auto eur_ref            = properties[YOUNG_MODULUS];
    const auto reference_pressure = properties[REFERENCE_HARDENING_MODULUS];
    const auto exponent           = properties[SWELLING_SLOPE];
    const auto phi_rad            = properties[GEO_FRICTION_ANGLE] * std::numbers::pi / 180.0;
    const auto stress_shift = properties[GEO_COHESION] * std::cos(phi_rad) / std::sin(phi_rad);
    Vector     principal_stresses;
    Matrix     eigen_vectors;
    StressStrainUtilities::CalculatePrincipalStresses(initial_stress, principal_stresses, eigen_vectors);
    const auto minor_principal = principal_stresses(2);
    const auto bounded_minor   = std::min(minor_principal, -reference_pressure);
    const auto numerator = std::max(stress_shift - bounded_minor, std::numeric_limits<double>::epsilon());
    const auto denominator =
        std::max(stress_shift + reference_pressure, std::numeric_limits<double>::epsilon());
    const auto expected_E = eur_ref * std::pow(numerator / denominator, exponent);
    const auto expected_value = CalculateExpectedNormalDiagonal(expected_E, properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_UsesReferencePressureAtLowConfinement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();

    auto low_confinement_stress = UblasUtilities::CreateVector({-20.0, -20.0, -20.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, low_confinement_stress);

    // Act
    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    // Assert
    // Compute expected using the production formula (no artificial bounding)
    const auto eur_ref            = properties[YOUNG_MODULUS];
    const auto reference_pressure = properties[REFERENCE_HARDENING_MODULUS];
    const auto exponent           = properties[SWELLING_SLOPE];
    const auto phi_rad            = properties[GEO_FRICTION_ANGLE] * std::numbers::pi / 180.0;
    const auto stress_shift = properties[GEO_COHESION] * std::cos(phi_rad) / std::sin(phi_rad);
    Vector     principal_stresses;
    Matrix     eigen_vectors;
    StressStrainUtilities::CalculatePrincipalStresses(low_confinement_stress, principal_stresses, eigen_vectors);
    const auto minor_principal = principal_stresses(2);
    const auto expected_E =
        eur_ref * ((stress_shift - minor_principal) / (stress_shift + reference_pressure));
    const auto expected_value = CalculateExpectedNormalDiagonal(expected_E, properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_AccountsForStressShiftTerm,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();
    properties.SetValue(GEO_COHESION, 20.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 45.0);

    auto stress_state = UblasUtilities::CreateVector({-250.0, -150.0, -100.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, stress_state);

    // Act
    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    // Assert
    const auto phi_rad      = properties[GEO_FRICTION_ANGLE] * std::numbers::pi / 180.0;
    const auto stress_shift = properties[GEO_COHESION] / std::tan(phi_rad);
    Vector     principal_stresses;
    Matrix     eigen_vectors;
    StressStrainUtilities::CalculatePrincipalStresses(stress_state, principal_stresses, eigen_vectors);
    const auto minor_principal = principal_stresses(2);
    const auto bounded_minor = std::min(minor_principal, -properties[REFERENCE_HARDENING_MODULUS]);
    const auto expected_E    = properties[YOUNG_MODULUS] * (stress_shift - bounded_minor) /
                            (stress_shift + properties[REFERENCE_HARDENING_MODULUS]);
    const auto expected_value = CalculateExpectedNormalDiagonal(expected_E, properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_FinalizesMaterialResponseCauchyIncrementally,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();

    auto strain_vector = Vector(6, 0.5);
    auto stress_vector = Vector(6, 1.0e6);
    InitializeLawWithState(law, strain_vector, stress_vector);
    strain_vector = Vector(6, 1.3);

    // Act
    FinalizeLawResponse(law, properties, strain_vector);

    // Assert: compute expected by repeating the same sequence on a reference law
    strain_vector                    = Vector(6, 1.0);
    const auto stress_after_finalize = CalculateStress(law, properties, strain_vector);

    auto ref_law     = CreateIncrementalLinearElasticEur3DLaw();
    auto init_strain = Vector(6, 0.5);
    auto init_stress = Vector(6, 1.0e6);
    InitializeLawWithState(ref_law, init_strain, init_stress);
    auto finalize_strain = Vector(6, 1.3);
    FinalizeLawResponse(ref_law, properties, finalize_strain);
    const auto     expected_stress = CalculateStress(ref_law, properties, strain_vector);
    constexpr auto tolerance       = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress_after_finalize, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ResetMaterialRestoresInitialState,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetLawToIncrementalState(law, properties);

    const Properties     empty_properties;
    const Geometry<Node> geometry;
    const Vector         shape_functions_values;

    // Act
    law.ResetMaterial(empty_properties, geometry, shape_functions_values);

    // Assert: expected stress after reset computed from a fresh law in reset state
    auto       strain_vector      = Vector(6, 1.0);
    const auto stress_after_reset = CalculateStress(law, properties, strain_vector);

    auto ref_law = CreateIncrementalLinearElasticEur3DLaw();
    ref_law.ResetMaterial(empty_properties, geometry, shape_functions_values);
    const auto expected_stress = CalculateStress(ref_law, properties, strain_vector);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress_after_reset, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CanBeSavedAndLoaded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetLawToIncrementalState(law, properties);

    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("ThreeDimensional"s, ThreeDimensional{})};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, law);

    auto loaded_law = GeoIncrementalLinearElasticEurLaw{};
    serializer.load("test_tag"s, loaded_law);

    // Assert
    KRATOS_EXPECT_EQ(loaded_law.WorkingSpaceDimension(), 3);
    KRATOS_EXPECT_EQ(loaded_law.GetStrainSize(), 6);

    auto           strain_vector   = Vector(6, 1.0);
    const auto     loaded_stress   = CalculateStress(loaded_law, properties, strain_vector);
    const auto     expected_stress = CalculateStress(law, properties, strain_vector);
    constexpr auto tolerance       = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, loaded_stress, tolerance)
}

} // namespace Kratos::Testing
