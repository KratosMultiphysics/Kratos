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
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/mat_variables.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <algorithm>
#include <limits>
#include <numbers>
#include <string>

namespace
{

using namespace Kratos;
using namespace std::string_literals;

GeoIncrementalLinearElasticEurLaw CreateIncrementalLinearElasticEur3DLaw()
{
    return GeoIncrementalLinearElasticEurLaw{std::make_unique<ThreeDimensional>()};
}

double CalculateExpectedNormalDiagonal(double PoissonsRatio, double YoungsModulus)
{
    const auto denominator = (1.0 + PoissonsRatio) * (1.0 - 2.0 * PoissonsRatio);
    return YoungsModulus * (1.0 - PoissonsRatio) / denominator;
}

double CalculateExpectedNormalDiagonal(const Properties& rProperties, const Vector& rStressVector)
{
    const auto reference_pressure = rProperties[GEO_PRESSURE_REFERENCE];
    const auto phi_rad            = rProperties[GEO_FRICTION_ANGLE] * std::numbers::pi / 180.0;
    const auto stress_shift       = rProperties[GEO_COHESION] / std::tan(phi_rad);

    Vector principal_stresses;
    Matrix eigen_vectors;
    StressStrainUtilities::CalculatePrincipalStresses(rStressVector, principal_stresses, eigen_vectors);

    const auto   minor_principal = principal_stresses(2);
    const double base = (stress_shift - minor_principal) / (stress_shift + reference_pressure);
    const auto   expected_youngs_modulus =
        rProperties[YOUNG_MODULUS] * std::pow(base, rProperties[GEO_STRESS_DEPENDENCY_EXPONENT]);

    return CalculateExpectedNormalDiagonal(rProperties[POISSON_RATIO], expected_youngs_modulus);
}

Properties CreateMaterialPropertiesForEurElasticLaw(IndexType Id = 0)
{
    Properties properties(Id);
    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED"s);
    properties.SetValue(GEO_PRESSURE_REFERENCE, 50.0);
    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 1.0);
    properties.SetValue(GEO_COHESION, 1000.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 20.0);
    return properties;
}

Properties CreateConstantYoungsModulusProperties(IndexType Id = 0)
{
    auto properties = CreateMaterialPropertiesForEurElasticLaw(Id);
    properties.SetValue(GEO_PRESSURE_REFERENCE, 1.0e12);
    return properties;
}

Matrix CalculateConstitutiveMatrixForEurElasticLaw(GeoIncrementalLinearElasticEurLaw& rLaw,
                                                   const Properties&                  rProperties,
                                                   Vector&                            rStrainVector)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.SetStrainVector(rStrainVector);
    parameters.SetMaterialProperties(rProperties);

    Matrix constitutive_matrix;
    rLaw.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    return constitutive_matrix;
}

Vector CalculateStressForEurElasticLaw(GeoIncrementalLinearElasticEurLaw& rLaw,
                                       const Properties&                  rProperties,
                                       Vector&                            rStrainVector)
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

double CalculateConstitutiveNormalDiagonalAtZeroStrain(GeoIncrementalLinearElasticEurLaw& rLaw,
                                                       const Properties& rProperties)
{
    auto strain_vector = Vector(6, 0.0);
    return CalculateConstitutiveMatrixForEurElasticLaw(rLaw, rProperties, strain_vector)(0, 0);
}

void InitializeEurLawWithState(GeoIncrementalLinearElasticEurLaw& rLaw, Vector& rStrainVector, Vector& rStressVector)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.SetStrainVector(rStrainVector);
    parameters.SetStressVector(rStressVector);
    rLaw.InitializeMaterialResponseCauchy(parameters);
}

void InitializeEurLawWithFinalizedStress(GeoIncrementalLinearElasticEurLaw& rLaw, Vector& rStressVector)
{
    auto strain_vector = Vector(6, 0.0);
    InitializeEurLawWithState(rLaw, strain_vector, rStressVector);
}

void FinalizeEurLawResponse(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties, Vector& rStrainVector)
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
    rLaw.FinalizeMaterialResponsePK2(parameters);
}

void SetEurLawToIncrementalState(GeoIncrementalLinearElasticEurLaw& rLaw, const Properties& rProperties)
{
    auto strain_vector = Vector(6, 0.5);
    auto stress_vector = Vector(6, 1.0e6);
    InitializeEurLawWithState(rLaw, strain_vector, stress_vector);
    strain_vector = Vector(6, 1.3);
    FinalizeEurLawResponse(rLaw, rProperties, strain_vector);
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
    SetEurLawToIncrementalState(law, properties);
    auto       strain = Vector(6, 1.0);
    const auto stress = CalculateStressForEurElasticLaw(law, properties, strain);

    auto copied_law = GeoIncrementalLinearElasticEurLaw{law};

    const auto empty_properties       = Properties{};
    const auto geometry               = Geometry<Node>{};
    const auto shape_functions_values = Vector{};
    law.ResetMaterial(empty_properties, geometry, shape_functions_values);
    const auto stress_after_reset = CalculateStressForEurElasticLaw(law, properties, strain);

    // Act
    const auto stress_from_copied_law = CalculateStressForEurElasticLaw(copied_law, properties, strain);

    // Assert
    constexpr auto tolerance = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(stress, stress_from_copied_law, tolerance)
    KRATOS_EXPECT_FALSE((std::abs((stress_after_reset[0] - stress[0]) / stress[0]) <= tolerance))
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CopyAssignmentCopiesInternalState,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Act
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetEurLawToIncrementalState(law, properties);
    auto strain = Vector(6, 1.0);
    auto stress = CalculateStressForEurElasticLaw(law, properties, strain);

    auto assigned_law = CreateIncrementalLinearElasticEur3DLaw();
    assigned_law      = law;

    const auto empty_properties       = Properties{};
    const auto geometry               = Geometry<Node>{};
    const auto shape_functions_values = Vector{};
    law.ResetMaterial(empty_properties, geometry, shape_functions_values);
    const auto stress_after_reset = CalculateStressForEurElasticLaw(law, properties, strain);

    // Act
    const auto stress_assigned_law = CalculateStressForEurElasticLaw(assigned_law, properties, strain);

    // Assert
    constexpr auto tolerance = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(stress, stress_assigned_law, tolerance)
    KRATOS_EXPECT_FALSE((std::abs((stress_after_reset[0] - stress[0]) / stress[0]) <= tolerance))
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CloneReturnsCopyOfCorrectType,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Act
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetEurLawToIncrementalState(law, properties);

    // Act and Assert
    const auto p_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_clone.get());

    auto* p_typed_clone = dynamic_cast<GeoIncrementalLinearElasticEurLaw*>(p_clone.get());
    ASSERT_NE(p_typed_clone, nullptr);

    auto strain_vector = Vector(6, 1.0);
    const auto stress_from_clone = CalculateStressForEurElasticLaw(*p_typed_clone, properties, strain_vector);
    // The clone should produce the same stress as the original prior to any resets
    auto expected_stress_from_original = CalculateStressForEurElasticLaw(law, properties, strain_vector);
    constexpr auto tolerance = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_from_original, stress_from_clone, tolerance)
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

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedSpaceAndStrainDimensions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange and Act
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    // Assert
    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 3);
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
        "GEO_PRESSURE_REFERENCE does not exist in the parameters of material with Id 3.")

    properties.SetValue(GEO_PRESSURE_REFERENCE, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "GEO_PRESSURE_REFERENCE in the parameters of material with Id 3 has an "
        "invalid value: 0 is out of the range (0, -).")

    properties.SetValue(GEO_PRESSURE_REFERENCE, 50.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "GEO_STRESS_DEPENDENCY_EXPONENT does not exist in the parameters of material with Id 3.")

    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, geometry, process_info),
        "GEO_STRESS_DEPENDENCY_EXPONENT in the parameters of material with Id 3 has an "
        "invalid value: 0 is out of the range (0, -).")

    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 1.0);
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
    const auto constitutive_matrix = CalculateConstitutiveMatrixForEurElasticLaw(law, properties, strain_vector);

    constexpr double expected_diagonal_value = 0.0369853;
    const auto       expected_constitutive_matrix =
        UblasUtilities::CreateMatrix({{expected_diagonal_value, 0.0, 0.0, 0.0, 0.0, 0.0},
                                      {0.0, expected_diagonal_value, 0.0, 0.0, 0.0, 0.0},
                                      {0.0, 0.0, expected_diagonal_value, 0.0, 0.0, 0.0},
                                      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}});
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(constitutive_matrix, expected_constitutive_matrix,
                                       Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedStressFromPK2Response,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law           = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties    = CreateConstantYoungsModulusProperties();
    auto       strain_vector = Vector(6, 1.0);

    // Act
    const auto calculated_stress = CalculateStressForEurElasticLaw(law, properties, strain_vector);

    // Assert
    constexpr auto tolerance = 1.0e-4;
    const auto     expected_stress =
        UblasUtilities::CreateVector({0.0686869, 0.0686869, 0.0686869, 0.0105672, 0.0105672, 0.0105672});
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, calculated_stress, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedCapabilityFlags,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    // Act and Assert
    KRATOS_EXPECT_TRUE(law.RequiresInitializeMaterialResponse())
    KRATOS_EXPECT_TRUE(law.RequiresFinalizeMaterialResponse())
    KRATOS_EXPECT_TRUE(law.IsIncremental())
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedDiagonalEntryForMultipleScenarios,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Scenario 1: at reference pressure the modulus should equal YOUNG_MODULUS
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateMaterialPropertiesForEurElasticLaw();

    ConstitutiveLaw::Parameters parameters;
    auto                        strain = Vector(6, 0.0);
    parameters.SetStrainVector(strain);
    parameters.SetMaterialProperties(properties);

    const auto reference_pressure = properties[GEO_PRESSURE_REFERENCE];
    auto       stress             = UblasUtilities::CreateVector(
        {-reference_pressure, -reference_pressure, -reference_pressure, 0.0, 0.0, 0.0});
    InitializeEurLawWithFinalizedStress(law, stress);

    Matrix constitutive_matrix;
    law.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    auto expected_value =
        CalculateExpectedNormalDiagonal(properties[POISSON_RATIO], properties[YOUNG_MODULUS]);
    KRATOS_EXPECT_NEAR(constitutive_matrix(0, 0), expected_value, Defaults::relative_tolerance);

    // Scenario 2: diagonal entry scales with confinement level
    law    = CreateIncrementalLinearElasticEur3DLaw();
    stress = UblasUtilities::CreateVector({-500.0, -300.0, -100.0, 0.0, 0.0, 0.0});
    InitializeEurLawWithFinalizedStress(law, stress);

    auto diagonal_entry = CalculateConstitutiveNormalDiagonalAtZeroStrain(law, properties);
    expected_value      = CalculateExpectedNormalDiagonal(properties, stress);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);

    // Scenario 3: low confinement uses the production formula directly
    law    = CreateIncrementalLinearElasticEur3DLaw();
    stress = UblasUtilities::CreateVector({-20.0, -20.0, -20.0, 0.0, 0.0, 0.0});
    InitializeEurLawWithFinalizedStress(law, stress);

    diagonal_entry = CalculateConstitutiveNormalDiagonalAtZeroStrain(law, properties);
    expected_value = CalculateExpectedNormalDiagonal(properties, stress);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);

    // Scenario 4: stress-shift term impact through cohesion and friction angle
    law = CreateIncrementalLinearElasticEur3DLaw();
    properties.SetValue(GEO_COHESION, 20.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 45.0);
    stress = UblasUtilities::CreateVector({-250.0, -150.0, -100.0, 0.0, 0.0, 0.0});
    InitializeEurLawWithFinalizedStress(law, stress);

    diagonal_entry = CalculateConstitutiveNormalDiagonalAtZeroStrain(law, properties);
    expected_value = CalculateExpectedNormalDiagonal(properties, stress);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ThrowsWhenPowBaseIsNonPositive,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: create a state with minor principal stress larger than stress_shift,
    // forcing base <= epsilon in CalculateStressDependentYoungsModulus.
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateMaterialPropertiesForEurElasticLaw();
    properties.SetValue(GEO_COHESION, 1.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 45.0);

    auto finalized_stress = UblasUtilities::CreateVector({100.0, 100.0, 100.0, 0.0, 0.0, 0.0});
    InitializeEurLawWithFinalizedStress(law, finalized_stress);

    auto strain_vector = Vector(6, 1.0);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CalculateStressForEurElasticLaw(law, properties, strain_vector), "Negative base for std::pow (-1.94118). Check GEO_COHESION, GEO_FRICTION_ANGLE, GEO_PRESSURE_REFERENCE and the finalized stress state.")
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_FinalizesMaterialResponseCauchyIncrementally,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();

    auto strain_vector = Vector(6, 0.5);
    auto stress_vector = Vector(6, 1.0e6);
    InitializeEurLawWithState(law, strain_vector, stress_vector);
    strain_vector = Vector(6, 1.3);

    // Act
    FinalizeEurLawResponse(law, properties, strain_vector);

    // Assert: compute expected by repeating the same sequence on a reference law
    strain_vector = Vector(6, 1.0);
    const auto stress_after_finalize = CalculateStressForEurElasticLaw(law, properties, strain_vector);
    const auto expected_stress = UblasUtilities::CreateVector({1e+06, 1e+06, 1e+06, 1e+06, 1e+06, 1e+06});
    constexpr auto tolerance = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress_after_finalize, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ResetMaterialRestoresInitialState,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetEurLawToIncrementalState(law, properties);

    const Properties     empty_properties;
    const Geometry<Node> geometry;
    const Vector         shape_functions_values;

    // Act
    law.ResetMaterial(empty_properties, geometry, shape_functions_values);

    // Assert: expected stress after reset computed from a fresh law in reset state
    auto       strain_vector      = Vector(6, 1.0);
    const auto stress_after_reset = CalculateStressForEurElasticLaw(law, properties, strain_vector);

    auto ref_law = CreateIncrementalLinearElasticEur3DLaw();
    const auto expected_stress = CalculateStressForEurElasticLaw(ref_law, properties, strain_vector);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress_after_reset, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CanBeSavedAndLoaded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetEurLawToIncrementalState(law, properties);
    auto       strain_vector   = Vector(6, 1.0);
    const auto expected_stress = CalculateStressForEurElasticLaw(law, properties, strain_vector);

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

    const auto loaded_stress = CalculateStressForEurElasticLaw(loaded_law, properties, strain_vector);
    constexpr auto tolerance = 1.0e-4;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, loaded_stress, tolerance)
}

} // namespace Kratos::Testing
