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

    auto       copied_law    = GeoIncrementalLinearElasticEurLaw{law};
    auto       strain_vector = Vector(6, 1.0);
    const auto copied_stress = CalculateStress(copied_law, properties, strain_vector);

    const Properties     empty_properties;
    const Geometry<Node> geometry;
    const Vector         shape_functions_values;
    law.ResetMaterial(empty_properties, geometry, shape_functions_values);

    strain_vector                          = Vector(6, 1.0);
    const auto original_stress_after_reset = CalculateStress(law, properties, strain_vector);
    const auto expected_copied_stress =
        UblasUtilities::CreateVector({1.35e7, 1.35e7, 1.35e7, 2.92308e6, 2.92308e6, 2.92308e6});
    const auto expected_reset_stress =
        UblasUtilities::CreateVector({2.5e7, 2.5e7, 2.5e7, 3.84615e6, 3.84615e6, 3.84615e6});

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_copied_stress, copied_stress, 1.0e-3)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_reset_stress, original_stress_after_reset, 1.0e-3)
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
    law.ResetMaterial(empty_properties, geometry, shape_functions_values);

    strain_vector                          = Vector(6, 1.0);
    const auto original_stress_after_reset = CalculateStress(law, properties, strain_vector);
    const auto expected_assigned_stress =
        UblasUtilities::CreateVector({1.35e7, 1.35e7, 1.35e7, 2.92308e6, 2.92308e6, 2.92308e6});
    const auto expected_reset_stress =
        UblasUtilities::CreateVector({2.5e7, 2.5e7, 2.5e7, 3.84615e6, 3.84615e6, 3.84615e6});

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_assigned_stress, assigned_stress, 1.0e-3)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_reset_stress, original_stress_after_reset, 1.0e-3)
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
    const auto expected_stress =
        UblasUtilities::CreateVector({1.35e7, 1.35e7, 1.35e7, 2.92308e6, 2.92308e6, 2.92308e6});
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, clone_stress, 1.0e-3)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsTrueForStenbergShearStabilizationSuitability,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law         = CreateIncrementalLinearElasticEur3DLaw();
    auto is_suitable = false;

    auto& r_value = law.GetValue(STENBERG_SHEAR_STABILIZATION_SUITABLE, is_suitable);

    KRATOS_EXPECT_EQ(&r_value, &is_suitable);
    KRATOS_EXPECT_TRUE(is_suitable);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedLawFeatures,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    ConstitutiveLaw::Features law_features;
    law.GetLawFeatures(law_features);

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
    auto law = CreateIncrementalLinearElasticEur3DLaw();
    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedStrainSize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law = CreateIncrementalLinearElasticEur3DLaw();
    KRATOS_EXPECT_EQ(law.GetStrainSize(), 6);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ChecksAdditionalMaterialParameters,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law          = CreateIncrementalLinearElasticEur3DLaw();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    auto valid_properties = CreateValidMaterialProperties(3);
    KRATOS_EXPECT_EQ(law.Check(valid_properties, geometry, process_info), 0);

    Properties missing_reference_pressure_properties(3);
    missing_reference_pressure_properties.SetValue(YOUNG_MODULUS, 1.0e7);
    missing_reference_pressure_properties.SetValue(POISSON_RATIO, 0.3);
    missing_reference_pressure_properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED"s);
    missing_reference_pressure_properties.SetValue(SWELLING_SLOPE, 1.0);
    EXPECT_THROW(law.Check(missing_reference_pressure_properties, geometry, process_info), Exception);

    auto invalid_swelling_slope_properties = CreateValidMaterialProperties(3);
    invalid_swelling_slope_properties.SetValue(SWELLING_SLOPE, 0.0);
    EXPECT_THROW(law.Check(invalid_swelling_slope_properties, geometry, process_info), Exception);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsDiagonalConstitutiveMatrixWhenOnlyDiagonalEntriesAreConsidered,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElasticEur3DLaw();
    law.SetConsiderDiagonalEntriesOnlyAndNoShear(true);

    const auto properties          = CreateConstantYoungsModulusProperties();
    auto       strain_vector       = Vector(6, 0.0);
    const auto constitutive_matrix = CalculateConstitutiveMatrix(law, properties, strain_vector);
    const auto expected_normal_diagonal =
        CalculateExpectedNormalDiagonal(properties[YOUNG_MODULUS], properties[POISSON_RATIO]);

    for (IndexType i = 0; i < 3; ++i) {
        KRATOS_EXPECT_NEAR(constitutive_matrix(i, i), expected_normal_diagonal, Defaults::relative_tolerance);
    }

    for (IndexType i = 0; i < 3; ++i) {
        for (IndexType j = 0; j < 3; ++j) {
            if (i != j) {
                KRATOS_EXPECT_NEAR(constitutive_matrix(i, j), 0.0, Defaults::absolute_tolerance);
            }
        }
    }

    for (IndexType i = 3; i < 6; ++i) {
        KRATOS_EXPECT_NEAR(constitutive_matrix(i, i), 0.0, Defaults::absolute_tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedStressFromPK2Response,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();

    auto       strain_vector = Vector(6, 1.0);
    const auto stress        = CalculateStress(law, properties, strain_vector);
    const auto expected_stress =
        UblasUtilities::CreateVector({2.5e7, 2.5e7, 2.5e7, 3.84615e6, 3.84615e6, 3.84615e6});
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1.0e-3)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_RequiresInitializeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElasticEur3DLaw();
    KRATOS_EXPECT_TRUE(law.RequiresInitializeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_RequiresFinalizeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElasticEur3DLaw();
    KRATOS_EXPECT_TRUE(law.RequiresFinalizeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ReturnsExpectedDiagonalEntryAtReferencePressure,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElasticEur3DLaw();

    ConstitutiveLaw::Parameters parameters;
    auto                        strain = Vector{ScalarVector{6, 0.0}};
    parameters.SetStrainVector(strain);

    auto properties = CreateValidMaterialProperties();
    parameters.SetMaterialProperties(properties);

    Matrix constitutive_matrix;
    law.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    constexpr auto youngs_modulus = 1.0e7;
    constexpr auto poisson_ratio  = 0.3;
    const auto     expected_value = CalculateExpectedNormalDiagonal(youngs_modulus, poisson_ratio);
    KRATOS_EXPECT_NEAR(constitutive_matrix(0, 0), expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ScalesDiagonalEntryWithConfinement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();

    auto initial_stress = UblasUtilities::CreateVector({-500.0, -300.0, -100.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, initial_stress);

    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    const auto eur_ref    = properties[YOUNG_MODULUS];
    const auto expected_E = eur_ref * 2.0;
    const auto expected_value = CalculateExpectedNormalDiagonal(expected_E, properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_UsesReferencePressureAtLowConfinement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();

    auto low_confinement_stress = UblasUtilities::CreateVector({-20.0, -20.0, -20.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, low_confinement_stress);

    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);
    const auto expected_value =
        CalculateExpectedNormalDiagonal(properties[YOUNG_MODULUS], properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_AccountsForStressShiftTerm,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law        = CreateIncrementalLinearElasticEur3DLaw();
    auto properties = CreateValidMaterialProperties();
    properties.SetValue(GEO_COHESION, 20.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 45.0);

    auto stress_state = UblasUtilities::CreateVector({-250.0, -150.0, -100.0, 0.0, 0.0, 0.0});
    InitializeLawWithFinalizedStress(law, stress_state);

    const auto diagonal_entry = CalculateConstitutiveNormalDiagonal(law, properties);

    const auto phi_rad      = properties[GEO_FRICTION_ANGLE] * std::numbers::pi / 180.0;
    const auto stress_shift = properties[GEO_COHESION] / std::tan(phi_rad);
    const auto expected_E   = properties[YOUNG_MODULUS] * (stress_shift - (-100.0)) /
                            (stress_shift + properties[REFERENCE_HARDENING_MODULUS]);
    const auto expected_value = CalculateExpectedNormalDiagonal(expected_E, properties[POISSON_RATIO]);
    KRATOS_EXPECT_NEAR(diagonal_entry, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_FinalizesMaterialResponseCauchyIncrementally,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();

    auto strain_vector = Vector(6, 0.5);
    auto stress_vector = Vector(6, 1.0e6);
    InitializeLawWithState(law, strain_vector, stress_vector);
    strain_vector = Vector(6, 1.3);
    FinalizeLawResponse(law, properties, strain_vector);

    strain_vector                    = Vector(6, 1.0);
    const auto stress_after_finalize = CalculateStress(law, properties, strain_vector);
    const auto expected_stress =
        UblasUtilities::CreateVector({1.35e7, 1.35e7, 1.35e7, 2.92308e6, 2.92308e6, 2.92308e6});
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress_after_finalize, 1.0e-3)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_ResetMaterialRestoresInitialState,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetLawToIncrementalState(law, properties);

    const Properties     empty_properties;
    const Geometry<Node> geometry;
    const Vector         shape_functions_values;
    law.ResetMaterial(empty_properties, geometry, shape_functions_values);

    auto       strain_vector      = Vector(6, 1.0);
    const auto stress_after_reset = CalculateStress(law, properties, strain_vector);
    const auto expected_stress =
        UblasUtilities::CreateVector({2.5e7, 2.5e7, 2.5e7, 3.84615e6, 3.84615e6, 3.84615e6});
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress_after_reset, 1.0e-3)
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElasticEur3DLaw_CanBeSavedAndLoaded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law        = CreateIncrementalLinearElasticEur3DLaw();
    const auto properties = CreateConstantYoungsModulusProperties();
    SetLawToIncrementalState(law, properties);

    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("ThreeDimensional"s, ThreeDimensional{})};
    auto serializer = StreamSerializer{};

    serializer.save("test_tag"s, law);

    auto loaded_law = GeoIncrementalLinearElasticEurLaw{};
    serializer.load("test_tag"s, loaded_law);

    KRATOS_EXPECT_EQ(loaded_law.WorkingSpaceDimension(), 3);
    KRATOS_EXPECT_EQ(loaded_law.GetStrainSize(), 6);

    auto       strain_vector = Vector(6, 1.0);
    const auto loaded_stress = CalculateStress(loaded_law, properties, strain_vector);
    const auto expected_stress =
        UblasUtilities::CreateVector({1.35e7, 1.35e7, 1.35e7, 2.92308e6, 2.92308e6, 2.92308e6});
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, loaded_stress, 1.0e-3)
}

} // namespace Kratos::Testing
