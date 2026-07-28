// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Uxue Chasco
//

#include "custom_constitutive/mohr_coulomb_with_tension_cutoff_elastoplastic_tangent_matrix.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombElastoPlasticTangentMatrixReturnsElasticMatrixWhenDisabled,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix(
        std::make_unique<PlaneStrain>());

    Properties properties;
    properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED");
    properties.SetValue(YOUNG_MODULUS, 1.0e8);
    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(GEO_USE_NUMERICAL_TANGENT_OPERATOR, false);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    Vector stress_vector = ZeroVector(4);
    parameters.SetStressVector(stress_vector);
    Matrix constitutive_matrix = ZeroMatrix(4, 4);
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);
    law.InitializeMaterialResponseCauchy(parameters);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    KRATOS_WATCH(stress_vector);
    KRATOS_WATCH(constitutive_matrix);

    // Assert
    const auto expected_constitutive_matrix = UblasUtilities::CreateMatrix(
        {{1.35E8, 5.77E7, 5.77E7, 0.0},
         {5.77E7, 1.35E8, 5.77E7, 0.0},
         {5.77E7, 5.77E7, 1.35E8, 0.0},
         {0.0, 0.0, 0.0, 3.85E7}});
    KRATOS_EXPECT_MATRIX_NEAR(constitutive_matrix, expected_constitutive_matrix,
                              1.0E6);
}
KRATOS_TEST_CASE_IN_SUITE(MohrCoulombElastoPlasticTangentMatrixReturnsElasticMatrixWhenEnabled,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix(
        std::make_unique<PlaneStrain>());

    Properties properties;
    properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED");
    properties.SetValue(YOUNG_MODULUS, 1.0e8);
    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(GEO_USE_NUMERICAL_TANGENT_OPERATOR, true);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    Vector stress_vector = ZeroVector(4);
    parameters.SetStressVector(stress_vector);
    Matrix constitutive_matrix = ZeroMatrix(4, 4);
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);
    law.InitializeMaterialResponseCauchy(parameters);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    KRATOS_WATCH(stress_vector);
    KRATOS_WATCH(constitutive_matrix);

    // Assert
    const auto expected_constitutive_matrix = UblasUtilities::CreateMatrix(
        {{1.35E8, 5.77E7, 5.77E7, 0.0},
         {5.77E7, 1.35E8, 5.77E7, 0.0},
         {5.77E7, 5.77E7, 1.35E8, 0.0},
         {0.0, 0.0, 0.0, 3.85E7}});
    KRATOS_EXPECT_MATRIX_NEAR(constitutive_matrix, expected_constitutive_matrix,
                              1.0E6);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombElastoPlasticTangentMatrixReturnsPlasticTangentMatrixWhenEnabled,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix(
        std::make_unique<PlaneStrain>());

    Properties properties;
    properties.SetValue(GEO_DRAINAGE_TYPE, "FULLY_COUPLED");
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.25);
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(GEO_USE_NUMERICAL_TANGENT_OPERATOR, true);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    parameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    Vector stress_vector = ZeroVector(4);
    parameters.SetStressVector(stress_vector);
    Matrix constitutive_matrix = ZeroMatrix(4, 4);
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);
    law.InitializeMaterialResponseCauchy(parameters);

    // This strain produces the elastic trial stress {8, 0, -12, 0}.
    strain_vector = UblasUtilities::CreateVector({1.1e-5, 1.0e-6, -1.4e-5, 0.0});

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    KRATOS_WATCH(stress_vector);
    KRATOS_WATCH(constitutive_matrix);

    // Assert
    const auto expected_stress_vector =
        UblasUtilities::CreateVector({7.338673315592010089, 0.0, -11.338673315592010089, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(stress_vector, expected_stress_vector,
                              Defaults::absolute_tolerance);

    int plasticity_status;
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);
    KRATOS_WATCH(plasticity_status);
    KRATOS_EXPECT_EQ(plasticity_status,
                     static_cast<int>(PlasticityStatus::MOHR_COULOMB_FAILURE));

    const auto elastic_constitutive_matrix = UblasUtilities::CreateMatrix(
        {{1.2e6, 4.0e5, 4.0e5, 0.0},
         {4.0e5, 1.2e6, 4.0e5, 0.0},
         {4.0e5, 4.0e5, 1.2e6, 0.0},
         {0.0, 0.0, 0.0, 4.0e5}});
    KRATOS_EXPECT_GT(norm_frobenius(constitutive_matrix - elastic_constitutive_matrix),
                     1.0e3);
}
} // namespace Kratos::Testing
