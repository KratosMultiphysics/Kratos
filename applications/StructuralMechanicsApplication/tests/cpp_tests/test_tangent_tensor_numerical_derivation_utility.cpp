// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "utilities/math_utils.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

namespace Kratos
{
namespace Testing
{

void SettingBasicCase(
    ConstitutiveLaw::Parameters& rCLConfigurationValues,
    Properties& rProperties,
    Vector& rStressVector,
    Vector& rStrainVector,
    Matrix& rTangentModuli,
    Matrix& rF,
    const bool ElementProvidedStrain = true
    )
{
    ProcessInfo current_process_info;

    Flags constitutive_law_options;
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, ElementProvidedStrain);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    rProperties.SetValue(YOUNG_MODULUS, 210e9);
    rProperties.SetValue(POISSON_RATIO, 0.22);

    rStressVector = ZeroVector(6);
    rStressVector[0] = 5.40984e+06;
    rStressVector[1] = 5.40984e+06;
    rStressVector[2] = 1.91803e+07;
    rStressVector[5] = 1.45804e-10;

    rStrainVector = ZeroVector(6);
    rStrainVector[2] = 8.0e-5;
    rStrainVector[5] = 1.6941e-21;

    rF = ZeroMatrix(3, 3);
    rF(0,0) = 1.0;
    rF(1,1) = 1.0;
    rF(2,2) = (1.0 - 8.0e-5);
    rF(2,0) = 0.5 * -1.6941e-21;
    rF(0,2) = 0.5 * -1.6941e-21;
    double detF = MathUtils<double>::DetMat(rF);

    rTangentModuli = ZeroMatrix(6, 6);

    rCLConfigurationValues.SetMaterialProperties(rProperties);
    rCLConfigurationValues.SetDeformationGradientF(rF);
    rCLConfigurationValues.SetDeterminantF(detF);
    rCLConfigurationValues.SetStrainVector(rStrainVector);
    rCLConfigurationValues.SetStressVector(rStressVector);
    rCLConfigurationValues.SetOptions(constitutive_law_options);
    rCLConfigurationValues.SetProcessInfo(current_process_info);
    rCLConfigurationValues.SetConstitutiveMatrix(rTangentModuli);
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the elastic linear tensor
 */
KRATOS_TEST_CASE_IN_SUITE(LinearElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties material_properties;
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    SettingBasicCase(cl_configuration_values, material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("LinearElastic3DLaw").Clone();

    Matrix C = ZeroMatrix(6, 6);
    C = p_constitutive_law->CalculateValue(cl_configuration_values,CONSTITUTIVE_MATRIX, C);

    TangentOperatorCalculatorUtility::CalculateTangentTensor(cl_configuration_values, p_constitutive_law.get());

    const double tolerance = 1.0e-4;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            if (std::abs(C(i, j)) > 0.0) {
                KRATOS_CHECK_LESS_EQUAL((tangent_moduli(i, j) - C(i, j))/C(i, j), tolerance);
            } else if (std::abs(tangent_moduli(i, j)) > 0.0) {
                KRATOS_WARNING("LinearElasticCasePertubationTensorUtility") << "Be careful tangent_moduli(" << i << " ," << j << ") is greater tha 0: " <<  tangent_moduli(i, j) << std::endl;
            }
        }
    }
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the hyperelastic tensor
 */
KRATOS_TEST_CASE_IN_SUITE(HyperElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties material_properties;
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    SettingBasicCase(cl_configuration_values, material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, false);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("HyperElastic3DLaw").Clone();

    Matrix C = ZeroMatrix(6, 6);
    C = p_constitutive_law->CalculateValue(cl_configuration_values,CONSTITUTIVE_MATRIX_PK2, C);

    TangentOperatorCalculatorUtility::CalculateTangentTensorFiniteDeformation(cl_configuration_values, p_constitutive_law.get());

    const double tolerance = 1.0e-4;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            if (std::abs(C(i, j)) > 0.0) {
                KRATOS_CHECK_LESS_EQUAL((tangent_moduli(i, j) - C(i, j))/C(i, j), tolerance);
            } else if (std::abs(tangent_moduli(i, j)) > 0.0) {
                KRATOS_WARNING("HyperElasticCasePertubationTensorUtility") << "Be careful tangent_moduli(" << i << " ," << j << ") is greater tha 0: " <<  tangent_moduli(i, j) << std::endl;
            }
        }
    }
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the elastic linear tensor
 */
KRATOS_TEST_CASE_IN_SUITE(QuadraticLinearElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite2)
{
    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties material_properties;
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    SettingBasicCase(cl_configuration_values, material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("LinearElastic3DLaw").Clone();

    // The delta parameters
    const Vector initial_strain_vector = strain_vector;
    const Vector delta_strain_vector = strain_vector;
    double alpha = 1.0;

    // First error computation
    Vector expected_stress = stress_vector;
    Vector expected_previous_stress = stress_vector;
    Vector expected_delta_stress = ZeroVector(6);
    Vector computed_delta_stress = ZeroVector(6);

    const double quadratic_threshold = 1.8;

    const std::size_t max_number_iters = 10;
    std::vector<double> vector_errors(max_number_iters);
    for (std::size_t iter = 0; iter < max_number_iters; ++iter) {
        alpha *= 0.5;
        const Vector final_increment_strain_vector = alpha * delta_strain_vector;
        noalias(strain_vector) =  initial_strain_vector + final_increment_strain_vector;

        noalias(expected_previous_stress) = expected_stress;
        expected_stress = p_constitutive_law->CalculateValue(cl_configuration_values,CAUCHY_STRESS_VECTOR, expected_stress);
        noalias(expected_delta_stress) = (expected_stress - expected_previous_stress);

        stress_vector = p_constitutive_law->CalculateValue(cl_configuration_values,CAUCHY_STRESS_VECTOR, stress_vector);
        TangentOperatorCalculatorUtility::CalculateTangentTensor(cl_configuration_values, p_constitutive_law.get());

        noalias(computed_delta_stress) = prod(tangent_moduli, final_increment_strain_vector);
        const Vector aux_error_vector = computed_delta_stress - expected_delta_stress;
        const double error = norm_2(aux_error_vector);

        vector_errors[iter] = error;

//         KRATOS_CHECK_GREATER_EQUAL(error, quadratic_threshold);
    }

    for (int i = 0; i < max_number_iters - 4; ++i) {
        const double slope = std::log((vector_errors[i + 3] - vector_errors[i + 2])/(vector_errors[i + 2] - vector_errors[i + 1]))/std::log((vector_errors[i + 2] - vector_errors[i + 1])/(vector_errors[i + 1] - vector_errors[i + 0]));
        KRATOS_WATCH(slope)
    }
}

} // namespace Testing
} // namespace Kratos
