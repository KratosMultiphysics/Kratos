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

/**
 * @brief This sets the basic case
 */
void SettingBasicCase(
    ConstitutiveLaw::Parameters& rCLConfigurationValues,
    Properties& rProperties,
    Vector& rStressVector,
    Vector& rStrainVector,
    Matrix& rTangentModuli,
    Matrix& rF,
    double& rDetF,
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

    rStrainVector = ZeroVector(6);
    rStrainVector[2] = 8.0e-5;
    rStrainVector[5] = 1.6941e-21;

    // Compute equivalent F
    rF = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i) {
        rF(i, i) = 1.0 + rStrainVector[i];
    }

    for (int i = 3; i < 6; ++i) {
        const int equivalent_i = (i == 3) ? 0 : (i == 4) ? 1 : 0;
        const int equivalent_j = (i == 3) ? 1 : 2;
        rF(equivalent_i, equivalent_j) = 0.5 * rStrainVector[i];
        rF(equivalent_j, equivalent_i) = 0.5 * rStrainVector[i];
    }
    rDetF = MathUtils<double>::DetMat(rF);

    rTangentModuli = ZeroMatrix(6, 6);

    rCLConfigurationValues.SetMaterialProperties(rProperties);
    rCLConfigurationValues.SetDeformationGradientF(rF);
    rCLConfigurationValues.SetDeterminantF(rDetF);
    rCLConfigurationValues.SetStrainVector(rStrainVector);
    rCLConfigurationValues.SetStressVector(rStressVector);
    rCLConfigurationValues.SetOptions(constitutive_law_options);
    rCLConfigurationValues.SetProcessInfo(current_process_info);
    rCLConfigurationValues.SetConstitutiveMatrix(rTangentModuli);
}

/**
 * This computes the convergence rate of the CL
 */
void ComputingConvergenceRate(
    ConstitutiveLaw::Pointer pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rCLConfigurationValues,
    Vector& rStressVector,
    Vector& rStrainVector,
    Matrix& rTangentModuli,
    Matrix& rF,
    double& rDetF,
    const bool FiniteDeformation = false,
    const bool Debug = false
    )
{
    // The delta parameters
    const Vector initial_strain_vector = rStrainVector;
    const Vector initial_stress_vector = rStressVector;
    const Matrix initial_deformation_gradient_F = rF;
//     const double initial_det_deformation_gradient_F = rDetF;
    Matrix delta_deformation_gradient_F = rF;
    double alpha = 1.0;

    // Ensure the proper flag
    Flags& cl_options = rCLConfigurationValues.GetOptions();
    const bool use_element_provided_strain = cl_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // First error computation
    Vector expected_delta_stress = ZeroVector(6);
    Vector computed_delta_stress = ZeroVector(6);

    Vector delta_strain_vector = rStrainVector;

    const double quadratic_threshold = 1.8;

    const std::size_t max_number_iters = 10;
    std::vector<double> vector_errors(max_number_iters);
    for (std::size_t iter = 0; iter < max_number_iters; ++iter) {
        alpha *= 0.5;

        noalias(delta_strain_vector) = alpha * initial_strain_vector;
        noalias(rStrainVector) =  initial_strain_vector + delta_strain_vector;
        if (!use_element_provided_strain) {
            for (int i = 0; i < 3; ++i) {
                delta_deformation_gradient_F(i, i) = 1.0 + delta_strain_vector[i];
            }

            for (int i = 3; i < 6; ++i) {
                const int equivalent_i = (i == 3) ? 0 : (i == 4) ? 1 : 0;
                const int equivalent_j = (i == 3) ? 1 : 2;
                delta_deformation_gradient_F(equivalent_i, equivalent_j) = 0.5 * delta_strain_vector[i];
                delta_deformation_gradient_F(equivalent_j, equivalent_i) = 0.5 * delta_strain_vector[i];
            }
            noalias(rF) = prod(delta_deformation_gradient_F, initial_deformation_gradient_F);
            rDetF = MathUtils<double>::DetMat(rF);
            rCLConfigurationValues.SetDeformationGradientF(rF);
            rCLConfigurationValues.SetDeterminantF(rDetF);
        }

        pConstitutiveLaw->CalculateMaterialResponse(rCLConfigurationValues, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
        noalias(expected_delta_stress) = (rStressVector - initial_stress_vector);

        if (FiniteDeformation) {
            TangentOperatorCalculatorUtility::CalculateTangentTensorFiniteDeformation(rCLConfigurationValues, pConstitutiveLaw.get(), ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
        } else {
            TangentOperatorCalculatorUtility::CalculateTangentTensor(rCLConfigurationValues, pConstitutiveLaw.get(), ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
        }

//         const Matrix C = rTangentModuli;
//
//         cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
//         pConstitutiveLaw->CalculateMaterialResponse(rCLConfigurationValues, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
//         cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
//
//         const Matrix aux_error = rTangentModuli - C;
//         const double error = norm_frobenius(aux_error);

        noalias(computed_delta_stress) = prod(rTangentModuli, delta_strain_vector);
        const Vector aux_error_vector = computed_delta_stress - expected_delta_stress;
        const double error = norm_2(aux_error_vector);

        vector_errors[iter] = error;

        if (Debug) {
            KRATOS_WATCH(error)
        }
    }

    for (int i = 1; i < max_number_iters - 4; ++i) { // We discard the first solution
        if ((vector_errors[i + 3] + vector_errors[i + 2] + vector_errors[i + 1] + vector_errors[i]) > 1.0e-6) { // If zero means exact solution
            const double slope = std::log((vector_errors[i + 3] - vector_errors[i + 2])/(vector_errors[i + 2] - vector_errors[i + 1]))/std::log((vector_errors[i + 2] - vector_errors[i + 1])/(vector_errors[i + 1] - vector_errors[i + 0]));

            if (Debug) {
                KRATOS_WATCH(slope)
            } else { // Check
                KRATOS_CHECK_GREATER_EQUAL(slope, quadratic_threshold);
            }
        }
    }
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
    double det_deformation_gradient_F;
    SettingBasicCase(cl_configuration_values, material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("LinearElastic3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_Cauchy);

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
    double det_deformation_gradient_F;
    SettingBasicCase(cl_configuration_values, material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("HyperElastic3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

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
KRATOS_TEST_CASE_IN_SUITE(QuadraticLinearElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties material_properties;
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(cl_configuration_values, material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("LinearElastic3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

    ComputingConvergenceRate(p_constitutive_law, cl_configuration_values, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false, false);
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the ehyper lastic tensor
 */
KRATOS_TEST_CASE_IN_SUITE(QuadraticHyperElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite2)
{
    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties material_properties;
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(cl_configuration_values, material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("HyperElastic3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

    ComputingConvergenceRate(p_constitutive_law, cl_configuration_values, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, true, true);
}

} // namespace Testing
} // namespace Kratos
