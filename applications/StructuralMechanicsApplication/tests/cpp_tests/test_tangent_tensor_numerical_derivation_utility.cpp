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
 * @brief This test tests that the perturbation utility is valid for computing the elastic linear tensor
 */
KRATOS_TEST_CASE_IN_SUITE(LinearElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties material_properties;
    Vector stress_vector, strain_vector;
    ProcessInfo CurrentProcessInfo;

    Flags constitutive_law_options;
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);

    stress_vector = ZeroVector(6);
    stress_vector[0] = 5.40984e+06;
    stress_vector[1] = 5.40984e+06;
    stress_vector[2] = 1.91803e+07;
    stress_vector[5] = 1.45804e-10;

    strain_vector = ZeroVector(6);
    strain_vector[2] = 8.0e-5;
    strain_vector[5] = 1.6941e-21;

    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient(0,0) = 1.0;
    deformation_gradient(1,1) = 1.0;
    deformation_gradient(2,2) = (1.0 - 8.0e-5);
    deformation_gradient(2,0) = 0.5 * -1.6941e-21;
    deformation_gradient(0,2) = 0.5 * -1.6941e-21;
    double detF = MathUtils<double>::DetMat(deformation_gradient);

    Matrix tangent_moduli = ZeroMatrix(6, 6);

    cl_configuration_values.SetMaterialProperties(material_properties);
    cl_configuration_values.SetDeformationGradientF(deformation_gradient);
    cl_configuration_values.SetDeterminantF(detF);
    cl_configuration_values.SetStrainVector(strain_vector);
    cl_configuration_values.SetStressVector(stress_vector);
    cl_configuration_values.SetOptions(constitutive_law_options);
    cl_configuration_values.SetProcessInfo(CurrentProcessInfo);
    cl_configuration_values.SetConstitutiveMatrix(tangent_moduli);

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
    ProcessInfo CurrentProcessInfo;

    Flags constitutive_law_options;
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);

    stress_vector = ZeroVector(6);
    stress_vector[0] = 5.40984e+06;
    stress_vector[1] = 5.40984e+06;
    stress_vector[2] = 1.91803e+07;
    stress_vector[5] = 1.45804e-10;

    strain_vector = ZeroVector(6);
    strain_vector[2] = 8.0e-5;
    strain_vector[5] = 1.6941e-21;

    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient(0,0) = 1.0;
    deformation_gradient(1,1) = 1.0;
    deformation_gradient(2,2) = (1.0 - 8.0e-5);
    deformation_gradient(2,0) = 0.5 * -1.6941e-21;
    deformation_gradient(0,2) = 0.5 * -1.6941e-21;
    const double detF = MathUtils<double>::DetMat(deformation_gradient);

    Matrix tangent_moduli = ZeroMatrix(6, 6);

    cl_configuration_values.SetMaterialProperties(material_properties);
    cl_configuration_values.SetDeformationGradientF(deformation_gradient);
    cl_configuration_values.SetDeterminantF(detF);
    cl_configuration_values.SetStrainVector(strain_vector);
    cl_configuration_values.SetStressVector(stress_vector);
    cl_configuration_values.SetOptions(constitutive_law_options);
    cl_configuration_values.SetProcessInfo(CurrentProcessInfo);
    cl_configuration_values.SetConstitutiveMatrix(tangent_moduli);

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

} // namespace Testing
} // namespace Kratos
