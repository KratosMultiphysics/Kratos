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

/* Contitutive Law */
#include "custom_utilities/tangent_operator_calculator_utility.h"

namespace Kratos
{
namespace Testing
{

/**
* Check the correct calculation of the uniaxial stress of the yield surfaces
*/
KRATOS_TEST_CASE_IN_SUITE(PertubationTensorTestUtility, KratosStructuralMechanicsFastSuite)
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

    Matrix F = IdentityMatrix(6);

    cl_configuration_values.SetMaterialProperties(material_properties);
    cl_configuration_values.SetDeformationGradientF(F);
    cl_configuration_values.SetStrainVector(strain_vector);
    cl_configuration_values.SetStressVector(stress_vector);
    cl_configuration_values.SetOptions(constitutive_law_options);
    cl_configuration_values.SetProcessInfo(CurrentProcessInfo);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("LinearElastic3DLaw").Clone();

    Matrix C = ZeroMatrix(6, 6);
    cl_configuration_values.SetConstitutiveMatrix(C);
    p_constitutive_law->CalculateValue(cl_configuration_values,CONSTITUTIVE_MATRIX, C);

    TangentOperatorCalculatorUtility::CalculateTangentTensor(cl_configuration_values, p_constitutive_law.get());
    Matrix& r_tangent_moduli = cl_configuration_values.GetConstitutiveMatrix();

    const double tolerance = 1.0e-6;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            if (std::abs(r_tangent_moduli(i, j)) > 0.0) {
                KRATOS_CHECK_NEAR(C(i, j), r_tangent_moduli(i, j), tolerance);
            }
        }
    }
}
} // namespace Testing
} // namespace Kratos
