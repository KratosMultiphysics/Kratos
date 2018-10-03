// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"

// Application includes

// Contitutive Law
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

namespace Kratos
{
namespace Testing
{

void CalculateElasticMatrix(
    Matrix &rElasticityTensor,
    const Properties &rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    const double mu = E / (2.0 + 2.0 * poisson_ratio);

    if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
        rElasticityTensor.resize(6, 6, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2.0 * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2.0 * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2.0 * mu;
    rElasticityTensor(3, 3) = mu;
    rElasticityTensor(4, 4) = mu;
    rElasticityTensor(5, 5) = mu;
}

/**
    * Check the correct calculation of the uniaxial stress of the yield surfaces
    */
KRATOS_TEST_CASE_IN_SUITE(PertubationTensorTestUtility, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters rValues;
    Properties rMaterialProperties;
    Vector rStressVector, rStrainVector;
    ProcessInfo CurrentProcessInfo;

    Flags ConstitutiveLawOptions;
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    rMaterialProperties.SetValue(YOUNG_MODULUS, 210e9);
    rMaterialProperties.SetValue(POISSON_RATIO, 0.22);

    rStressVector = ZeroVector(6);
    rStressVector[0] = 5.40984e+06;
    rStressVector[1] = 5.40984e+06;
    rStressVector[2] = 1.91803e+07;
    rStressVector[3] = 0.0;
    rStressVector[4] = 0.0;
    rStressVector[5] = 1.45804e-10;

    rStrainVector = ZeroVector(6);
    rStrainVector[0] = 0.0;
    rStrainVector[1] = 0.0;
    rStrainVector[2] = 8.0e-5;
    rStrainVector[3] = 0.0;
    rStrainVector[4] = 0.0;
    rStrainVector[5] = 1.6941e-21;

    Matrix F = IdentityMatrix(6);

    rValues.SetMaterialProperties(rMaterialProperties);
    rValues.SetDeformationGradientF(F);
    rValues.SetStrainVector(rStrainVector);
    rValues.SetStressVector(rStressVector);
    rValues.SetOptions(ConstitutiveLawOptions);
    rValues.SetProcessInfo(CurrentProcessInfo);

    ElasticIsotropic3D ConstitutiveLaw = ElasticIsotropic3D();
    ElasticIsotropic3D *pConstitutiveLaw = &ConstitutiveLaw;

    Matrix C = ZeroMatrix(6, 6);
    rValues.SetConstitutiveMatrix(C);
    CalculateElasticMatrix(C, rMaterialProperties);

    TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, pConstitutiveLaw);
    Matrix &Tangent = rValues.GetConstitutiveMatrix();

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            KRATOS_CHECK_NEAR(C(i, j), Tangent(i, j), 1.0e-3);
        }
    }
}
} // namespace Testing
} // namespace Kratos