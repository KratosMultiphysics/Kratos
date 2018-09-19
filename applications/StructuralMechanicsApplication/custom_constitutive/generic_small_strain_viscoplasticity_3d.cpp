// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/viscous_generalized_maxwell.h"
#include "generic_small_strain_viscoplasticity_3d.h"

namespace Kratos
{
ConstitutiveLaw::Pointer GenericSmallStrainViscoplasticity3D::Create(Kratos::Parameters NewParameters) const
{
    ConstitutiveLaw::Pointer p_plasticity_cl = SmallStrainIsotropicPlasticityFactory().Create(NewParameters);
    ConstitutiveLaw::Pointer p_viscous_cl = Kratos::make_shared<ViscousGeneralizedMaxwell<ElasticIsotropic3D>>();
    return Kratos::make_shared<GenericSmallStrainViscoplasticity3D>(p_plasticity_cl, p_viscous_cl);

}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    ConstitutiveLaw::Pointer plaw     = this->GetPlasticityConstitutiveLaw();
    ConstitutiveLaw::Pointer viscolaw = this->GetViscousConstitutiveLaw();

    Vector plastic_strain = ZeroVector(6);
    plaw->GetValue(PLASTIC_STRAIN_VECTOR, plastic_strain);

    Vector& strain_vector = rValues.GetStrainVector();
    const Vector strain_for_visco = strain_vector - plastic_strain;
    const Vector initial_strain_vector = strain_vector;

    strain_vector = strain_for_visco;
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    viscolaw->CalculateMaterialResponseCauchy(rValues); // Relaxes the Stress...

    strain_vector = initial_strain_vector;
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Auxiliar! to be discussed in the comittee -> Do not be afraid VICENTE
    rValues.GetOptions().Set(ConstitutiveLaw::U_P_LAW, true);
    // ********************************************************************


    plaw->CalculateMaterialResponseCauchy(rValues); // Plastification occurs...

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::FinalizeSolutionStep(
    const Properties &rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Update the int vars of each SubConstitutiveLaw
    mpPlasticityConstitutiveLaw->FinalizeSolutionStep(rMaterialProperties, rElementGeometry,
                                                      rShapeFunctionsValues, rCurrentProcessInfo);

    mpViscousConstitutiveLaw->FinalizeSolutionStep(rMaterialProperties, rElementGeometry,
                                                   rShapeFunctionsValues, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::CalculateElasticMatrix(
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

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericSmallStrainViscoplasticity3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

double &GenericSmallStrainViscoplasticity3D::GetValue(
    const Variable<double> &rThisVariable,
    double &rValue)
{
    if (rThisVariable == UNIAXIAL_STRESS)
    {
        rValue = mpPlasticityConstitutiveLaw->GetValue(UNIAXIAL_STRESS, rValue);
    }
    else if (rThisVariable == PLASTIC_DISSIPATION)
    {
        rValue = mpPlasticityConstitutiveLaw->GetValue(PLASTIC_DISSIPATION, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector &GenericSmallStrainViscoplasticity3D::GetValue(
    const Variable<Vector> &rThisVariable,
    Vector &rValue)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR)
    {
        rValue = mpPlasticityConstitutiveLaw->GetValue(PLASTIC_STRAIN_VECTOR, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericSmallStrainViscoplasticity3D::Has(const Variable<double> &rThisVariable)
{
    if (rThisVariable == UNIAXIAL_STRESS)
    {
        return true;
    }
    else if (rThisVariable == PLASTIC_DISSIPATION)
    {
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

double &GenericSmallStrainViscoplasticity3D::CalculateValue(
    Parameters &rParameterValues,
    const Variable<double> &rThisVariable,
    double &rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

} // namespace Kratos
