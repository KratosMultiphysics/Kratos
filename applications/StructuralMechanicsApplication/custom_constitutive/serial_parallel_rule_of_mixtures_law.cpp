// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo 
//                   Vicente Mataix
//                   Fernando Rastellini
//  Collaborator:    Lucia Barbu
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
ConstitutiveLaw::Pointer SerialParallelRuleOfMixturesLaw::Create(Kratos::Parameters NewParameters) const
{


}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
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
    viscolaw->CalculateMaterialResponseCauchy(rValues); // Viscous Process

    strain_vector = initial_strain_vector;
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Flag to tell the plasticity to take the predictor from the viscous law
    rValues.GetOptions().Set(ConstitutiveLaw::U_P_LAW, true);

    plaw->CalculateMaterialResponseCauchy(rValues); // Plastic Process

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeSolutionStep(
    const Properties &rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Deprecated
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateElasticMatrix(
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

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    mpViscousConstitutiveLaw->FinalizeMaterialResponseCauchy(rValues);
    mpPlasticityConstitutiveLaw->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

double &SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<double> &rThisVariable,
    double &rValue)
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        rValue = mpPlasticityConstitutiveLaw->GetValue(UNIAXIAL_STRESS, rValue);
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mpPlasticityConstitutiveLaw->GetValue(PLASTIC_DISSIPATION, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector &SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<Vector> &rThisVariable,
    Vector &rValue)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mpPlasticityConstitutiveLaw->GetValue(PLASTIC_STRAIN_VECTOR, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<double> &rThisVariable)
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        return true;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

double &SerialParallelRuleOfMixturesLaw::CalculateValue(
    Parameters &rParameterValues,
    const Variable<double> &rThisVariable,
    double &rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

} // namespace Kratos
