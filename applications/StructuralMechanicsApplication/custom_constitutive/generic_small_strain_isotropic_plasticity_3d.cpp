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

// Project includes
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_constitutive/generic_small_strain_isotropic_plasticity_3d.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    // Integrate Stress plasticity
    const Properties &rMaterialProperties = rValues.GetMaterialProperties();
    const int VoigtSize = this->GetStrainSize();
    Vector &IntegratedStressVector = rValues.GetStressVector();
    Matrix &TangentTensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const double CharacteristicLength = rValues.GetElementGeometry().Length();
    const Flags &ConstitutiveLawOptions = rValues.GetOptions();

    // Elastic Matrix
    Matrix C;
    this->CalculateElasticMatrix(C, rMaterialProperties);

    double Threshold, PlasticDissipation;
    Vector PlasticStrain = ZeroVector(this->GetStrainSize());

    Threshold = this->GetThreshold();
    PlasticDissipation = this->GetPlasticDissipation();
    PlasticStrain = this->GetPlasticStrain();

    // S0 = C:(E-Ep)
    Vector PredictiveStressVector = prod(C, rValues.GetStrainVector() - PlasticStrain);

    // Initialize Plastic Parameters
    double UniaxialStress = 0.0, PlasticDenominator = 0.0;
    Vector Fflux = ZeroVector(VoigtSize), Gflux = ZeroVector(VoigtSize); // DF/DS & DG/DS
    Vector PlasticStrainIncrement = ZeroVector(VoigtSize);

    ConstLawIntegratorType::CalculatePlasticParameters(PredictiveStressVector, rValues.GetStrainVector(),
                                                       UniaxialStress, Threshold, PlasticDenominator, Fflux, Gflux, PlasticDissipation,
                                                       PlasticStrainIncrement, C, rMaterialProperties, CharacteristicLength);

    const double F = UniaxialStress - Threshold;

    if (F <= std::abs(1.0e-4 * Threshold))
    { // Elastic case
        noalias(IntegratedStressVector) = PredictiveStressVector;
        this->SetNonConvPlasticDissipation(PlasticDissipation);
        this->SetNonConvPlasticStrain(PlasticStrain);
        this->SetNonConvThreshold(Threshold);

        if (ConstitutiveLawOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR) == true)
        {
            noalias(TangentTensor) = C;
            this->SetValue(UNIAXIAL_STRESS, UniaxialStress, rValues.GetProcessInfo());
        }
    }
    else // Plastic case
    {
        // while loop backward euler
        /* Inside "IntegrateStressVector" the PredictiveStressVector
            is updated to verify the yield criterion */
        ConstLawIntegratorType::IntegrateStressVector(PredictiveStressVector, rValues.GetStrainVector(),
                                                      UniaxialStress, Threshold, PlasticDenominator, Fflux, Gflux,
                                                      PlasticDissipation, PlasticStrainIncrement, C, PlasticStrain,
                                                      rMaterialProperties, CharacteristicLength);
        noalias(IntegratedStressVector) = PredictiveStressVector;

        this->SetNonConvPlasticDissipation(PlasticDissipation);
        this->SetNonConvPlasticStrain(PlasticStrain);
        this->SetNonConvThreshold(Threshold);

        if (ConstitutiveLawOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR) == true)
        {
            this->CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
            noalias(TangentTensor) = rValues.GetConstitutiveMatrix();
            this->SetValue(UNIAXIAL_STRESS, UniaxialStress, rValues.GetProcessInfo());
        }
    }
    // this->SetValue(INTEGRATED_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(IntegratedStressVector), rValues.GetProcessInfo());
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues)
{
    // Calculates the Tangent Constitutive Tensor by perturbation
    TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties &rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    this->SetPlasticDissipation(this->GetNonConvPlasticDissipation());
    this->SetThreshold(this->GetNonConvThreshold());
    this->SetPlasticStrain(this->GetNonConvPlasticStrain());
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponse(
    ConstitutiveLaw::Parameters &rValues,
    const StressMeasure &rStressMeasure)
{
} // FinalizeMaterialResponse

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateElasticMatrix(
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

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
bool GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::Has(const Variable<double> &rThisVariable)
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

template <class ConstLawIntegratorType>
bool GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::Has(const Variable<Vector> &rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR)
    {
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::SetValue(
    const Variable<double> &rThisVariable,
    const double &rValue,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rThisVariable == UNIAXIAL_STRESS)
    {
        mUniaxialStress = rValue;
    }
    else if (rThisVariable == PLASTIC_DISSIPATION)
    {
        mPlasticDissipation = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::SetValue(
    const Variable<Vector> &rThisVariable,
    const Vector &rValue,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR)
    {
        mPlasticStrain = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
double &GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::GetValue(
    const Variable<double> &rThisVariable,
    double &rValue)
{
    if (rThisVariable == UNIAXIAL_STRESS)
    {
        rValue = mUniaxialStress;
    }
    else if (rThisVariable == PLASTIC_DISSIPATION)
    {
        rValue = mPlasticDissipation;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
Vector &GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::GetValue(
    const Variable<Vector> &rThisVariable,
    Vector &rValue)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR)
    {
        rValue = mPlasticStrain;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
Matrix &GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::GetValue(
    const Variable<Matrix> &rThisVariable,
    Matrix &rValue)
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR)
    {
        rValue = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
double &GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<double> &rThisVariable,
    double &rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
Matrix &GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<Matrix> &rThisVariable,
    Matrix &rValue)
{
    if (rThisVariable == INTEGRATED_STRESS_TENSOR)
    {
        //1.-Compute total deformation gradient
        const Matrix &DeformationGradientF = rParameterValues.GetDeformationGradientF();
        //2.-Right Cauchy-Green tensor C
        Matrix right_cauchy_green = prod(trans(DeformationGradientF), DeformationGradientF);
        Vector strain_vector = ZeroVector(6);

        //E= 0.5*(FT*F-1) or E = 0.5*(C-1)
        strain_vector[0] = 0.5 * (right_cauchy_green(0, 0) - 1.00);
        strain_vector[1] = 0.5 * (right_cauchy_green(1, 1) - 1.00);
        strain_vector[2] = 0.5 * (right_cauchy_green(2, 2) - 1.00);
        strain_vector[3] = right_cauchy_green(0, 1); // xy
        strain_vector[4] = right_cauchy_green(1, 2); // yz
        strain_vector[5] = right_cauchy_green(0, 2); // xz

        Matrix C;
        const Properties &MaterialProperties = rParameterValues.GetMaterialProperties();
        this->CalculateElasticMatrix(C, MaterialProperties);

        rValue = MathUtils<double>::StressVectorToTensor(prod(C, strain_vector - mPlasticStrain));
        return rValue;
    }
    else
    {
        return this->GetValue(rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential>>>;

} // namespace Kratos
