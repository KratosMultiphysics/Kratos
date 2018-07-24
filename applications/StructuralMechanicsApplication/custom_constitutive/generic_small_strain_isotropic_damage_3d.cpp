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
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/generic_small_strain_isotropic_damage_3d.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

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
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    // Integrate Stress Damage
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Vector& integrated_stress_vector = rValues.GetStressVector();
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // Elastic Matrix
    Matrix C;
    this->CalculateElasticMatrix(C, r_material_properties);

    double threshold, damage;
    // In the 1st step Threshold must be set
    if (std::abs(this->GetThreshold()) < tolerance) {
        ConstLawIntegratorType::YieldSurfaceType::GetInitialUniaxialThreshold(r_material_properties, threshold);
        this->SetThreshold(threshold);
    }

    // Converged values
    threshold = this->GetThreshold();
    damage = this->GetDamage();

    // S0 = C:(E-Ep)
    Vector predictive_stress_vector = prod(C, rValues.GetStrainVector());

    // Initialize Plastic Parameters
    double uniaxial_stress;
    ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector,
                                                                        rValues.GetStrainVector(), uniaxial_stress, r_material_properties);

    const double F = uniaxial_stress - threshold;

    if (F <= 0.0) { // Elastic case

        this->SetNonConvDamage(damage);
        this->SetNonConvThreshold(threshold);
        noalias(integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;

        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR) == true) {
            noalias(r_tangent_tensor) = (1.0 - damage) * C;

            ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(integrated_stress_vector,
                                                                                rValues.GetStrainVector(), uniaxial_stress, r_material_properties);

            this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());
        }
    } else { // Damage case
        const double characteristic_length = rValues.GetElementGeometry().Length();
        // This routine updates the PredictiveStress to verify the yield surf
        ConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, uniaxial_stress,
                                                      damage, threshold, r_material_properties,
                                                      characteristic_length);

        // Updated Values
        noalias(integrated_stress_vector) = predictive_stress_vector;
        this->SetNonConvDamage(damage);
        this->SetNonConvThreshold(uniaxial_stress);

        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR) == true) {
            this->CalculateTangentTensor(rValues);
            noalias(r_tangent_tensor) = rValues.GetConstitutiveMatrix();

            ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(integrated_stress_vector,
                                         rValues.GetStrainVector(), uniaxial_stress, r_material_properties);

            this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());
        }
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues)
{
    // Calculates the Tangent Constitutive Tensor by perturbation
    TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    this->SetDamage(this->GetNonConvDamage());
    this->SetThreshold(this->GetNonConvThreshold());
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateElasticMatrix(
    Matrix& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda = E * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
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
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
bool GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == DAMAGE) {
        return true;
    } else if (rThisVariable == THRESHOLD) {
        return true;
    } else if (rThisVariable == UNIAXIAL_STRESS) {
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == DAMAGE) {
        mDamage = rValue;
    } else if (rThisVariable == THRESHOLD) {
        mThreshold = rValue;
    } else if (rThisVariable == UNIAXIAL_STRESS) {
        mUniaxialStress = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
double& GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if (rThisVariable == DAMAGE) {
        rValue = mDamage;
    } else if (rThisVariable == THRESHOLD) {
        rValue = mThreshold;
    } else if (rThisVariable == UNIAXIAL_STRESS) {
        rValue = mUniaxialStress;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
double& GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateValue(
    Parameters &rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
Matrix& GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<Matrix> &rThisVariable,
    Matrix& rValue)
{
    if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        //1.-Compute total deformation gradient
        const Matrix& deformation_gradient_F = rParameterValues.GetDeformationGradientF();
        //2.-Right Cauchy-Green tensor C
        Matrix right_cauchy_green = prod(trans(deformation_gradient_F), deformation_gradient_F);
        Vector strain_vector = ZeroVector(6);

        //E= 0.5*(FT*F-1) or E = 0.5*(C-1)
        strain_vector[0] = 0.5 * (right_cauchy_green(0, 0) - 1.00);
        strain_vector[1] = 0.5 * (right_cauchy_green(1, 1) - 1.00);
        strain_vector[2] = 0.5 * (right_cauchy_green(2, 2) - 1.00);
        strain_vector[3] = right_cauchy_green(0, 1); // xy
        strain_vector[4] = right_cauchy_green(1, 2); // yz
        strain_vector[5] = right_cauchy_green(0, 2); // xz

        Matrix C;
        const Properties& MaterialProperties = rParameterValues.GetMaterialProperties();
        this->CalculateElasticMatrix(C, MaterialProperties);

        Vector stress = prod(C, strain_vector);
        stress *= (1.0 - mDamage);
        rValue =  MathUtils<double>::StressVectorToTensor(stress);
        return rValue;
    }
	return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicDamage3D<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential>>>;

} // namespace Kratos
