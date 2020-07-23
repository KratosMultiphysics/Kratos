// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jimenez/Alejandro Cornejo/Lucia Barbu
//  Collaborator:
//

// System includes

// External includes

// Project includes
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "structural_mechanics_application_variables.h"
#include "custom_advanced_constitutive/generic_small_strain_high_cycle_fatigue_law.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/high_cycle_fatigue_law_integrator.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const Flags& r_constitutive_law_options = rValues.GetOptions();
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }

    // We compute the stress
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        // Converged values
        double threshold = this->GetThreshold();
        double damage = this->GetDamage();

        // S0 = C:E
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);

        // Initialize Plastic Parameters
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        double max_stress = this->GetMaxStress();
        double min_stress = this->GetMinStress();
        double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(predictive_stress_vector);
        uniaxial_stress *= sign_factor;
        unsigned int global_number_of_cycles = this->GetNumberOfCyclesGlobal();
        unsigned int local_number_of_cycles = this->GetNumberOfCyclesLocal();
        bool max_indicator = this->GetMaxDetected();
        bool min_indicator = this->GetMinDetected();
        double fatigue_reduction_factor = this->GetFatigueReductionFactor();
        double B0 = this->GetFatigueReductionParameter();
        double previous_max_stress = this->GetPreviousMaxStress();
        double previous_min_stress = this->GetPreviousMinStress();
        double wohler_stress = this->GetWohlerStress();


        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(
                uniaxial_stress,
                max_stress,
                min_stress,
                this->GetPreviousStresses(),
                max_indicator,
                min_indicator);
            this->SetMaxStress(max_stress);
            this->SetMinStress(min_stress);

            if (max_indicator && min_indicator) {
                double previous_reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(previous_max_stress, previous_min_stress);
                double reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(max_stress, min_stress);

                double s_th, alphat;
                HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(
                    max_stress,
                    reversion_factor,
                    rValues.GetMaterialProperties(),
                    B0,
                    s_th,
                    alphat);

                double betaf = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];
                double reversion_factor_relative_error = std::abs((previous_reversion_factor - reversion_factor) / previous_reversion_factor);
                double max_stress_relative_error = std::abs((previous_reversion_factor - reversion_factor) / previous_reversion_factor);

                if (global_number_of_cycles > 2 && (reversion_factor_relative_error > 0.001 || max_stress_relative_error > 0.001)) {
                    local_number_of_cycles = std::trunc(std::pow(10, std::pow(-(std::log(fatigue_reduction_factor) / B0), 1.0 / (betaf * betaf)))) + 1;
                }
                global_number_of_cycles++;
                local_number_of_cycles++;
                max_indicator = false;
                min_indicator = false;
                previous_max_stress = max_stress;
                previous_min_stress = min_stress;

                HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(rValues.GetMaterialProperties(),
                                                                                                max_stress,
                                                                                                local_number_of_cycles,
                                                                                                global_number_of_cycles,
                                                                                                B0,
                                                                                                s_th,
                                                                                                alphat,
                                                                                                fatigue_reduction_factor,
                                                                                                wohler_stress);
            }
            this->SetNumberOfCyclesGlobal(global_number_of_cycles);
            this->SetNumberOfCyclesLocal(local_number_of_cycles);
            this->SetMaxDetected(max_indicator);
            this->SetMinDetected(min_indicator);
            this->SetFatigueReductionParameter(B0);
            this->SetFatigueReductionFactor(fatigue_reduction_factor);
            this->SetPreviousMaxStress(previous_max_stress);
            this->SetPreviousMinStress(previous_min_stress);
            this->SetWohlerStress(wohler_stress);
        }

        uniaxial_stress *= sign_factor;
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());
        }

        uniaxial_stress /= fatigue_reduction_factor;  // Fatigue contribution
        const double F = uniaxial_stress - threshold;

        if (F <= tolerance) { // Elastic case
            noalias(r_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                r_constitutive_matrix *= (1.0 - damage);
                this->SetStressVector(r_integrated_stress_vector);
            }
        } else { // Damage case
            const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
            // This routine updates the PredictiveStress to verify the yield surf
            TConstLawIntegratorType::IntegrateStressVector(
                predictive_stress_vector,
                uniaxial_stress,
                damage,
                threshold,
                rValues,
                characteristic_length);

            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector,uniaxial_stress, rValues);
                this->SetStressVector(r_integrated_stress_vector);
                this->CalculateTangentTensor(rValues);
            }
        }
    }
}


/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    const Flags& r_constitutive_law_options = rValues.GetOptions();
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }
    array_1d<double, VoigtSize> predictive_stress_vector;
    // We compute the stress
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        // Converged values
        double threshold = this->GetThreshold();
        double damage = this->GetDamage();

        // S0 = C:E
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);

        // Initialize Plastic Parameters
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);
        this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());

        double fatigue_reduction_factor = this->GetFatigueReductionFactor();

        uniaxial_stress /= fatigue_reduction_factor;  // Fatigue contribution
        const double F = uniaxial_stress - threshold;

        if (F > tolerance) {
            const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
            // This routine updates the PredictiveStress to verify the yield surf
            TConstLawIntegratorType::IntegrateStressVector(
                predictive_stress_vector,
                uniaxial_stress,
                damage,
                threshold,
                rValues,
                characteristic_length);
            this->SetDamage(damage);
            this->SetThreshold(uniaxial_stress);
        }
    }
    double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(predictive_stress_vector);
    Vector previous_stresses = ZeroVector(2);
    const Vector& r_aux_stresses = this->GetPreviousStresses();
    previous_stresses[1] = this->GetValue(UNIAXIAL_STRESS, previous_stresses[1])*sign_factor;
    previous_stresses[0] = r_aux_stresses[1];
    this->SetPreviousStresses(previous_stresses);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == DAMAGE || rThisVariable == UNIAXIAL_STRESS || rThisVariable == THRESHOLD) {
        GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::Has(rThisVariable);
    } else if (rThisVariable == FATIGUE_REDUCTION_FACTOR) {
        return true;
    } else if (rThisVariable == WOHLER_STRESS) {
        return true;
    } else {
        return GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::Has(rThisVariable);
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::Has(const Variable<int>& rThisVariable)
{
    if (rThisVariable == NUMBER_OF_CYCLES) {
        return true;
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == DAMAGE || rThisVariable == UNIAXIAL_STRESS || rThisVariable == THRESHOLD) {
        GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    } else if (rThisVariable == FATIGUE_REDUCTION_FACTOR) {
        mFatigueReductionFactor = rValue;
    } else if (rThisVariable == WOHLER_STRESS) {
        mWohlerStress = rValue;
    } else {
        return GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue
    )
{
    if (rThisVariable == NUMBER_OF_CYCLES) {
        mNumberOfCyclesGlobal = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == DAMAGE || rThisVariable == UNIAXIAL_STRESS || rThisVariable == THRESHOLD) {
        GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::GetValue(rThisVariable, rValue);
    } else if (rThisVariable == FATIGUE_REDUCTION_FACTOR) {
        rValue = mFatigueReductionFactor;
    } else if (rThisVariable == WOHLER_STRESS) {
        rValue = mWohlerStress;
    } else {
        return GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int& GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    if (rThisVariable == NUMBER_OF_CYCLES) {
        rValue = mNumberOfCyclesGlobal;
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        rValue = MathUtils<double>::StressVectorToTensor(this->GetStressVector());
    } else if (rThisVariable == CONSTITUTIVE_MATRIX) {
        this->CalculateElasticMatrix(rValue, rParameterValues);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>>;

} // namespace Kratos
