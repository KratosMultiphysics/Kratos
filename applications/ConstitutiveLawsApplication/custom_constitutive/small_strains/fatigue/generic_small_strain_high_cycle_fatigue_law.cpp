// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
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
#include "constitutive_laws_application_variables.h"
#include "constitutive_laws_application_variables.h"
#include "generic_small_strain_high_cycle_fatigue_law.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/high_cycle_fatigue_law_integrator.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"

namespace Kratos
{

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    if (rElementGeometry.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rElementGeometry.GetValue(REFERENCE_TEMPERATURE);
    } else if (rMaterialProperties.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rMaterialProperties[REFERENCE_TEMPERATURE];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::InitializeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::InitializeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::InitializeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const bool current_load_type = rValues.GetProcessInfo()[CURRENT_LOAD_TYPE];
    const double max_stress = mMaxStress;
    const double min_stress = mMinStress;
    bool max_indicator = mMaxDetected;
    bool min_indicator = mMinDetected;
    double fatigue_reduction_factor = mFatigueReductionFactor;
    double reversion_factor_relative_error = mReversionFactorRelativeError;
    double max_stress_relative_error = mMaxStressRelativeError;
    unsigned int global_number_of_cycles = mNumberOfCyclesGlobal;
    unsigned int local_number_of_cycles = mNumberOfCyclesLocal;
    double B0 = mFatigueReductionParameter;
    double previous_max_stress = mPreviousMaxStress;
    double previous_min_stress = mPreviousMinStress;
    double wohler_stress = mWohlerStress;
    bool new_cycle = false;
    double s_th = mThresholdStress;
    double cycles_to_failure = mCyclesToFailure;
    bool advance_in_time_process_applied = rValues.GetProcessInfo()[ADVANCE_STRATEGY_APPLIED];
    // bool no_linearity_activation = rValues.GetProcessInfo()[NO_LINEARITY_ACTIVATION];
    double c_factor = mCFactor;

    const bool new_model_part = rValues.GetProcessInfo()[NEW_MODEL_PART];
    if (new_model_part) {
        // KRATOS_WATCH("HERE")
        max_indicator = false;
        min_indicator = false;
        mFirstCycleOfANewLoad = true;
    }
    const double damage = this->GetDamage();
    const double reference_damage = mReferenceDamage;    //Threshold is used here to define Nf. This is required for those cases that damage has started
    double threshold = this->GetThreshold();

    if (max_indicator && min_indicator && current_load_type) {
        if (mFirstCycleOfANewLoad) {
            const SizeType fatigue_parameters_size = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS].size();
            if (fatigue_parameters_size == 8) {
                c_factor = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][7];
            } else if (fatigue_parameters_size == 11) {
                c_factor = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][7] * ((1.0 - reference_damage) * max_stress) + rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][8];
                const double c_factor_min = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][9];
                const double c_factor_max = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][10];

                KRATOS_ERROR_IF(c_factor_min > c_factor_max) << "The min and max C factor order is not correct: first C_min and then C_max" << std::endl;

                c_factor = (c_factor < c_factor_min) ? c_factor_min : c_factor;
                c_factor = (c_factor > c_factor_max) ? c_factor_max : c_factor;
            }
        }

        const double previous_reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(previous_max_stress, previous_min_stress);
        const double reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(max_stress, min_stress); // Does not depend on the reference damage
        double alphat;
        const double ultimate_stress = HighCycleFatigueLawIntegrator<6>::UltimateStressDamage(rValues, mReferenceTemperature);
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(
            (1.0 - reference_damage) * max_stress,
            reversion_factor,
            rValues,
            B0,
            s_th,
            alphat,
            cycles_to_failure,
            ultimate_stress,
            c_factor,
            mReferenceTemperature);
        cycles_to_failure = HighCycleFatigueLawIntegrator<6>::NumberOfCyclesToFailure(
            cycles_to_failure,
            (1.0 - reference_damage) * max_stress, //Damaged maximum stress considering the reference state
            rValues.GetMaterialProperties(),
            (1.0 - damage) * threshold, //Current damage threshold with no influence of fatigue, only damage no-linearity
            s_th,
            ultimate_stress,
            c_factor);

        double betaf = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];
        if (std::abs(min_stress) < 0.001) {
            reversion_factor_relative_error = std::abs(reversion_factor - previous_reversion_factor);
        } else {
            reversion_factor_relative_error = std::abs((reversion_factor - previous_reversion_factor) / reversion_factor);
        }
        max_stress_relative_error = std::abs((max_stress - previous_max_stress) / max_stress);

        // if (!damage_activation && global_number_of_cycles > 2 && !advance_in_time_process_applied && (reversion_factor_relative_error > 0.001 || max_stress_relative_error > 0.001)) {
        if (mFirstCycleOfANewLoad && global_number_of_cycles > 2 && !advance_in_time_process_applied && (reversion_factor_relative_error > 0.001 || max_stress_relative_error > 0.1) && ((1.0 - reference_damage) * max_stress >= s_th)) {
            local_number_of_cycles = std::trunc(std::pow(10, std::pow(-(std::log(fatigue_reduction_factor) / B0), 1.0 / (betaf * betaf * c_factor)))) + 1;
        }
        global_number_of_cycles++;
        local_number_of_cycles++;
        new_cycle = true;
        max_indicator = false;
        min_indicator = false;
        previous_max_stress = max_stress;
        previous_min_stress = min_stress;
        mCyclesToFailure = cycles_to_failure;

        HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(rValues.GetMaterialProperties(),
                                                                                        (1.0 - reference_damage) * max_stress,
                                                                                        local_number_of_cycles,
                                                                                        global_number_of_cycles,
                                                                                        B0,
                                                                                        s_th,
                                                                                        alphat,
                                                                                        fatigue_reduction_factor,
                                                                                        wohler_stress,
                                                                                        ultimate_stress,
                                                                                        c_factor);
        mFirstCycleOfANewLoad = false;
    }
    if (advance_in_time_process_applied && current_load_type) {
        // KRATOS_WATCH("ADV IN TIME")
        const double reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(max_stress, min_stress);
        double alphat;
        const double ultimate_stress = HighCycleFatigueLawIntegrator<6>::UltimateStressDamage(rValues, mReferenceTemperature);
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(
            (1.0 - reference_damage) * max_stress,
            reversion_factor,
            rValues,
            B0,
            s_th,
            alphat,
            cycles_to_failure,
            ultimate_stress,
            c_factor,
            mReferenceTemperature);
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(rValues.GetMaterialProperties(),
                                                                                        (1.0 - reference_damage) * max_stress,
                                                                                        local_number_of_cycles,
                                                                                        global_number_of_cycles,
                                                                                        B0,
                                                                                        s_th,
                                                                                        alphat,
                                                                                        fatigue_reduction_factor,
                                                                                        wohler_stress,
                                                                                        ultimate_stress,
                                                                                        c_factor);
    }
    mNumberOfCyclesGlobal = global_number_of_cycles;
    mNumberOfCyclesLocal = local_number_of_cycles;
    mReversionFactorRelativeError = reversion_factor_relative_error;
    mMaxStressRelativeError = max_stress_relative_error;
    mMaxDetected = max_indicator;
    mMinDetected = min_indicator;
    mFatigueReductionParameter = B0;
    mPreviousMaxStress = previous_max_stress;
    mPreviousMinStress = previous_min_stress;
    mFatigueReductionFactor = fatigue_reduction_factor;
    mWohlerStress = wohler_stress;
    mNewCycleIndicator = new_cycle;
    mThresholdStress = s_th;
    mCFactor = c_factor;
    if (new_model_part) {
        mReferenceDamage = this->GetDamage();   //Updating the damage reference values. This needs to be changed by the end of the method because the calculations
                                                //done here are built using the values of the previous step. This should not have a big effect in this CL but is consistent.
    }
}


/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
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
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // Elastic Matrix
    // if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
    //     this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    // }

    // We compute the stress
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        // this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        const double E = AdvCLutils::GetMaterialPropertyThroughAccessor(YOUNG_MODULUS, rValues);
        const double poisson_ratio = AdvCLutils::GetMaterialPropertyThroughAccessor(POISSON_RATIO, rValues);

        if constexpr (Dimension == 2) {
            CLutils::CalculateElasticMatrixPlaneStrain(r_constitutive_matrix, E, poisson_ratio);
            AdvCLutils::SubstractThermalStrain(r_strain_vector, mReferenceTemperature, rValues, true);
        } else if constexpr (Dimension == 3) {
            CLutils::CalculateElasticMatrix(r_constitutive_matrix, E, poisson_ratio);
            AdvCLutils::SubstractThermalStrain(r_strain_vector, mReferenceTemperature, rValues);
        }

        // Converged values
        double threshold = this->GetThreshold();
        double damage = this->GetDamage();

        // S0 = C:E
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);

        // Initialize Plastic Parameters
        double fatigue_reduction_factor = mFatigueReductionFactor;
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        const double ref_yield = AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS, rValues, mReferenceTemperature);
        const double current_yield = AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS, rValues);
        const double temperature_reduction_factor = current_yield / ref_yield;

        uniaxial_stress /= fatigue_reduction_factor;      // Fatigue contribution
        uniaxial_stress /= temperature_reduction_factor;  // Thermal contribution
        const double F = uniaxial_stress - threshold;

        if (F <= threshold_tolerance) { // Elastic case
            noalias(r_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                r_constitutive_matrix *= (1.0 - damage);
                this->SetStressVector(r_integrated_stress_vector);
                rValues.SetStressVector(r_integrated_stress_vector);
            }
        } else { // Damage case
            const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
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
                rValues.SetStressVector(r_integrated_stress_vector);
                this->SetStressVector(rValues.GetStressVector());
                this->CalculateTangentTensor(rValues);
            }
        }
    }
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
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // Elastic Matrix
    // if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
    //     this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    // }
    array_1d<double, VoigtSize> predictive_stress_vector;
    // We compute the stress
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {

        // Elastic Matrix
        // this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        const double E = AdvCLutils::GetMaterialPropertyThroughAccessor(YOUNG_MODULUS, rValues);
        const double poisson_ratio = AdvCLutils::GetMaterialPropertyThroughAccessor(POISSON_RATIO, rValues);

        if constexpr (Dimension == 2) {
            CLutils::CalculateElasticMatrixPlaneStrain(r_constitutive_matrix, E, poisson_ratio);
            AdvCLutils::SubstractThermalStrain(r_strain_vector, mReferenceTemperature, rValues, true);
        } else if constexpr (Dimension == 3) {
            CLutils::CalculateElasticMatrix(r_constitutive_matrix, E, poisson_ratio);
            AdvCLutils::SubstractThermalStrain(r_strain_vector, mReferenceTemperature, rValues);
        }

        // Converged values
        double threshold = this->GetThreshold();
        double damage = this->GetDamage();

        // S0 = C:E
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);

        // Initialize Plastic Parameters
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(predictive_stress_vector);
        uniaxial_stress *= sign_factor;
        double max_stress = mMaxStress;
        double min_stress = mMinStress;
        bool max_indicator = mMaxDetected;
        bool min_indicator = mMinDetected;
        double fatigue_reduction_factor = mFatigueReductionFactor;

        HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(
            uniaxial_stress,
            max_stress,
            min_stress,
            mPreviousStresses,
            max_indicator,
            min_indicator);
        mMaxStress = max_stress;
        mMinStress = min_stress;
        mMaxDetected = max_indicator;
        mMinDetected = min_indicator;

        uniaxial_stress *= sign_factor;

        const double ref_yield = AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS, rValues, mReferenceTemperature);
        const double current_yield = AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS, rValues);
        const double temperature_reduction_factor = current_yield / ref_yield;

        uniaxial_stress /= fatigue_reduction_factor;      // Fatigue contribution
        uniaxial_stress /= temperature_reduction_factor;  // Thermal contribution

        const double F = uniaxial_stress - threshold;

        if (F > threshold_tolerance) {
                const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
                // This routine updates the PredictiveStress to verify the yield surface
                TConstLawIntegratorType::IntegrateStressVector(
                    predictive_stress_vector,
                    uniaxial_stress,
                    damage,
                    threshold,
                    rValues,
                    characteristic_length);
                this->SetDamage(damage);
                this->SetThreshold(uniaxial_stress);
        } else {
            predictive_stress_vector *= (1.0 - this->GetDamage());
            TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);
        }

        Vector previous_stresses = ZeroVector(2);
        const Vector& r_aux_stresses = mPreviousStresses;
        previous_stresses[1] = this->CalculateValue(rValues, UNIAXIAL_STRESS, previous_stresses[1]) * sign_factor / (1.0 - this->GetDamage());
        previous_stresses[0] = r_aux_stresses[1];
        mPreviousStresses = previous_stresses;
    }
}
/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::Has(const Variable<bool>& rThisVariable)
{
    if (rThisVariable == CYCLE_INDICATOR) {
        return true;
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
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES) {
        return true;
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{

    if (rThisVariable == FATIGUE_REDUCTION_FACTOR) {
        return true;
    } else if (rThisVariable == WOHLER_STRESS) {
        return true;
    } else if (rThisVariable == CYCLES_TO_FAILURE) {
        return true;
    } else if (rThisVariable == REVERSION_FACTOR_RELATIVE_ERROR) {
        return true;
    } else if (rThisVariable == MAX_STRESS_RELATIVE_ERROR) {
        return true;
    } else if (rThisVariable == MAX_STRESS) {
        return true;
    } else if (rThisVariable == THRESHOLD_STRESS) {
        return true;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        return true;
    } else if (rThisVariable == CYCLE_PERIOD) {
        return true;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == CYCLE_INDICATOR) {
        mNewCycleIndicator = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == NUMBER_OF_CYCLES) {
        mNumberOfCyclesGlobal = rValue;
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES) {
        mNumberOfCyclesLocal = rValue;
    }
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

    if (rThisVariable == FATIGUE_REDUCTION_FACTOR) {
        mFatigueReductionFactor = rValue;
    } else if (rThisVariable == WOHLER_STRESS) {
        mWohlerStress = rValue;
    } else if (rThisVariable == CYCLES_TO_FAILURE) {
        mCyclesToFailure = rValue;
    } else if (rThisVariable == REVERSION_FACTOR_RELATIVE_ERROR) {
        mReversionFactorRelativeError = rValue;
    } else if (rThisVariable == MAX_STRESS_RELATIVE_ERROR) {
        mMaxStressRelativeError = rValue;
    } else if (rThisVariable == MAX_STRESS) {
        mMaxStress = rValue;
    } else if (rThisVariable == THRESHOLD_STRESS) {
        mThresholdStress = rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        mPreviousCycleTime = rValue;
    } else if (rThisVariable == CYCLE_PERIOD) {
        mPeriod = rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        mPreviousCycleDamage = rValue;
    } else {
        return BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool& GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    if (rThisVariable == CYCLE_INDICATOR) {
        rValue = mNewCycleIndicator;
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
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES) {
        rValue = mNumberOfCyclesLocal;
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{

    if (rThisVariable == FATIGUE_REDUCTION_FACTOR) {
        rValue = mFatigueReductionFactor;
    } else if (rThisVariable == WOHLER_STRESS) {
        rValue = mMinStress/mMaxStress;
    } else if (rThisVariable == CYCLES_TO_FAILURE) {
        rValue = mCyclesToFailure;
    } else if (rThisVariable == REVERSION_FACTOR_RELATIVE_ERROR) {
        rValue = mReversionFactorRelativeError;
    } else if (rThisVariable == MAX_STRESS_RELATIVE_ERROR) {
        rValue = mMaxStressRelativeError;
    } else if (rThisVariable == MAX_STRESS) {
        rValue = mMaxStress;
    } else if (rThisVariable == THRESHOLD_STRESS) {
        rValue = mThresholdStress;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        rValue = mPreviousCycleTime;
    } else if (rThisVariable == CYCLE_PERIOD) {
        rValue = mPeriod;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        rValue = mPreviousCycleDamage;
    } else if (rThisVariable == INFINITY_YIELD_STRESS) {
        rValue = mMinStress / mMaxStress;
    } else if (rThisVariable == DAMAGE_MATRIX) {
        rValue = mCFactor;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (this->Has(rThisVariable))
        return this->GetValue(rThisVariable, rValue);
    else
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
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

template class GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;

} // namespace Kratos