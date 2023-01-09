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
//  Main authors:    Sergio Jimenez
//                   Alejandro Cornejo
//  Collaborator:    Lucia Barbu
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "unified_fatigue_rule_of_mixtures_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_constitutive/constitutive_laws_integrators/high_cycle_fatigue_law_integrator.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity_with_fatigue.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template<class TConstLawIntegratorType>
ConstitutiveLaw::Pointer UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::Create(Kratos::Parameters NewParameters) const
{
    const double isotropic_damage_initial_volumetric_participation = NewParameters["combination_factors"][0].GetDouble();
    return Kratos::make_shared<UnifiedFatigueRuleOfMixturesLaw>(isotropic_damage_initial_volumetric_participation);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::InitializeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::InitializeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::InitializeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // //FATIGUE CASES
    // const double max_stress = mMaxStress;
    // const double min_stress = mMinStress;
    // bool max_indicator = mMaxDetected;
    // bool min_indicator = mMinDetected;
    // double fatigue_reduction_factor = mFatigueReductionFactor;
    // double reversion_factor_relative_error = mReversionFactorRelativeError;
    // double max_stress_relative_error = mMaxStressRelativeError;
    // unsigned int global_number_of_cycles = mNumberOfCyclesGlobal;
    // unsigned int local_number_of_cycles = mNumberOfCyclesLocal;
    // double B0 = mFatigueReductionParameter;
    // double previous_max_stress = mPreviousMaxStress;
    // double previous_min_stress = mPreviousMinStress;
    // double wohler_stress = mWohlerStress;
    // bool new_cycle = false;
    // double s_th = mFatigueLimit;
    // double cycles_to_failure = mCyclesToFailure;
    // bool advance_in_time_process_applied = rValues.GetProcessInfo()[ADVANCE_STRATEGY_APPLIED];
    // double c_factor = mCFactor;

    // // bool no_linearity_activation = rValues.GetProcessInfo()[NO_LINEARITY_ACTIVATION];
    // const bool current_load_type = rValues.GetProcessInfo()[CURRENT_LOAD_TYPE];
    // const bool new_model_part = rValues.GetProcessInfo()[NEW_MODEL_PART];
    // const double damage = mDamage;
    // const double plastic_dissipation = mPlasticDissipation;
    // const double new_load_bock_damage = mNewLoadBlockDamage;
    // const double new_load_bock_plastic_dissipation = mNewLoadBlockPlasticDissipation;
    // const double new_load_block_fatigue_reduction_factor = mNewLoadBlockFatigueReductionFactor;

    // if (new_model_part) {
    //     max_indicator = false; //New model part indicates restarting the cycle counter.
    //     min_indicator = false;
    //     cycles_to_failure = 0.0; //Each time that a new load block is detected the most restrictive case (ULCF) is considered.
    //                                                 //The behaviour will be updated depending on the Nf value after a cycle.
    //     // mFirstCycleOfANewLoad = true;
    // }

    // auto& r_material_properties = rValues.GetMaterialProperties();
    // const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    // const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
    // const auto& r_props_plasticity_cl = *(it_cl_begin + 1);
    // ConstitutiveLaw::Parameters values_fatigue  = rValues;
    // ConstitutiveLaw::Parameters values_damage_component  = rValues;
    // ConstitutiveLaw::Parameters values_plasticity_component = rValues;
    // values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
    // values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);

    // //Checking which material has the fatigue properties
    // if (r_props_isotropic_damage_cl.Has(HIGH_CYCLE_FATIGUE_COEFFICIENTS)) {
    //     values_fatigue.SetMaterialProperties(r_props_isotropic_damage_cl);
    // } else if (r_props_plasticity_cl.Has(HIGH_CYCLE_FATIGUE_COEFFICIENTS)) {
    //     values_fatigue.SetMaterialProperties(r_props_plasticity_cl);
    // } else {
    //     KRATOS_ERROR << "Fatigue properties not defined" << std::endl;
    // }

    // if (max_indicator && min_indicator && current_load_type) {
    //     if (mFirstCycleOfANewLoad) {
    //         const SizeType fatigue_parameters_size = values_fatigue.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS].size();
    //         if (fatigue_parameters_size == 8) {
    //             c_factor = values_fatigue.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][7];
    //         } else if (fatigue_parameters_size == 9) {
    //             c_factor = values_fatigue.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][7] * max_stress + values_fatigue.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][8];
    //         }
    //     }
    //     //Checking if there is any no-linearity accumulation
    //     mAcumulatedDamageCurrentCycle = (damage - mPreviousCycleDamage > machine_tolerance) ? true : false;
    //     mAcumulatedPlasticityCurrentCycle = (plastic_dissipation - mPreviousCyclePlasticDissipation > machine_tolerance) ? true : false;

    //     const double previous_reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(previous_max_stress, previous_min_stress);
    //     const double reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(max_stress, min_stress);
    //     double alphat;

    //     double ultimate_stress = CalculateUltimateStress(values_damage_component.GetMaterialProperties(), values_plasticity_component.GetMaterialProperties());
    //     //Check to ensure that we are working in a feasible stress regime
    //     // KRATOS_ERROR_IF(mMaxStressDamageBranch > ultimate_stress || mMaxStressPlasticityBranch > ultimate_stress) << "Component cycle stresses outside the composite stress regime" << std::endl;

    //     HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(
    //         max_stress,
    //         reversion_factor,
    //         values_fatigue.GetMaterialProperties(),
    //         B0,
    //         s_th,
    //         alphat,
    //         cycles_to_failure,
    //         ultimate_stress,
    //         c_factor);

    //     double betaf = values_fatigue.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];
    //     if (std::abs(min_stress) < 0.001) {
    //         reversion_factor_relative_error = std::abs(reversion_factor - previous_reversion_factor);
    //     } else {
    //         reversion_factor_relative_error = std::abs((reversion_factor - previous_reversion_factor) / reversion_factor);
    //     }
    //     if (std::abs(max_stress) < 0.001) {
    //         max_stress_relative_error = std::abs(max_stress - previous_max_stress);
    //     } else {
    //         max_stress_relative_error = std::abs((max_stress - previous_max_stress) / max_stress);
    //     }
    //     if (mFirstCycleOfANewLoad && global_number_of_cycles > 2 && !advance_in_time_process_applied && (reversion_factor_relative_error > 0.001 || max_stress_relative_error > 0.1) && (max_stress >= s_th)) {
    //         local_number_of_cycles = std::trunc(std::pow(10, std::pow(-(std::log(fatigue_reduction_factor) / B0), 1.0 / (betaf * betaf * c_factor)))) + 1;
    //     }

    //     global_number_of_cycles++;

    //     //The local_number_of_cycles are only updated when fred degradation is expected, i.e. HCF and LCF. ULCF showing no plasticity neither damage accumulation are also counting Nlocal cycles
    //     local_number_of_cycles = ((!mAcumulatedDamageCurrentCycle && !mAcumulatedPlasticityCurrentCycle) || cycles_to_failure > 1.0e3) ? local_number_of_cycles + 1.0 : local_number_of_cycles;

    //     //Idfentifying the number of cycles that iniciate the no-linearities (either damage or plasticity)
    //     const double cycles_to_failure_damage = HighCycleFatigueLawIntegrator<6>::NumberOfCyclesToFailure(
    //         cycles_to_failure,
    //         mMaxStressDamageBranch,
    //         values_fatigue.GetMaterialProperties(),
    //         (1.0 - damage) * mIsotropicDamageThreshold, //Current damage threshold with no influence of fatigue, only damage no-linearity
    //         s_th,
    //         ultimate_stress, //Composite ultimate stress which is the one used to compute the original number of cycles to failure
    //         c_factor);
    //     const double cycles_to_failure_plasticity = HighCycleFatigueLawIntegrator<6>::NumberOfCyclesToFailure(
    //         cycles_to_failure,
    //         mMaxStressPlasticityBranch,
    //         values_fatigue.GetMaterialProperties(),
    //         mPlasticityThreshold / fatigue_reduction_factor, //Plasticity threshold with no influence of fatigue, only damage no-linearity
    //         s_th,
    //         ultimate_stress, //Composite ultimate stress which is the one used to compute the original number of cycles to failure
    //         c_factor);
    //     cycles_to_failure = std::min(cycles_to_failure_damage, cycles_to_failure_plasticity);

    //     new_cycle = true;
    //     max_indicator = false;
    //     min_indicator = false;
    //     previous_max_stress = max_stress;
    //     previous_min_stress = min_stress;

    //     HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(values_fatigue.GetMaterialProperties(),
    //                                                                                     max_stress,
    //                                                                                     local_number_of_cycles,
    //                                                                                     global_number_of_cycles,
    //                                                                                     B0,
    //                                                                                     s_th,
    //                                                                                     alphat,
    //                                                                                     fatigue_reduction_factor,
    //                                                                                     wohler_stress,
    //                                                                                     ultimate_stress,
    //                                                                                     c_factor);
    //     mFirstCycleOfANewLoad = false;
    // }
    // mCyclesToFailure = cycles_to_failure;

    // if (advance_in_time_process_applied && current_load_type) {
    //     const double reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(max_stress, min_stress);
    //     double alphat;
    //     double ultimate_stress = CalculateUltimateStress(values_damage_component.GetMaterialProperties(), values_plasticity_component.GetMaterialProperties());
    //     HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(
    //         max_stress,
    //         reversion_factor,
    //         values_fatigue.GetMaterialProperties(),
    //         B0,
    //         s_th,
    //         alphat,
    //         cycles_to_failure,
    //         ultimate_stress,
    //         c_factor);
    //     HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(values_fatigue.GetMaterialProperties(),
    //                                                                                     max_stress,
    //                                                                                     local_number_of_cycles,
    //                                                                                     global_number_of_cycles,
    //                                                                                     B0,
    //                                                                                     s_th,
    //                                                                                     alphat,
    //                                                                                     fatigue_reduction_factor,
    //                                                                                     wohler_stress,
    //                                                                                     ultimate_stress,
    //                                                                                     c_factor);
    // }
    // mNumberOfCyclesGlobal = global_number_of_cycles;
    // mNumberOfCyclesLocal = local_number_of_cycles;
    // mReversionFactorRelativeError = reversion_factor_relative_error;
    // mMaxStressRelativeError = max_stress_relative_error;
    // mMaxDetected = max_indicator;
    // mMinDetected = min_indicator;
    // mFatigueReductionParameter = B0;
    // mPreviousMaxStress = previous_max_stress;
    // mPreviousMinStress = previous_min_stress;
    // mFatigueReductionFactor = fatigue_reduction_factor;
    // mWohlerStress = wohler_stress;
    // mNewCycleIndicator = new_cycle;
    // mFatigueLimit = s_th;
    // mCFactor = c_factor;

    // //Adapting the volumetric participation to the type of load that is being applied
    // mIsotrpicDamageVolumetricParticipation = VolumetricParticipationUpdate(
    //     current_load_type, cycles_to_failure,
    //     damage, plastic_dissipation, fatigue_reduction_factor,
    //     new_load_bock_damage, new_load_bock_plastic_dissipation, new_load_block_fatigue_reduction_factor,
    //     values_plasticity_component);


    // if (new_model_part) {   //Updating the reference values. This needs to be changed by the end of the method because the calculation of the volumetric participation here is done
    //                         //with results of the non-linear indicators from previous step and so the volumetric participation computed here corresponds to the one for previous step
    //     mReferenceDamage = damage;
    //     mReferencePlasticStrain = mPlasticStrain;

    //     mNewLoadBlockDamage = damage;
    //     mNewLoadBlockPlasticDissipation = plastic_dissipation;
    //     mNewLoadBlockFatigueReductionFactor = fatigue_reduction_factor;
    //     mNewLoadBlockVolumetricParticipation = mIsotrpicDamageVolumetricParticipation; // is not necessary to prevent here the d=0 && kp=0 case.
    // }


    //LOADING-UNLOADING CASE

}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateUltimateStress(
    const Properties&  rValuesDamageComponentProperties,
    const Properties&  rValuesPlasticityComponentProperties
    )
{
    const double ultimate_stress_damage_component = HighCycleFatigueLawIntegrator<6>::UltimateStressDamage(rValuesDamageComponentProperties);
    const double ultimate_stress_plasticity_component = GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::UltimateStressPlasticity(rValuesDamageComponentProperties);
    return mIsotrpicDamageVolumetricParticipation * ultimate_stress_damage_component + (1.0 - mIsotrpicDamageVolumetricParticipation) * ultimate_stress_plasticity_component;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::VolumetricParticipationUpdate(
    const bool CurrentLoadType,
    const double GlobalCyclesToFailure,
    const double CurrentDamage,
    const double CurrentPlasticDissipation,
    const double CurrentFatigueReductionFactor,
    const double NewLoadBlockDamage,
    const double NewLoadBlockPlasticDissipation,
    const double NewLoadBlockFatigueReductionFactor,
    ConstitutiveLaw::Parameters& rValuesPlasticityComponent
    )
{
    if (!CurrentLoadType) {       //Monotonic load being applied, i.e., no more damage is accumulated and plasticity is isotropic.
        const double new_load_block_volumetric_participation = (CurrentDamage == 0.0 && CurrentPlasticDissipation == 0.0) ? 0.0 : mNewLoadBlockVolumetricParticipation; //Fully plasticity behaviour while non-linear process has not started.
        return (CurrentDamage > 0.0) ? new_load_block_volumetric_participation * NewLoadBlockDamage / CurrentDamage : 0.0;

    } else if (CurrentLoadType) { //Cyclic load being applied. three cases depending on Nf.
        // if (GlobalCyclesToFailure < 1.0e3) {               // 1 - ULCF case = no more damage is accumulated and plasticity is kinematic.
        //     //Isotropic plasticity is used meanwhile; using kinematic plasticity is only interesting when
        //     const double new_load_block_volumetric_participation = (CurrentDamage == 0.0 && CurrentPlasticDissipation == 0.0) ? 0.0 : mNewLoadBlockVolumetricParticipation; //Fully plasticity behaviour while non-linear process has not started.
        //     return (CurrentDamage > 0.0) ? new_load_block_volumetric_participation * NewLoadBlockDamage / CurrentDamage : 0.0;

        // // } else if (Nf < 1.0e5) {        // 2 - LCF case = intermidiate case. Volumetric participation computed through linear transition.

        // } else {                        // 3 - HCF case = no more plasticity is accumulated.
            double new_load_block_equivalent_plastic_strain, equivalent_plastic_strain;
            const double new_load_block_volumetric_participation = (CurrentDamage == 0.0 && CurrentPlasticDissipation == 0.0) ? 1.0 : mNewLoadBlockVolumetricParticipation; //Fully damage behaviour while non-linear process has not started.
            const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValuesPlasticityComponent.GetElementGeometry());

            GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::EquivalencyPlasticDissipationUniaxialPlasticStrain(new_load_block_equivalent_plastic_strain,
                                                                                                                                                                rValuesPlasticityComponent, NewLoadBlockPlasticDissipation,
                                                                                                                                                                NewLoadBlockFatigueReductionFactor,
                                                                                                                                                                characteristic_length);
            GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::EquivalencyPlasticDissipationUniaxialPlasticStrain(equivalent_plastic_strain,
                                                                                                                                                                rValuesPlasticityComponent, CurrentPlasticDissipation,
                                                                                                                                                                CurrentFatigueReductionFactor,
                                                                                                                                                                characteristic_length);

            // This function has an indeterminancy when
            //      (div by 0.0 indeterminancy) NewLoadBlockDamage == 1.0 && new_load_block_volumetric_participation == 1.0. When this happens, CurrentDamage == 1.0. This situation is not possible because mDamage < 0.99999 as set in the damage integrator
            //      (0.0 / 0.0) equivalent_plastic_strain == 0.0 && new_load_block_equivalent_plastic_strain == 0.0. This occurs while plasticity is not activated and so the plastic dissipation is 0.
            return  (CurrentPlasticDissipation > 0.0) ? ((1.0 - new_load_block_volumetric_participation) * new_load_block_equivalent_plastic_strain - (1.0 - new_load_block_volumetric_participation * NewLoadBlockDamage) * equivalent_plastic_strain)
                    / (CurrentDamage * (1.0 - new_load_block_volumetric_participation) * new_load_block_equivalent_plastic_strain - (1.0 - new_load_block_volumetric_participation * NewLoadBlockDamage) * equivalent_plastic_strain) : 1.0;
        // }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliar values
    const SizeType dimension = WorkingSpaceDimension();
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

    const Properties& r_material_properties  = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }
    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        Matrix F_deformation_gradient(dimension, dimension);
        this->CalculateValue(rValues, DEFORMATION_GRADIENT, F_deformation_gradient);
        const Matrix B_matrix = prod(F_deformation_gradient, trans(F_deformation_gradient));
        // Doing resize in case is needed
        if (r_strain_vector.size() != voigt_size)
            r_strain_vector.resize(voigt_size);

         // Identity matrix
        Matrix identity_matrix(dimension, dimension);
        for (IndexType i = 0; i < dimension; ++i) {
            for (IndexType j = 0; j < dimension; ++j) {
                if (i == j) identity_matrix(i, j) = 1.0;
                else identity_matrix(i, j) = 0.0;
            }
        }

        // Calculating the inverse of the left Cauchy tensor
        Matrix inverse_B_tensor(dimension, dimension);
        double aux_det_b = 0;
        MathUtils<double>::InvertMatrix(B_matrix, inverse_B_tensor, aux_det_b);

        // Calculate E matrix
        const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
        // Almansi Strain Calculation
        r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false); //Tangent tensor is not computed at the component level (only at the composite level)
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector hich is equal to component CLs strains
        Vector& r_strain_vector = rValues.GetStrainVector();

        // This method integrates the stress according to each simple material CL
        Vector isotropic_damage_stress_vector, plasticity_stress_vector;
        this->IntegrateStressesOfIsotropicDamageAndPlasticityModels(rValues, r_strain_vector, r_strain_vector, isotropic_damage_stress_vector, plasticity_stress_vector);

        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        noalias(r_integrated_stress_vector) = mIsotrpicDamageVolumetricParticipation * isotropic_damage_stress_vector
                                     + (1.0 - mIsotrpicDamageVolumetricParticipation) * plasticity_stress_vector;
        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->CalculateTangentTensor(rValues);
        }
    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::IntegrateStressesOfIsotropicDamageAndPlasticityModels(
    ConstitutiveLaw::Parameters& rValues,
    Vector rIsotropicDamgeStrainVector,
    Vector rPlasticityStrainVector,
    Vector& rIsotropicDamageStressVector,
    Vector& rPlasticityStressVector
)
{
    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
    const auto& r_props_plasticity_cl = *(it_cl_begin + 1);

    ConstitutiveLaw::Parameters values_damage_component  = rValues;
    ConstitutiveLaw::Parameters values_plasticity_component = rValues;

    values_damage_component.SetStrainVector(rIsotropicDamgeStrainVector);
    values_plasticity_component.SetStrainVector(rPlasticityStrainVector);


    // Integrate Stress of the HCF part
    values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
    this->CalculateMaterialResponseIsotropicDamageModelModel(values_damage_component);
    // mpIsotropicDamageConstitutiveLaw->CalculateMaterialResponseCauchy(values_damage_component);
    rIsotropicDamageStressVector = values_damage_component.GetStressVector();

    // Integrate Stress of the UÑCF part
    values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);
    this->CalculateMaterialResponsePlasticityModel(values_plasticity_component);
    // mpPlasticityConstitutiveLaw->CalculateMaterialResponseCauchy(values_plasticity_component);
    rPlasticityStressVector = values_plasticity_component.GetStressVector();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateMaterialResponseIsotropicDamageModelModel(
        ConstitutiveLaw::Parameters& values_damage_component
)
{
    Vector& integrated_stress_vector = values_damage_component.GetStressVector();
    array_1d<double, VoigtSize> auxiliar_integrated_stress_vector = integrated_stress_vector;
    Matrix& r_tangent_tensor = values_damage_component.GetConstitutiveMatrix(); // todo modify after integration
    const Flags& r_constitutive_law_options = values_damage_component.GetOptions();
    Matrix& r_constitutive_matrix = values_damage_component.GetConstitutiveMatrix();
    mpIsotropicDamageConstitutiveLaw->CalculateValue(values_damage_component, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

    Vector& r_strain_vector = values_damage_component.GetStrainVector();
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // Converged values
    double isotropic_damage_threshold = this->GetIsotropicDamageThreshold();
    double damage = this->GetDamage();

    // S0 = C:(E-E0) + S0
    array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
    this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_stress_vector);

    // Initialize Plastic Parameters
    double fatigue_reduction_factor = mFatigueReductionFactor;
    double uniaxial_stress;
    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, values_damage_component);
    uniaxial_stress /= fatigue_reduction_factor;  // Fatigue contribution
    const double F = uniaxial_stress - isotropic_damage_threshold;

    if (F <= threshold_tolerance) { // Elastic case
        noalias(auxiliar_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;
		noalias(integrated_stress_vector) = auxiliar_integrated_stress_vector;

        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            noalias(r_tangent_tensor) = (1.0 - damage) * r_constitutive_matrix;
        }
    } else { // Damage case
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_damage_component.GetElementGeometry());
        // This routine updates the PredictiveStress to verify the yield surf
        TConstLawIntegratorType::IntegrateStressVector(
            predictive_stress_vector,
            uniaxial_stress, damage,
            isotropic_damage_threshold, values_damage_component,
            characteristic_length);

        // Updated Values
        noalias(auxiliar_integrated_stress_vector) = predictive_stress_vector;

        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->CalculateTangentTensor(values_damage_component);
        }
        noalias(integrated_stress_vector) = auxiliar_integrated_stress_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateMaterialResponsePlasticityModel(
        ConstitutiveLaw::Parameters& values_plasticity_component
)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = values_plasticity_component.GetOptions();
    // We get the strain vector
    Vector& r_strain_vector = values_plasticity_component.GetStrainVector();
    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = values_plasticity_component.GetConstitutiveMatrix();

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = values_plasticity_component.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ) ) {
            BaseType::CalculateCauchyGreenStrain( values_plasticity_component, r_strain_vector);
        }
        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = values_plasticity_component.GetStressVector();
            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateElasticMatrix( r_constitutive_matrix, values_plasticity_component);
                noalias(r_stress_vector) = prod( r_constitutive_matrix, r_strain_vector);
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            } else {
                BaseType::CalculatePK2Stress( r_strain_vector, r_stress_vector, values_plasticity_component);
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& r_integrated_stress_vector = values_plasticity_component.GetStressVector();
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_plasticity_component.GetElementGeometry());

        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( values_plasticity_component, r_strain_vector);
        }

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double plasticity_threshold = this->GetPlasticityThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain = this->GetPlasticStrain();

            array_1d<double, VoigtSize> predictive_stress_vector;
            if (r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW)) {
                noalias(predictive_stress_vector) = values_plasticity_component.GetStressVector();
            } else {
                // S0 = Elastic stress with strain (E-Ep) + S0
                Vector aux_stress = ZeroVector(VoigtSize);
                BaseType::CalculatePK2Stress(r_strain_vector - plastic_strain, aux_stress, values_plasticity_component);
                this->template AddInitialStressVectorContribution<Vector>(aux_stress);
                noalias(predictive_stress_vector) = aux_stress;
            }

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            array_1d<double, VoigtSize> f_flux = ZeroVector(VoigtSize); // DF/DS
            array_1d<double, VoigtSize> g_flux = ZeroVector(VoigtSize); // DG/DS
            array_1d<double, VoigtSize> plastic_strain_increment = ZeroVector(VoigtSize);

            // Elastic Matrix
            BaseType::CalculateElasticMatrix(r_constitutive_matrix, values_plasticity_component);

            double fatigue_reduction_factor = mFatigueReductionFactor;
            // Compute the plastic parameters
            const double F = GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                plasticity_threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, values_plasticity_component, characteristic_length,
                plastic_strain, fatigue_reduction_factor);

            const double yield_criteria = uniaxial_stress - plasticity_threshold;

            if (yield_criteria <= std::abs(1.0e-4 * plasticity_threshold)) { // Elastic case
                noalias(r_integrated_stress_vector) = predictive_stress_vector;
            } else { // Plastic case
                // While loop backward euler
                /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
                GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::IntegrateStressVector(
                    predictive_stress_vector, r_strain_vector, uniaxial_stress,
                    plasticity_threshold, plastic_denominator, f_flux, g_flux,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, plastic_strain, values_plasticity_component,
                    characteristic_length, fatigue_reduction_factor);
                noalias(r_integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(values_plasticity_component); // this modifies the ConstitutiveMatrix
                } else {
                    BaseType::CalculateElasticMatrix( r_constitutive_matrix, values_plasticity_component);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Deprecated
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& r_strain_vector = rValues.GetStrainVector();

    // Recalculation to obtain the serial_strain_matrix and store the value
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector& r_strain_vector = rValues.GetStrainVector();

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        // We call the FinalizeMaterialResponse of the HCF and ULCF CL
        auto& r_material_properties = rValues.GetMaterialProperties();
        const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
        const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
        const auto& r_props_plasticity_cl = *(it_cl_begin + 1);

        ConstitutiveLaw::Parameters values_damage_component  = rValues;
        ConstitutiveLaw::Parameters values_plasticity_component = rValues;

        values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
        values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);

        values_damage_component.SetStrainVector(r_strain_vector);
        values_plasticity_component.SetStrainVector(r_strain_vector);

        double isotropic_damage_predictive_uniaxial_stress, plasticity_predictive_uniaxial_stress;
        Vector isotropic_damage_stress_vector, plasticity_stress_vector;

        this->FinalizeMaterialResponseIsotropicDamageModel(isotropic_damage_predictive_uniaxial_stress, isotropic_damage_stress_vector, values_damage_component);
        this->FinalizeMaterialResponsePlastitictyModel(plasticity_predictive_uniaxial_stress, plasticity_stress_vector, values_plasticity_component);

        //Is this the correct way to compute it?? f(S) = kd * (C0 * E) + kp (C0 * (E-Ep)) tendría que ser directamente C0 * E?
        double predictive_uniaxial_stress;
        predictive_uniaxial_stress = mIsotrpicDamageVolumetricParticipation * isotropic_damage_predictive_uniaxial_stress
                                     + (1.0 - mIsotrpicDamageVolumetricParticipation) * plasticity_predictive_uniaxial_stress;

        Vector r_integrated_stress_vector;
        r_integrated_stress_vector = mIsotrpicDamageVolumetricParticipation * isotropic_damage_stress_vector
                                     + (1.0 - mIsotrpicDamageVolumetricParticipation) * plasticity_stress_vector;

        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(r_integrated_stress_vector, r_strain_vector, uniaxial_stress, rValues);
        double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(r_integrated_stress_vector);
        predictive_uniaxial_stress *= sign_factor;
        double max_stress = mMaxStress;
        double min_stress = mMinStress;
        bool max_indicator = mMaxDetected;
        bool min_indicator = mMinDetected;

        HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(
            predictive_uniaxial_stress,
            max_stress,
            min_stress,
            mPreviousStresses,
            max_indicator,
            min_indicator);

        mMaxStress = max_stress;
        mMinStress = min_stress;
        if (!mMaxDetected && max_indicator) {
            mMaxStressDamageBranch = isotropic_damage_predictive_uniaxial_stress;
            mMaxStressPlasticityBranch = plasticity_predictive_uniaxial_stress;
        }
        // if (!mMinDetected && min_indicator) {
        //     mMinStressDamageBranch = isotropic_damage_predictive_uniaxial_stress;
        //     mMinStressPlasticityBranch = plasticity_predictive_uniaxial_stress;
        // }
        mMaxDetected = max_indicator;
        mMinDetected = min_indicator;
        Vector previous_stresses = ZeroVector(2);
        const Vector& r_aux_stresses = mPreviousStresses;
        previous_stresses[1] = predictive_uniaxial_stress;
        previous_stresses[0] = r_aux_stresses[1];
        mPreviousStresses = previous_stresses;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::FinalizeMaterialResponseIsotropicDamageModel(
        double& rIsotropicDamagePredictiveUniaxialStress,
        Vector& rIsotropicDamageStressVector,
        ConstitutiveLaw::Parameters& values_damage_component
)
{
    // Integrate Stress Damage
    const Flags& r_constitutive_law_options = values_damage_component.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = values_damage_component.GetStrainVector();

    // Elastic Matrix
    Matrix& r_constitutive_matrix = values_damage_component.GetConstitutiveMatrix();
    mpIsotropicDamageConstitutiveLaw->CalculateValue(values_damage_component, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain( values_damage_component, r_strain_vector);
    }

    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // Converged values
    double isotropic_damage_threshold = this->GetIsotropicDamageThreshold();
    double damage = this->GetDamage();

    // S0 = C:(E-E0) + S0
    array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
    this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_stress_vector);

    // Initialize Plastic Parameters
    double uniaxial_stress;
    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, values_damage_component);
    // rIsotropicDamagePredictiveUniaxialStress = uniaxial_stress;

    double fatigue_reduction_factor = mFatigueReductionFactor;
    uniaxial_stress /= fatigue_reduction_factor;  // Fatigue contribution

    const double F = uniaxial_stress - isotropic_damage_threshold;

    if (F >= threshold_tolerance) { // Plastic case
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_damage_component.GetElementGeometry());
        // This routine updates the PredictiveStress to verify the yield surf
        TConstLawIntegratorType::IntegrateStressVector(
            predictive_stress_vector,
            uniaxial_stress, damage,
            isotropic_damage_threshold, values_damage_component,
            characteristic_length);
        mDamage = damage;
        mIsotropicDamageThreshold = uniaxial_stress;

    } else {
        predictive_stress_vector *= (1.0 - mDamage);
    }

    mReferenceDamage = (mCyclesToFailure < 1.0e5) ? mDamage : mReferenceDamage; //Updating the reference value for ULCF or LCF case
    rIsotropicDamagePredictiveUniaxialStress = (1.0 - mReferenceDamage) * uniaxial_stress;
    rIsotropicDamageStressVector = predictive_stress_vector;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePlastitictyModel(
        double& rPlasticityPredictiveUniaxialStress,
        Vector& rPlasticityStressVector,
        ConstitutiveLaw::Parameters& values_plasticity_component
)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = values_plasticity_component.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = values_plasticity_component.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = values_plasticity_component.GetConstitutiveMatrix();

    // Integrate Stress plasticity
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_plasticity_component.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain( values_plasticity_component, r_strain_vector);
    }

    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // We compute the stress
    // Elastic Matrix
    this->CalculateElasticMatrix(r_constitutive_matrix, values_plasticity_component);

    // We get some variables
    double plasticity_threshold = this->GetPlasticityThreshold();
    double plastic_dissipation = this->GetPlasticDissipation();
    Vector plastic_strain = this->GetPlasticStrain();

    array_1d<double, VoigtSize> predictive_stress_vector;
    if (r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW)) {
        noalias(predictive_stress_vector) = values_plasticity_component.GetStressVector();
    } else {
        // Spred = r_constitutive_matrix:(E-Ep) + S0
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_stress_vector);
    }

    // Initialize Plastic Parameters
    double uniaxial_stress = 0.0, plastic_denominator = 0.0;
    array_1d<double, VoigtSize> f_flux = ZeroVector(VoigtSize); // DF/DS
    array_1d<double, VoigtSize> g_flux = ZeroVector(VoigtSize); // DG/DS
    array_1d<double, VoigtSize> plastic_strain_increment = ZeroVector(VoigtSize);
    double fatigue_reduction_factor = mFatigueReductionFactor;

    const double F = GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::CalculatePlasticParameters(
        predictive_stress_vector, r_strain_vector, uniaxial_stress,
        plasticity_threshold, plastic_denominator, f_flux, g_flux,
        plastic_dissipation, plastic_strain_increment,
        r_constitutive_matrix, values_plasticity_component, characteristic_length,
        plastic_strain, fatigue_reduction_factor);

    const double yield_criteria = uniaxial_stress - plasticity_threshold;

    if (yield_criteria > std::abs(1.0e-4 * plasticity_threshold)) { // Plastic case
        // While loop backward euler
        /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
        GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::IntegrateStressVector(
            predictive_stress_vector, r_strain_vector, uniaxial_stress,
            plasticity_threshold, plastic_denominator, f_flux, g_flux,
            plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, plastic_strain, values_plasticity_component,
            characteristic_length, fatigue_reduction_factor);
        BaseType::CalculateElasticMatrix( r_constitutive_matrix, values_plasticity_component);
    }
    mReferencePlasticStrain = (mCyclesToFailure < 1.0e5) ? plastic_strain : mReferencePlasticStrain;    //ReferencePlasticStrain will be updated always except of those cases where HCF is expected,
                                                                                                        //where cycles should be counted consistently to guarantee fred soft evolution
    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(prod(r_constitutive_matrix, r_strain_vector - mReferencePlasticStrain), r_strain_vector, rPlasticityPredictiveUniaxialStress, values_plasticity_component);

    mPlasticDissipation = plastic_dissipation;
    mPlasticStrain = plastic_strain;
    mPlasticityThreshold = plasticity_threshold;
    rPlasticityStressVector = predictive_stress_vector;
}
/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::Has(const Variable<bool>& rThisVariable)
{
    if (rThisVariable == CYCLE_INDICATOR) {
        return true;
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::Has(const Variable<int>& rThisVariable)
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
bool UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
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
    } else if (rThisVariable == PREVIOUS_CYCLE_PLASTIC_DISSIPATION) {
        return true;
    } else if (rThisVariable == DAMAGE) {
        return true;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;

}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::SetValue(
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
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::SetValue(
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
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
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
        mFatigueLimit = rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        mPreviousCycleTime = rValue;
    } else if (rThisVariable == CYCLE_PERIOD) {
        mPeriod = rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        mPreviousCycleDamage = rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE_PLASTIC_DISSIPATION) {
        mPreviousCyclePlasticDissipation = rValue;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else if (rThisVariable == DAMAGE) {
        mDamage = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        mPlasticStrain = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::GetValue(
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
int& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::GetValue(
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
double& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == FATIGUE_REDUCTION_FACTOR) {
        rValue = mFatigueReductionFactor;
    } else if (rThisVariable == WOHLER_STRESS) {
        rValue = mWohlerStress;
    } else if (rThisVariable == CYCLES_TO_FAILURE) {
        rValue = mCyclesToFailure;
    } else if (rThisVariable == REVERSION_FACTOR_RELATIVE_ERROR) {
        rValue = mReversionFactorRelativeError;
    } else if (rThisVariable == MAX_STRESS_RELATIVE_ERROR) {
        rValue = mMaxStressRelativeError;
    } else if (rThisVariable == MAX_STRESS) {
        rValue = mMaxStress;
    } else if (rThisVariable == THRESHOLD_STRESS) {
        rValue = mFatigueLimit;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        rValue = mPreviousCycleTime;
    } else if (rThisVariable == CYCLE_PERIOD) {
        rValue = mPeriod;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        rValue = mPreviousCycleDamage;
    } else if (rThisVariable == PREVIOUS_CYCLE_PLASTIC_DISSIPATION) {
        rValue = mPreviousCyclePlasticDissipation;
    } else if (rThisVariable == DAMAGE) {
        rValue = mDamage;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
        return rValue;
    } else {
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        rValue = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    // const Vector& r_strain_vector = rParameterValues.GetStrainVector();

    // Recalculation to obtain the serial_strain_matrix and store the value
    // const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    // Flags& r_flags = rParameterValues.GetOptions();

    // Previous flags saved
    // const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);


    auto& r_material_properties = rParameterValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
    const auto& r_props_plasticity_cl = *(it_cl_begin + 1);
    ConstitutiveLaw::Parameters values_damage_component  = rParameterValues;
    ConstitutiveLaw::Parameters values_plasticity_component = rParameterValues;
    values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
    values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);

    if (rThisVariable == UNIAXIAL_STRESS) {
        // Get Values to compute the constitutive law:
        // Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        // const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        // const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        // r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        // r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // Calculate the stress vector
        CalculateMaterialResponseCauchy(rParameterValues);

        // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector& r_integrated_stress_vector = rParameterValues.GetStressVector();
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( r_integrated_stress_vector, r_strain_vector, rValue, rParameterValues);

        return rValue;

    } else if (rThisVariable == UNIAXIAL_STRESS_HCF) {
        this->CalculateMaterialResponseIsotropicDamageModelModel(values_damage_component);
        Vector& r_strain_vector = values_damage_component.GetStrainVector();
        Vector& r_integrated_stress_vector = values_damage_component.GetStressVector();
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( r_integrated_stress_vector, r_strain_vector, rValue, values_damage_component);
        return rValue;
    } else if (rThisVariable == UNIAXIAL_STRESS_ULCF) {
        this->CalculateMaterialResponsePlasticityModel(values_plasticity_component);
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector& r_integrated_stress_vector = values_plasticity_component.GetStressVector();
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( r_integrated_stress_vector, r_strain_vector, rValue, values_plasticity_component);
        return rValue;
    // } else if (rThisVariable == EQUIVALENT_PLASTIC_STRAIN) {
        // // Get Values to compute the constitutive law:
        // Flags& r_flags = rParameterValues.GetOptions();

        // // Previous flags saved
        // const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        // const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        // r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        // r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // // Calculate the stress vector
        // CalculateMaterialResponseCauchy(rParameterValues);
        // const Vector& r_stress_vector = rParameterValues.GetStressVector();

        // // Previous flags restored
        // r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        // r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );

        // // Compute the equivalent plastic strain
        // double uniaxial_stress;
        // this->CalculateValue(rParameterValues, UNIAXIAL_STRESS, uniaxial_stress);
        // TConstLawIntegratorType::CalculateEquivalentPlasticStrain(r_stress_vector, uniaxial_stress, mPlasticStrain, 0.0, rParameterValues, rValue);
        // return rValue;
    } else {
        return this->GetValue(rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We do some special operations for constitutive matrices
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, false);

         // We compute the constitutive matrix
        if (rThisVariable == CONSTITUTIVE_MATRIX) {
            this->CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_PK2) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        }

        noalias(rValue) = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    } else if (rThisVariable == DEFORMATION_GRADIENT) { // TODO: Make in the future modifications for take into account different layers combinations
        noalias(rValue) = rParameterValues.GetDeformationGradientF();
    } else if (rThisVariable == CAUCHY_STRESS_TENSOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        this->CalculateMaterialResponseCauchy(rParameterValues);
        rValue = MathUtils<double>::StressVectorToTensor(rParameterValues.GetStressVector());

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
        return rValue;
    } else {
        Matrix aux_value;
        Properties material_properties  = rParameterValues.GetMaterialProperties();
        Properties& r_prop = material_properties.GetSubProperties(0);

        rValue.clear();
        rParameterValues.SetMaterialProperties(r_prop);
        mpIsotropicDamageConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += mIsotrpicDamageVolumetricParticipation * aux_value;

        r_prop = material_properties.GetSubProperties(1);
        rParameterValues.SetMaterialProperties(r_prop);
        mpPlasticityConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mIsotrpicDamageVolumetricParticipation) * aux_value;

        // Reset properties
        rParameterValues.SetMaterialProperties(material_properties);
    }
    return(rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto r_props_isotropic_damage_cl = *(it_cl_begin);
    const auto r_props_plasticity_cl  = *(it_cl_begin + 1);

    KRATOS_ERROR_IF_NOT(r_props_isotropic_damage_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    KRATOS_ERROR_IF_NOT(r_props_plasticity_cl.Has(CONSTITUTIVE_LAW))  << "No constitutive law set" << std::endl;

    mpIsotropicDamageConstitutiveLaw = r_props_isotropic_damage_cl[CONSTITUTIVE_LAW]->Clone();
    mpPlasticityConstitutiveLaw  = r_props_plasticity_cl[CONSTITUTIVE_LAW]->Clone();
    mpIsotropicDamageConstitutiveLaw->InitializeMaterial(r_props_isotropic_damage_cl, rElementGeometry, rShapeFunctionsValues);
    mpPlasticityConstitutiveLaw ->InitializeMaterial(r_props_plasticity_cl, rElementGeometry, rShapeFunctionsValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}
/***********************************************************************************/
/***********************************************************************************/
template class UnifiedFatigueRuleOfMixturesLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;

} // namespace Kratos
