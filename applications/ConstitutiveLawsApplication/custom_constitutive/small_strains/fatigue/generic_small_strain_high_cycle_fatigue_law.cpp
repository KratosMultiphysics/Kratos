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
#include "custom_utilities/advanced_constitutive_law_utilities.h"

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
    double reference_damage = mReferenceDamage;
    double previous_max_stress = mPreviousMaxStress;
    double previous_min_stress = mPreviousMinStress;
    double max_stress;
    double min_stress;
    double uniaxial_stress;
    double ultimate_stress;
    // double initial_threshold = mInitialTherhold;
    // double first_cycle_relaxation_factor = mFirstCycleRelaxationFactor;
    double relaxation_factor = mRelaxationFactor;
    bool first_nonlinearity = mFirstNonlinearity;
    bool first_max_indicator = mFirstMaxDetected;
    bool first_min_indicator = mFirstMinDetected;
    bool max_indicator = mMaxDetected;
    bool min_indicator = mMinDetected;
    double fatigue_reduction_factor = mFatigueReductionFactor;
    double reversion_factor_relative_error = mReversionFactorRelativeError;
    double max_stress_relative_error = mMaxStressRelativeError;
    unsigned int global_number_of_cycles = mNumberOfCyclesGlobal;
    unsigned int local_number_of_cycles = mNumberOfCyclesLocal;
    double B0 = mFatigueReductionParameter;
    bool new_cycle = false;
    double s_th = mThresholdStress;
    double cycles_to_failure = mCyclesToFailure;
    bool advance_strategy_applied = rValues.GetProcessInfo()[ADVANCE_STRATEGY_APPLIED];
    const bool new_model_part = rValues.GetProcessInfo()[NEW_MODEL_PART];
    int time = rValues.GetProcessInfo()[TIME];
    const int load_increments_per_cycle = rValues.GetProcessInfo()[LOAD_INCREMENTS_PER_CYCLE];
    int time_offset = rValues.GetProcessInfo()[NEW_MODEL_PART_START_TIME];

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();
    
    // We compute the residual stress
    array_1d<double, VoigtSize> residual_stress_vector;
    double uniaxial_residual_stress;
    
    residual_stress_vector = ZeroVector(VoigtSize);
    this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(residual_stress_vector);
    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(residual_stress_vector, r_strain_vector, uniaxial_residual_stress, rValues);
    uniaxial_residual_stress = (rValues.GetElementGeometry().Has(INITIAL_STRESS_VECTOR)) ? uniaxial_residual_stress : 0.0;

    if (((time - time_offset) % load_increments_per_cycle) == 0 && (time - time_offset) > 0) {

        if (std::abs(mMaxStress) > std::abs(mMinStress)){
            uniaxial_stress = std::abs(mMaxStress);
        } else {
            uniaxial_stress = std::abs(mMinStress);
        }

        double threshold = this->GetThreshold() * (1 - this->GetDamage());

        if (uniaxial_stress / threshold > 1.0 && first_nonlinearity) {
            reference_damage = this->GetDamage();
            mFirstNonlinearity = false;  
        }

        max_stress = (1 - reference_damage) * mMaxStress;
        min_stress = (1 - reference_damage) * mMinStress;

        const double previous_reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(previous_max_stress, previous_min_stress);
        double reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(max_stress, min_stress);  

        if (std::abs(min_stress) < 0.001) {
            reversion_factor_relative_error = std::abs(reversion_factor - previous_reversion_factor);
        } else {
            reversion_factor_relative_error = std::abs((reversion_factor - previous_reversion_factor) / reversion_factor);
        }
        max_stress_relative_error = std::abs((max_stress - previous_max_stress) / max_stress);

        // const double betaf = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];
        // const double fatigue_reduction_factor_smoothness = rValues.GetMaterialProperties()[HIGH_CYCLE_FATIGUE_COEFFICIENTS][7];
        // if (global_number_of_cycles > 2 && !advance_strategy_applied && (reversion_factor_relative_error > 0.001 || max_stress_relative_error > 0.001)) {
        //     local_number_of_cycles = std::trunc(std::pow(10, std::pow(-(std::log(fatigue_reduction_factor) / B0), 1.0 / (fatigue_reduction_factor_smoothness* betaf * betaf)))) + 1;
        // }
        
        global_number_of_cycles++;
        local_number_of_cycles++;

        previous_max_stress = max_stress;
        previous_min_stress = min_stress;

        residual_stress_vector *= relaxation_factor;
        const double residual_stress_sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(residual_stress_vector);
        
        if (residual_stress_sign_factor > 0.0) {
           uniaxial_residual_stress *= relaxation_factor;
        } else {
           uniaxial_residual_stress = 0.0;
        }

        double alphat;
        HighCycleFatigueLawIntegrator<6>::CalculateUltimateStress(ultimate_stress, rValues.GetMaterialProperties());

        HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(
            max_stress,
            uniaxial_residual_stress,
            ultimate_stress,
            reversion_factor,
            threshold,
            local_number_of_cycles,
            rValues.GetMaterialProperties(),
            rValues.GetElementGeometry(),
            B0,
            s_th,
            alphat,
            cycles_to_failure);
             
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactor(rValues.GetMaterialProperties(),
                                                                                        max_stress,
                                                                                        local_number_of_cycles,
                                                                                        B0,
                                                                                        s_th,
                                                                                        fatigue_reduction_factor);


        new_cycle = true;
        first_max_indicator = true;
        first_min_indicator = true;
        mCyclesToFailure = cycles_to_failure;
        mPreviousCycleDamage = this->GetDamage();
        mFirstNonlinearity = false;  
    }
    
    if (advance_strategy_applied) {

        max_stress = (1 - reference_damage) * mMaxStress;
        min_stress = (1 - reference_damage) * mMinStress;
               
        residual_stress_vector *= relaxation_factor;
        const double residual_stress_sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(residual_stress_vector);

        if (residual_stress_sign_factor > 0.0) {
           uniaxial_residual_stress *= relaxation_factor;
        } else {
           uniaxial_residual_stress = 0.0;
        }
               
        double alphat;
        double threshold = this->GetThreshold() * (1 - this->GetDamage());
        double reversion_factor = HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(max_stress, min_stress);
        HighCycleFatigueLawIntegrator<6>::CalculateUltimateStress(ultimate_stress, rValues.GetMaterialProperties());

        HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(
            max_stress,
            uniaxial_residual_stress,
            ultimate_stress,
            reversion_factor,
            threshold,
            local_number_of_cycles,
            rValues.GetMaterialProperties(),
            rValues.GetElementGeometry(),
            B0,
            s_th,
            alphat,
            cycles_to_failure);
      
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactor(rValues.GetMaterialProperties(),
                                                                                        max_stress,
                                                                                        local_number_of_cycles,
                                                                                        B0,
                                                                                        s_th,
                                                                                        fatigue_reduction_factor);

        first_max_indicator = true;
        first_min_indicator = true;
      
    }
   
    if (new_model_part) {
        mCyclesToFailure = std::numeric_limits<double>::infinity();
        mReferenceDamage = 0.0;
        mFirstNonlinearity = true;
    }

    max_indicator = false;
    min_indicator = false;

    if (uniaxial_residual_stress > mInitialTherhold){
        this->SetThreshold(uniaxial_residual_stress);
        this->SetInitialThreshold(uniaxial_residual_stress);
    }

    mNumberOfCyclesGlobal = global_number_of_cycles;
    mNumberOfCyclesLocal = local_number_of_cycles;
    mReversionFactorRelativeError = reversion_factor_relative_error;
    mMaxStressRelativeError = max_stress_relative_error;
    mFirstMaxDetected = first_max_indicator;
    mFirstMinDetected = first_min_indicator;
    mMaxDetected = max_indicator;
    mMinDetected = min_indicator;
    mFatigueReductionParameter = B0;
    mPreviousMaxStress = previous_max_stress;
    mPreviousMinStress = previous_min_stress;
    mFatigueReductionFactor = fatigue_reduction_factor;
    mNewCycleIndicator = new_cycle;
    mThresholdStress = s_th;
    mRelaxationFactor = relaxation_factor;
    mReferenceDamage = reference_damage;  
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
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }

    array_1d<double, VoigtSize> predictive_residual_stress_vector;
    array_1d<double, VoigtSize> predictive_stress_vector;
    array_1d<double, VoigtSize> aux_predictive_stress_vector;
    // We compute the stress
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
                
        // Elastic Matrix
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // Converged values
        double damage = this->GetDamage();
        double threshold = this->GetThreshold() * (1 - damage);
        unsigned int local_number_of_cycles = mNumberOfCyclesLocal;
        double first_cycle_relaxation_factor = mFirstCycleRelaxationFactor;
        double relaxation_factor = mRelaxationFactor;
        bool first_nonlinearity = mFirstNonlinearity;

        // S00       
        noalias(predictive_residual_stress_vector) = ZeroVector(VoigtSize);
        double uniaxial_residual_stress;
        this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_residual_stress_vector);
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_residual_stress_vector, r_strain_vector, uniaxial_residual_stress, rValues);

        // S0 = C:(E-E0) + S00
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);
        noalias(aux_predictive_stress_vector) = predictive_stress_vector;
        double nominal_uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, nominal_uniaxial_stress, rValues);

        predictive_stress_vector += relaxation_factor * predictive_residual_stress_vector;
        
        // Initialize Plastic Parameters
        double fatigue_reduction_factor = mFatigueReductionFactor;
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        const double plasticity_tol = 0.05;
        // int time = rValues.GetProcessInfo()[TIME];
        
        double F = uniaxial_stress - threshold;

        if (F > threshold_tolerance && first_nonlinearity && (1.0 / fatigue_reduction_factor - 1.0) <= plasticity_tol) {
            HighCycleFatigueLawIntegrator<6>::CalculateRelaxationFactor(nominal_uniaxial_stress,
                                                                        uniaxial_residual_stress,
                                                                        threshold,          
                                                                        local_number_of_cycles,
                                                                        first_cycle_relaxation_factor,
                                                                        relaxation_factor,
                                                                        first_nonlinearity);
            mFirstCycleRelaxationFactor = first_cycle_relaxation_factor;
            mRelaxationFactor = relaxation_factor;
            
            // if (relaxation_factor < 0.3319){
            //     KRATOS_WATCH("RelaxationFactor")
            //     KRATOS_WATCH(relaxation_factor)
            // }
        }
        
        threshold /= (1 - damage);
        uniaxial_stress /= fatigue_reduction_factor;  // Fatigue contribution
        F = uniaxial_stress - threshold;

        double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(predictive_stress_vector);        

        if (F <= threshold_tolerance || (sign_factor < 0.0 && !first_nonlinearity)) { // Elastic case
            noalias(r_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                r_constitutive_matrix *= (1.0 - damage);
                this->SetStressVector(r_integrated_stress_vector);
                rValues.SetStressVector(r_integrated_stress_vector);
            }
        } else { // Damage case
            const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
            double initial_threshold = mInitialTherhold;
            // This routine updates the PredictiveStress to verify the yield surf
            TConstLawIntegratorType::IntegrateStressVector(
                predictive_stress_vector,
                uniaxial_stress,
                damage,
                threshold,
                rValues,
                characteristic_length,
                initial_threshold);          
            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);
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
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    BaseType::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);

    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TConstLawIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    this->SetInitialThreshold(initial_threshold);
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
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }
    
    array_1d<double, VoigtSize> predictive_residual_stress_vector;
    array_1d<double, VoigtSize> predictive_stress_vector;
    // We compute the stress
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {

        // Elastic Matrix
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // Converged values
        double threshold = this->GetThreshold();
        double damage = this->GetDamage();
        // double first_cycle_relaxation_factor = mFirstCycleRelaxationFactor;
        double relaxation_factor = mRelaxationFactor;
        // unsigned int local_number_of_cycles = mNumberOfCyclesLocal;

        // S00
        noalias(predictive_residual_stress_vector) = ZeroVector(VoigtSize);
        this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_residual_stress_vector);

        // S0 = C:(E-E0) + S00
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);
        predictive_stress_vector += relaxation_factor * predictive_residual_stress_vector;
        
        // Initialize Plastic Parameters
        double fatigue_reduction_factor = mFatigueReductionFactor;
        double uniaxial_stress;

        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(predictive_stress_vector);
        uniaxial_stress *= sign_factor;
        double max_stress = mMaxStress;
        double min_stress = mMinStress;
        bool first_max_indicator = mFirstMaxDetected;
        bool first_min_indicator = mFirstMinDetected;
        bool max_indicator = mMaxDetected;
        bool min_indicator = mMinDetected;
        bool first_nonlinearity = mFirstNonlinearity;

        HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(
            uniaxial_stress,
            mPreviousStresses,
            max_stress,
            min_stress,
            first_max_indicator,
            first_min_indicator,
            max_indicator,
            min_indicator);    
    
        mMaxStress = max_stress;
        mMinStress = min_stress;
        mFirstMaxDetected = first_max_indicator;
        mFirstMaxDetected = first_min_indicator;
        mMaxDetected = max_indicator;
        mMinDetected = min_indicator;
                
        uniaxial_stress *= sign_factor;
        uniaxial_stress /= fatigue_reduction_factor;  // Fatigue contribution

        const double F = uniaxial_stress - threshold;

        if (F > threshold_tolerance && (sign_factor > 0.0 || first_nonlinearity)) {
                const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
                double initial_threshold = mInitialTherhold;
                // This routine updates the PredictiveStress to verify the yield surface
                TConstLawIntegratorType::IntegrateStressVector(
                    predictive_stress_vector,
                    uniaxial_stress,
                    damage,
                    threshold,
                    rValues,
                    characteristic_length,
                    initial_threshold);
                this->SetDamage(damage);
                this->SetThreshold(uniaxial_stress);
        } else {
            predictive_stress_vector *= (1.0 - damage);
            TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);
        }

        Vector previous_stresses = ZeroVector(2);
        const Vector& r_aux_stresses = mPreviousStresses;
        previous_stresses[1] = this->CalculateValue(rValues, UNIAXIAL_STRESS, previous_stresses[1]) * sign_factor / (1.0 - damage);
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
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        return true;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        return true;
    } else if (rThisVariable == CYCLE_PERIOD) {
        return true;
    } else if (rThisVariable == INFINITY_YIELD_STRESS) {
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
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        mPreviousCycleDamage = rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        mThresholdStress = rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        mPreviousCycleTime = rValue;
    } else if (rThisVariable == CYCLE_PERIOD) {
        mPeriod = rValue;
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
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        rValue = mPreviousCycleDamage;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        rValue = mPreviousCycleTime;
    } else if (rThisVariable == CYCLE_PERIOD) {
        rValue = mPeriod;
    } else if (rThisVariable == INFINITY_YIELD_STRESS) {
        rValue = mMinStress/mMaxStress;
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
    if (this->Has(rThisVariable)){
        return this->GetValue(rThisVariable, rValue);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
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