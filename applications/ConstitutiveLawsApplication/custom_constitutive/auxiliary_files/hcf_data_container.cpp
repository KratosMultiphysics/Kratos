// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:                 BSD License
//                               license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo
//
// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/auxiliary_files/hcf_data_container.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/high_cycle_fatigue_law_integrator.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include <limits>

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::CalculateSminAndSmax(const double CurrentStress,
                        HCFDataContainer::FatigueVariables &rFatigueVariables)
{
    HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(CurrentStress,
                                                                        rFatigueVariables.MaxStress,
                                                                        rFatigueVariables.MinStress,
                                                                        rFatigueVariables.PreviousStresses,
                                                                        rFatigueVariables.MaxIndicator,
                                                                        rFatigueVariables.MinIndicator);
    // KRATOS_WATCH(mPreviousStresses);
    // KRATOS_WATCH(CurrentStress);
    const Vector& r_aux_stresses = mPreviousStresses;
    rFatigueVariables.PreviousStresses[1] = CurrentStress;
    rFatigueVariables.PreviousStresses[0] = r_aux_stresses[1];
    // KRATOS_WATCH(rFatigueVariables.PreviousStresses);
}

/***********************************************************************************/
/***********************************************************************************/

double HCFDataContainer::CalculateTensionOrCompressionIdentifier(const Vector& rStressVector)
{
    return HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(rStressVector);
}

/***********************************************************************************/
/***********************************************************************************/

double HCFDataContainer::CalculateReversionFactor(const double MaxStress, const double MinStress)
{
    return HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(MaxStress, MinStress);
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::CalculateFatigueParameters(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables, const Variable<double>& rVariable)
{
    // const Vector& r_fatigue_coefficients = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS];

    double ultimate_stress = 0.0;

    if (rVariable == YIELD_STRESS) {
        ultimate_stress = rMaterialParameters.Has(YIELD_STRESS) ? rMaterialParameters[YIELD_STRESS] : rMaterialParameters[YIELD_STRESS_TENSION];
    } else {
        ultimate_stress = rMaterialParameters[rVariable];
    }
    const double yield_stress = ultimate_stress;

    //These variables have been defined following the model described by S. Oller et al. in A continuum mechanics model for mechanical fatigue analysis (2005), equation 13 on page 184.
    // const double Se = r_fatigue_coefficients[0] * ultimate_stress;
    // const double STHR1 = r_fatigue_coefficients[1];
    // const double STHR2 = r_fatigue_coefficients[2];
    // const double ALFAF = r_fatigue_coefficients[3];
    // const double BETAF = r_fatigue_coefficients[4];
    // const double AUXR1 = r_fatigue_coefficients[5];
    // const double AUXR2 = r_fatigue_coefficients[6];
    // const double c_factor = r_fatigue_coefficients[7];

    const double Se = rFatigueVariables.hcf_coefficients[0] * ultimate_stress;
    const double STHR1 = rFatigueVariables.hcf_coefficients[1];
    const double STHR2 = rFatigueVariables.hcf_coefficients[2];
    const double ALFAF = rFatigueVariables.hcf_coefficients[3];
    const double BETAF = rFatigueVariables.hcf_coefficients[4];
    const double AUXR1 = rFatigueVariables.hcf_coefficients[5];
    const double AUXR2 = rFatigueVariables.hcf_coefficients[6];
    const double c_factor = rFatigueVariables.hcf_coefficients[7];

    if (std::abs(rFatigueVariables.ReversionFactor) < 1.0) {
        rFatigueVariables.Sth = Se + (rFatigueVariables.UltimateStress - Se) * std::pow((0.5 + 0.5 * rFatigueVariables.ReversionFactor), STHR1);
        rFatigueVariables.Alphat = ALFAF + (0.5 + 0.5 * rFatigueVariables.ReversionFactor) * AUXR1;
    } else {
        rFatigueVariables.Sth = Se + (rFatigueVariables.UltimateStress - Se) * std::pow((0.5 + 0.5 / rFatigueVariables.ReversionFactor), STHR2);
        rFatigueVariables.Alphat = ALFAF - (0.5 + 0.5 / rFatigueVariables.ReversionFactor) * AUXR2;
    }

    const double square_betaf = std::pow(BETAF, 2.0);
    double MaxStressRD = (1 - rFatigueVariables.ReferenceDamage) * rFatigueVariables.MaxStress;

    if (MaxStressRD > rFatigueVariables.Sth && MaxStressRD <= ultimate_stress) {
        rFatigueVariables.CyclesToFailure = std::pow(10.0,std::pow(-std::log((MaxStressRD - rFatigueVariables.Sth) / (ultimate_stress - rFatigueVariables.Sth))/rFatigueVariables.Alphat,(1.0/BETAF)));
        rFatigueVariables.B0 = -(std::log(MaxStressRD / ultimate_stress) / std::pow((std::log10(rFatigueVariables.CyclesToFailure)), c_factor * square_betaf));

        // if (softening_type == curve_by_points) {
        //     rFatigueVariables.CyclesToFailure = std::pow(rFatigueVariables.CyclesToFailure, std::pow(std::log(rFatigueVariables.MaxStress / yield_stress) / std::log(rFatigueVariables.MaxStress / ultimate_stress), 1.0 / square_betaf));
        // }
        
    } else {
        rFatigueVariables.CyclesToFailure = std::numeric_limits<double>::infinity();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables, const Variable<double>& rVariable)
{
    const double BETAF = rFatigueVariables.hcf_coefficients[4];
    const double c_factor = rFatigueVariables.hcf_coefficients[7];

    double ultimate_stress = 0.0;

    if (rVariable == YIELD_STRESS) {
        ultimate_stress = rMaterialParameters.Has(YIELD_STRESS) ? rMaterialParameters[YIELD_STRESS] : rMaterialParameters[YIELD_STRESS_TENSION];
    } else {
        ultimate_stress = rMaterialParameters[rVariable];
    }

    // The calculation is prepared to update the rN_f value when using a softening curve which initiates with hardening.
    // The jump in the advance in time process is done in these cases to the Syield rather to Sult.
    const int softening_type = rMaterialParameters[SOFTENING_TYPE];
    const int curve_by_points = static_cast<int>(SofteningType::CurveFittingDamage);
    if (softening_type == curve_by_points) {
        const Vector& stress_damage_curve = rMaterialParameters[STRESS_DAMAGE_CURVE]; //Integrated_stress points of the fitting curve
        const SizeType curve_points = stress_damage_curve.size() - 1;

        ultimate_stress = 0.0;
        for (IndexType i = 1; i <= curve_points; ++i) {
            ultimate_stress = std::max(ultimate_stress, stress_damage_curve[i-1]);
        }
    }

    if (rFatigueVariables.GlobalNumberOfCycles > 2){
        rFatigueVariables.WohlerStress = (rFatigueVariables.Sth + (ultimate_stress - rFatigueVariables.Sth) * std::exp(-rFatigueVariables.Alphat * (std::pow(std::log10(static_cast<double>(rFatigueVariables.LocalNumberOfCycles)), BETAF)))) / ultimate_stress;
    }

    double MaxStressRD = (1 - rFatigueVariables.ReferenceDamage) * rFatigueVariables.MaxStress;
    if (MaxStressRD > rFatigueVariables.Sth) {
        if (MaxStressRD <= ultimate_stress) {
        rFatigueVariables.FatigueReductionFactor = std::min(rFatigueVariables.FatigueReductionFactor, std::exp(-rFatigueVariables.B0 * std::pow(std::log10(static_cast<double>(rFatigueVariables.LocalNumberOfCycles)), c_factor * (BETAF * BETAF))));
        } else {
            if (!rFatigueVariables.is_initiated) {
                rFatigueVariables.nc_initiation = rFatigueVariables.GlobalNumberOfCycles;
                rFatigueVariables.fred_initiation = rFatigueVariables.FatigueReductionFactor;
                rFatigueVariables.is_initiated = true;
            }
            // double etha = (rFatigueVariables.GlobalNumberOfCycles - rFatigueVariables.nc_initiation) / 10;
            // rFatigueVariables.FatigueReductionFactor = (-rFatigueVariables.fred_initiation) * etha + rFatigueVariables.fred_initiation;

            double alfa_fred = rFatigueVariables.hcf_coefficients[8];
            double beta_fred = rFatigueVariables.hcf_coefficients[9];
            double etha = rFatigueVariables.GlobalNumberOfCycles - rFatigueVariables.nc_initiation;
            
            // rFatigueVariables.FatigueReductionFactor = rFatigueVariables.fred_initiation - ((1/(alfa_fred*(1-beta_fred)))*(std::pow(rFatigueVariables.GlobalNumberOfCycles,(1-beta_fred))-std::pow(rFatigueVariables.nc_initiation,(1-beta_fred))));
            rFatigueVariables.FatigueReductionFactor = rFatigueVariables.fred_initiation - ((1/(alfa_fred*(1-beta_fred)))*(std::pow(etha,(1-beta_fred))));
        }
        rFatigueVariables.FatigueReductionFactor = (rFatigueVariables.FatigueReductionFactor < 0.01) ? 0.01 : rFatigueVariables.FatigueReductionFactor;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::InitializeFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables, ConstitutiveLaw::Parameters& rValues)
{
    rFatigueVariables.MaxStress = mMaxStress;
    rFatigueVariables.MinStress = mMinStress;
    rFatigueVariables.MaxIndicator = mMaxDetected;
    rFatigueVariables.MinIndicator = mMinDetected;
    rFatigueVariables.PreviousStresses = mPreviousStresses;
    rFatigueVariables.FatigueReductionFactor = mFatigueReductionFactor;
    rFatigueVariables.ReversionFactorRelativeError = mReversionFactorRelativeError;
    rFatigueVariables.MaxStressRelativeError = mMaxStressRelativeError;
    rFatigueVariables.GlobalNumberOfCycles = mNumberOfCyclesGlobal;
    rFatigueVariables.LocalNumberOfCycles = mNumberOfCyclesLocal;
    rFatigueVariables.PreviousMaxStress = mPreviousMaxStress;
    rFatigueVariables.PreviousMinStress = mPreviousMinStress;
    rFatigueVariables.WohlerStress = mWohlerStress;
    rFatigueVariables.NewCycle = false;
    rFatigueVariables.Sth = mThresholdStress;
    rFatigueVariables.ReferenceDamage = mReferenceDamage;
    rFatigueVariables.CyclesToFailure = mCyclesToFailure;
    rFatigueVariables.DamageActivation = rValues.GetProcessInfo()[DAMAGE_ACTIVATION];
    rFatigueVariables.AdvanceStrategyApplied = rValues.GetProcessInfo()[ADVANCE_STRATEGY_APPLIED];
    rFatigueVariables.current_load_type = rValues.GetProcessInfo()[CURRENT_LOAD_TYPE];
    rFatigueVariables.new_model_part = rValues.GetProcessInfo()[NEW_MODEL_PART];
    rFatigueVariables.is_initiated = mIsInitiated;
    rFatigueVariables.nc_initiation = mNcInitiation;
    rFatigueVariables.fred_initiation = mFredInitiation;
}

/***********************************************************************************/
/***********************************************************************************/

double HCFDataContainer::UltimateStressDamage(const Properties& rMaterialParameters)
{
    return HighCycleFatigueLawIntegrator<6>::UltimateStressDamage(rMaterialParameters);
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::UpdateFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables)
{
    mMaxStress = rFatigueVariables.MaxStress;
    mMinStress = rFatigueVariables.MinStress;
    mMaxDetected = rFatigueVariables.MaxIndicator;
    mMinDetected = rFatigueVariables.MinIndicator;
    mNumberOfCyclesGlobal = rFatigueVariables.GlobalNumberOfCycles;
    mNumberOfCyclesLocal = rFatigueVariables.LocalNumberOfCycles;
    mNewCycleIndicator = rFatigueVariables.NewCycle;
    mPreviousMaxStress = rFatigueVariables.PreviousMaxStress;
    mPreviousMinStress = rFatigueVariables.PreviousMinStress;
    mPreviousStresses = rFatigueVariables.PreviousStresses;
    mFatigueReductionFactor = rFatigueVariables.FatigueReductionFactor;
    mWohlerStress = rFatigueVariables.WohlerStress;
    mThresholdStress = rFatigueVariables.Sth;
    mReferenceDamage = rFatigueVariables.ReferenceDamage;
    mReversionFactorRelativeError = rFatigueVariables.ReversionFactorRelativeError;
    mMaxStressRelativeError = rFatigueVariables.MaxStressRelativeError;
    mIsInitiated = rFatigueVariables.is_initiated;
    mNcInitiation = rFatigueVariables.nc_initiation;
    mFredInitiation = rFatigueVariables.fred_initiation;
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::FinalizeSolutionStep(HCFDataContainer::FatigueVariables &rFatigueVariables,
                        const Properties& rMaterialProperties,
                        const ProcessInfo& rCurrentProcessInfo,
                        ConstitutiveLaw::StressVectorType stress_vector,
                        double uniaxial_stress,
                        double damage,
                        double threshold,
                        const Variable<double>& rVariable,
                        const Variable<bool>& rMethodVariable)
{
    if (rFatigueVariables.current_load_type) {

    double ultimate_stress = 0.0;

    if (rVariable == YIELD_STRESS) {
        ultimate_stress = rMaterialProperties.Has(YIELD_STRESS) ? rMaterialProperties[YIELD_STRESS] : rMaterialProperties[YIELD_STRESS_TENSION];
    } else {
        ultimate_stress = rMaterialProperties[rVariable];
    }
    const double yield_stress = ultimate_stress;

    double sign_factor = 0;

    if (rMethodVariable == DEFAULT_METHOD) {
        sign_factor = CalculateTensionOrCompressionIdentifier(stress_vector);
    } else if (rMethodVariable == USER_DEFINED_METHOD) {
        double aux_sum = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets(stress_vector[4]) + AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets(stress_vector[5]);
        double resultant_shear_stress = std::sqrt(std::pow(stress_vector[4],2.0)+std::pow(stress_vector[5],2.0));
        const double pre_indicator = aux_sum / resultant_shear_stress;

        if (pre_indicator < 0.7) {
            sign_factor = 1.0;
        } else {
            sign_factor = -1.0;
        }
    }
    uniaxial_stress *= sign_factor;

    if (rFatigueVariables.new_model_part) {
        rFatigueVariables.MaxIndicator = false;
        rFatigueVariables.MinIndicator = false;
        rFatigueVariables.ReferenceDamage = damage;
    }

    CalculateSminAndSmax(uniaxial_stress, rFatigueVariables);
    
    rFatigueVariables.AdvanceStrategyApplied = rCurrentProcessInfo.Has(ADVANCE_STRATEGY_APPLIED) ? rCurrentProcessInfo[ADVANCE_STRATEGY_APPLIED] : false;
    rFatigueVariables.DamageActivation = rCurrentProcessInfo.Has(DAMAGE_ACTIVATION) ? rCurrentProcessInfo[DAMAGE_ACTIVATION] : false;
    const double c_factor = rFatigueVariables.hcf_coefficients[7];

    if (rFatigueVariables.MaxIndicator && rFatigueVariables.MinIndicator) {
        rFatigueVariables.PreviousReversionFactor = CalculateReversionFactor(rFatigueVariables.PreviousMaxStress, rFatigueVariables.PreviousMinStress);
        rFatigueVariables.ReversionFactor = CalculateReversionFactor(rFatigueVariables.MaxStress, rFatigueVariables.MinStress);

        mReversionFactor = rFatigueVariables.ReversionFactor;

        // if (rFatigueVariables.LocalNumberOfCycles < 2.0 && std::abs(uniaxial_stress) > ultimate_stress) {
        //     rFatigueVariables.ReferenceDamage = damage;
        // }

        CalculateFatigueParameters(rMaterialProperties, rFatigueVariables, rVariable);

        // const double betaf = rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];
        const double betaf = rFatigueVariables.hcf_coefficients[4];

        const double square_betaf = std::pow(betaf, 2.0);
        double MaxStressRD = (1 - rFatigueVariables.ReferenceDamage) * rFatigueVariables.MaxStress;
        threshold *= (1 - damage);

        if (threshold - MaxStressRD >= -1 * tolerance) {
            if (MaxStressRD > rFatigueVariables.Sth) {
                rFatigueVariables.CyclesToFailure = std::pow(rFatigueVariables.CyclesToFailure, std::pow(std::log(MaxStressRD / threshold) / std::log(MaxStressRD / ultimate_stress), 1.0 / (c_factor * square_betaf)));
            }
        }

        if (MaxStressRD > rFatigueVariables.Sth) { //&& MaxStressRD < ultimate_stress
            if (std::abs(rFatigueVariables.MinStress) < tolerance) {
                rFatigueVariables.ReversionFactorRelativeError = std::abs(rFatigueVariables.ReversionFactor - rFatigueVariables.PreviousReversionFactor);
            } else {
                rFatigueVariables.ReversionFactorRelativeError = std::abs((rFatigueVariables.ReversionFactor - rFatigueVariables.PreviousReversionFactor) / rFatigueVariables.ReversionFactor);
            }
            rFatigueVariables.MaxStressRelativeError = std::abs((rFatigueVariables.MaxStress - rFatigueVariables.PreviousMaxStress) / rFatigueVariables.MaxStress);

            if (rFatigueVariables.GlobalNumberOfCycles > 2 && !rFatigueVariables.AdvanceStrategyApplied && (rFatigueVariables.ReversionFactorRelativeError > tolerance || rFatigueVariables.MaxStressRelativeError > tolerance)) {
                rFatigueVariables.LocalNumberOfCycles = std::trunc(std::pow(10, std::pow(-(std::log(rFatigueVariables.FatigueReductionFactor) / rFatigueVariables.B0), 1.0 / (betaf * betaf)))) + 1;
            }

        // if (MaxStressRD > rFatigueVariables.Sth) { //&& MaxStressRD < ultimate_stress
            rFatigueVariables.GlobalNumberOfCycles++;
            rFatigueVariables.LocalNumberOfCycles++;
            // rFatigueVariables.NewCycle = true;
            // rFatigueVariables.MaxIndicator = false;
            // rFatigueVariables.MinIndicator = false;
            rFatigueVariables.PreviousMaxStress = rFatigueVariables.MaxStress;
            rFatigueVariables.PreviousMinStress = rFatigueVariables.MinStress;
            mCyclesToFailure = rFatigueVariables.CyclesToFailure;

            CalculateFatigueReductionFactorAndWohlerStress(rMaterialProperties, rFatigueVariables, rVariable);
        }
        rFatigueVariables.MaxIndicator = false;
        rFatigueVariables.MinIndicator = false;
        rFatigueVariables.NewCycle = true;
    }
    if (rFatigueVariables.AdvanceStrategyApplied) {
        double MaxStressRD = (1 - rFatigueVariables.ReferenceDamage) * rFatigueVariables.MaxStress;
        if (MaxStressRD > rFatigueVariables.Sth) {
            rFatigueVariables.GlobalNumberOfCycles += mIncrementInNumberOfCycles; 
            rFatigueVariables.LocalNumberOfCycles += mIncrementInNumberOfCycles; 
        }
        rFatigueVariables.ReversionFactor = CalculateReversionFactor(rFatigueVariables.MaxStress, rFatigueVariables.MinStress);

        CalculateFatigueParameters(rMaterialProperties, rFatigueVariables, rVariable);

        // double MaxStressRD = (1 - rFatigueVariables.ReferenceDamage) * rFatigueVariables.MaxStress;

        if (MaxStressRD > rFatigueVariables.Sth) {
            CalculateFatigueReductionFactorAndWohlerStress(rMaterialProperties, rFatigueVariables, rVariable);
        }
        mAITControlParameter = rFatigueVariables.LocalNumberOfCycles;
        // KRATOS_WATCH(mAITControlParameter);
        // rFatigueVariables.PreviousStresses = ZeroVector(2);
        rFatigueVariables.MaxIndicator = false;
        rFatigueVariables.MinIndicator = false;
    }

    // mAITControlParameter += 1.0;
}
}

/***********************************************************************************/
/***********************************************************************************/

int HCFDataContainer::Check(
    const Properties& rMaterialProperties)
{
    // Check if input parameters are well defined
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not a defined vector" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO)) << "POISSON_RATIO is not a defined vector" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRACTURE_ENERGY)) << "FRACTURE_ENERGY is not a defined vector" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS)) << "YIELD_STRESS is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SOFTENING_TYPE)) << "SOFTENING_TYPE is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(HIGH_CYCLE_FATIGUE_COEFFICIENTS)) << "HIGH_CYCLE_FATIGUE_COEFFICIENTS is not a defined vector" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS].size() != 7) << "HIGH_CYCLE_FATIGUE_COEFFICIENTS badly defined" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][0] < 0.0 || rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][0] > 1.0) << "Se/Su is not properly defined." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][1] < 0.0) << "SthR1 is not properly defined." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][3] < std::max(0.0,-rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][5])) << "alfaf is not properly defined." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4] < 1.0) << "bataf is not properly defined." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][5] < -rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][3]) << "AUXR1 is not properly defined." << std::endl;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
