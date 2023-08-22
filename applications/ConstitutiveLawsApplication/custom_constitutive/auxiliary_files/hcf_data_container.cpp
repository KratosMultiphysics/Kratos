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
    const Vector& r_aux_stresses = mPreviousStresses;
    rFatigueVariables.PreviousStresses[1] = CurrentStress;
    rFatigueVariables.PreviousStresses[0] = r_aux_stresses[1];
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

void HCFDataContainer::CalculateFatigueParameters(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables)
{
    HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(rFatigueVariables.MaxStress,
                                                                rFatigueVariables.ReversionFactor,
                                                                rMaterialParameters,
                                                                rFatigueVariables.B0,
                                                                rFatigueVariables.Sth,
                                                                rFatigueVariables.Alphat,
                                                                rFatigueVariables.CyclesToFailure);
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables)
{
    HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(rMaterialParameters,
                                                                                    rFatigueVariables.MaxStress,
                                                                                    rFatigueVariables.LocalNumberOfCycles,
                                                                                    rFatigueVariables.GlobalNumberOfCycles,
                                                                                    rFatigueVariables.B0,
                                                                                    rFatigueVariables.Sth,
                                                                                    rFatigueVariables.Alphat,
                                                                                    rFatigueVariables.FatigueReductionFactor,
                                                                                    rFatigueVariables.WohlerStress);
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::InitializeFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables)
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
    rFatigueVariables.CyclesToFailure = mCyclesToFailure;
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
    mReversionFactorRelativeError = rFatigueVariables.ReversionFactorRelativeError;
    mMaxStressRelativeError = rFatigueVariables.MaxStressRelativeError;
}

/***********************************************************************************/
/***********************************************************************************/

void HCFDataContainer::FinalizeSolutionStep(HCFDataContainer::FatigueVariables &rFatigueVariables,
                        const Properties& rMaterialProperties,
                        const ProcessInfo& rCurrentProcessInfo,
                        ConstitutiveLaw::StressVectorType stress_vector,
                        double uniaxial_stress)
{
    double sign_factor = CalculateTensionOrCompressionIdentifier(stress_vector);
    uniaxial_stress *= sign_factor;

    CalculateSminAndSmax(uniaxial_stress, rFatigueVariables);

    rFatigueVariables.AdvanceStrategyApplied = rCurrentProcessInfo[ADVANCE_STRATEGY_APPLIED];
    rFatigueVariables.DamageActivation = rCurrentProcessInfo[DAMAGE_ACTIVATION];

    if (rFatigueVariables.MaxIndicator && rFatigueVariables.MinIndicator) {
        rFatigueVariables.PreviousReversionFactor = CalculateReversionFactor(rFatigueVariables.PreviousMaxStress, rFatigueVariables.PreviousMinStress);
        rFatigueVariables.ReversionFactor = CalculateReversionFactor(rFatigueVariables.MaxStress, rFatigueVariables.MinStress);

        CalculateFatigueParameters(rMaterialProperties, rFatigueVariables);

        const double betaf = rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];

        if (std::abs(rFatigueVariables.MinStress) < tolerance) {
            rFatigueVariables.ReversionFactorRelativeError = std::abs(rFatigueVariables.ReversionFactor - rFatigueVariables.PreviousReversionFactor);
        } else {
            rFatigueVariables.ReversionFactorRelativeError = std::abs((rFatigueVariables.ReversionFactor - rFatigueVariables.PreviousReversionFactor) / rFatigueVariables.ReversionFactor);
        }
        rFatigueVariables.MaxStressRelativeError = std::abs((rFatigueVariables.MaxStress - rFatigueVariables.PreviousMaxStress) / rFatigueVariables.MaxStress);

        if (!rFatigueVariables.DamageActivation && rFatigueVariables.GlobalNumberOfCycles > 2 && !rFatigueVariables.AdvanceStrategyApplied && (rFatigueVariables.ReversionFactorRelativeError > tolerance || rFatigueVariables.MaxStressRelativeError > 0.001)) {
            rFatigueVariables.LocalNumberOfCycles = std::trunc(std::pow(10, std::pow(-(std::log(rFatigueVariables.FatigueReductionFactor) / rFatigueVariables.B0), 1.0 / (betaf * betaf)))) + 1;
        }

        rFatigueVariables.GlobalNumberOfCycles++;
        rFatigueVariables.LocalNumberOfCycles++;
        rFatigueVariables.NewCycle = true;
        rFatigueVariables.MaxIndicator = false;
        rFatigueVariables.MinIndicator = false;
        rFatigueVariables.PreviousMaxStress = rFatigueVariables.MaxStress;
        rFatigueVariables.PreviousMinStress = rFatigueVariables.MinStress;
        mCyclesToFailure = rFatigueVariables.CyclesToFailure;

        CalculateFatigueReductionFactorAndWohlerStress(rMaterialProperties, rFatigueVariables);
    }
    if (rFatigueVariables.AdvanceStrategyApplied) {
    rFatigueVariables.ReversionFactor = CalculateReversionFactor(rFatigueVariables.MaxStress, rFatigueVariables.MinStress);

    CalculateFatigueParameters(rMaterialProperties, rFatigueVariables);

    CalculateFatigueReductionFactorAndWohlerStress(rMaterialProperties, rFatigueVariables);
    }
}

    /***********************************************************************************/
    /***********************************************************************************/
} // Namespace Kratos
