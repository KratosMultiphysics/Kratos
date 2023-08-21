//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo
//
//

# pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/auxiliary_files/cl_integrators/high_cycle_fatigue_law_integrator.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class HCFDataContainer
 * @ingroup ConstitutiveLawsApplication
 * @brief Defining all the methods and variables required in fatigue simulations
 * @details Give access to methods and restores all member variables in fatigue simulations
 * @author Alireza Taherzadeh-Fard, Alejandro Cornejo
 */
class HCFDataContainer
{
ConstitutiveLaw::Parameters rValues = ConstitutiveLaw::Parameters();
public:

    struct FatigueVariables {
        double MaxStress = 0.0;
        double MinStress = 0.0;
        bool MaxIndicator = false;
        bool MinIndicator = false;
        Vector PreviousStresses = ZeroVector(2);
        double FatigueReductionFactor = 1.0;
        double ReversionFactorRelativeError = 0.0;
        double MaxStressRelativeError = 0.0;
        unsigned int GlobalNumberOfCycles = 1;
        unsigned int LocalNumberOfCycles = 1;
        double B0 = 0.0;
        double PreviousMaxStress = 0.0;
        double PreviousMinStress = 0.0;
        double WohlerStress = 1.0;
        double Sth = 0.0;
        double CyclesToFailure = 0.0;
        bool NewCycle = false;
        double Alphat = 0.0;
        double PreviousReversionFactor = 0.0;
        double ReversionFactor = 0.0;
        bool AdnvanceStrategyApplied;
        bool DamageActivation;
    };

	HCFDataContainer()
    {};

    // Defining fatigue methods

    void CalculateSminAndSmax(const double CurrentStress,
                            HCFDataContainer::FatigueVariables &rFatigueVariables
                            )
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
        // mPreviousStresses = rFatigueVariables.PreviousStresses;
    }

    double CalculateTensionOrCompressionIdentifier(const Vector& rStressVector)
    {
        return HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(rStressVector);
    }

    double CalculateReversionFactor(const double MaxStress, const double MinStress)
    {
        return HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(MaxStress, MinStress);
    }

    void CalculateFatigueParameters(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables)
    {
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(rFatigueVariables.MaxStress,
                                                                    rFatigueVariables.ReversionFactor,
                                                                    rMaterialParameters,
                                                                    rFatigueVariables.B0,
                                                                    rFatigueVariables.Sth,
                                                                    rFatigueVariables.Alphat,
                                                                    rFatigueVariables.CyclesToFailure);
    }

    void CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables)
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

    /**
     * @brief This method initializes all the values
     * in the FatigueVariables
     */
    void InitializeFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables)
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
        rFatigueVariables.B0 = mFatigueReductionParameter;
        rFatigueVariables.PreviousMaxStress = mPreviousMaxStress;
        rFatigueVariables.PreviousMinStress = mPreviousMinStress;
        rFatigueVariables.WohlerStress = mWohlerStress;
        rFatigueVariables.NewCycle = false;
        rFatigueVariables.Sth = mThresholdStress;
        rFatigueVariables.CyclesToFailure = mCyclesToFailure;
    }

    /**
     * @brief This method updates all member variables
     */
    void UpdateFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables)
    {
        mMaxStress = rFatigueVariables.MaxStress;
        mMinStress = rFatigueVariables.MinStress;
        mMaxDetected = rFatigueVariables.MaxIndicator;
        mMinDetected = rFatigueVariables.MinIndicator;
        mNumberOfCyclesGlobal = rFatigueVariables.GlobalNumberOfCycles;
        mNumberOfCyclesLocal = rFatigueVariables.LocalNumberOfCycles;
        mNewCycleIndicator = rFatigueVariables.NewCycle;
        mFatigueReductionParameter = rFatigueVariables.B0;
        mPreviousMaxStress = rFatigueVariables.PreviousMaxStress;
        mPreviousMinStress = rFatigueVariables.PreviousMinStress;
        mPreviousStresses = rFatigueVariables.PreviousStresses;
        mFatigueReductionFactor = rFatigueVariables.FatigueReductionFactor;
        mWohlerStress = rFatigueVariables.WohlerStress;
        mThresholdStress = rFatigueVariables.Sth;
        mReversionFactorRelativeError = rFatigueVariables.ReversionFactorRelativeError;
        mMaxStressRelativeError = rFatigueVariables.MaxStressRelativeError;
    }

    /**
     * @brief This method computes fatigue-related quantities
     */
    void FinalizeSolutionStep(HCFDataContainer::FatigueVariables &rFatigueVariables,
                            const Properties& rMaterialProperties,
                            const ProcessInfo& rCurrentProcessInfo,
                            ConstitutiveLaw::StressVectorType stress_vector,
                            double uniaxial_stress)
    {
        double sign_factor = CalculateTensionOrCompressionIdentifier(stress_vector);
        uniaxial_stress *= sign_factor;

        CalculateSminAndSmax(uniaxial_stress, rFatigueVariables);

        rFatigueVariables.AdnvanceStrategyApplied = rCurrentProcessInfo[ADVANCE_STRATEGY_APPLIED];
        rFatigueVariables.DamageActivation = rCurrentProcessInfo[DAMAGE_ACTIVATION];

        if (rFatigueVariables.MaxIndicator && rFatigueVariables.MinIndicator) {
            rFatigueVariables.PreviousReversionFactor = CalculateReversionFactor(rFatigueVariables.PreviousMaxStress, rFatigueVariables.PreviousMinStress);
            rFatigueVariables.ReversionFactor = CalculateReversionFactor(rFatigueVariables.MaxStress, rFatigueVariables.MinStress);

            CalculateFatigueParameters(rMaterialProperties, rFatigueVariables);

            double betaf = rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];

            if (std::abs(rFatigueVariables.MinStress) < 0.001) {
                rFatigueVariables.ReversionFactorRelativeError = std::abs(rFatigueVariables.ReversionFactor - rFatigueVariables.PreviousReversionFactor);
            } else {
                rFatigueVariables.ReversionFactorRelativeError = std::abs((rFatigueVariables.ReversionFactor - rFatigueVariables.PreviousReversionFactor) / rFatigueVariables.ReversionFactor);
            }
            rFatigueVariables.MaxStressRelativeError = std::abs((rFatigueVariables.MaxStress - rFatigueVariables.PreviousMaxStress) / rFatigueVariables.MaxStress);

            if (!rFatigueVariables.DamageActivation && rFatigueVariables.GlobalNumberOfCycles > 2 && !rFatigueVariables.AdnvanceStrategyApplied && (rFatigueVariables.ReversionFactorRelativeError > 0.001 || rFatigueVariables.MaxStressRelativeError > 0.001)) {
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
        if (rFatigueVariables.AdnvanceStrategyApplied) {
        rFatigueVariables.ReversionFactor = CalculateReversionFactor(rFatigueVariables.MaxStress, rFatigueVariables.MinStress);

        CalculateFatigueParameters(rMaterialProperties, rFatigueVariables);

        CalculateFatigueReductionFactorAndWohlerStress(rMaterialProperties, rFatigueVariables);
        }
    }

    // Defining fatigue variables

    ///@name Member Variables
    ///@{
    double mMaxStress = 0.0;
    double mMinStress = 0.0;
    bool mMaxDetected = false;
    bool mMinDetected = false;
    Vector mPreviousStresses = ZeroVector(2);
    double mFatigueReductionFactor = 1.0;
    double mReversionFactorRelativeError = 0.0;
    double mMaxStressRelativeError = 0.0;
    unsigned int mNumberOfCyclesGlobal = 1;
    unsigned int mNumberOfCyclesLocal = 1;
    double mFatigueReductionParameter = 0.0;
    double mPreviousMaxStress = 0.0;
    double mPreviousMinStress = 0.0;
    double mWohlerStress = 1.0;
    double mThresholdStress = 0.0;
    double mCyclesToFailure = 0.0;
    bool mNewCycleIndicator = false;
    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}

friend class Serializer;

void save(Serializer& rSerializer) const
{
    rSerializer.save("MaxStress",mMaxStress);
    rSerializer.save("MinStress",mMinStress);
    rSerializer.save("MaxDetected",mMaxDetected);
    rSerializer.save("MinDetected",mMinDetected);
    rSerializer.save("PreviousStresses",mPreviousStresses);
    rSerializer.save("FatigueReductionFactor",mFatigueReductionFactor);
    rSerializer.save("ReversionFactorRelativeError",mReversionFactorRelativeError);
    rSerializer.save("MaxStressRelativeError",mMaxStressRelativeError);
    rSerializer.save("NumberOfCyclesGlobal",mNumberOfCyclesGlobal);
    rSerializer.save("NumberOfCyclesLocal",mNumberOfCyclesLocal);
    rSerializer.save("FatigueReductionParameter",mFatigueReductionParameter);
    rSerializer.save("PreviousMaxStress",mPreviousMaxStress);
    rSerializer.save("PreviousMinStress",mPreviousMinStress);
    rSerializer.save("WohlerStress",mWohlerStress);
    rSerializer.save("ThresholdStress",mThresholdStress);
    rSerializer.save("CyclesToFailure",mCyclesToFailure);
    rSerializer.save("NewCycleIndicator",mNewCycleIndicator);
}

void load(Serializer& rSerializer)
{
    rSerializer.load("MaxStress",mMaxStress);
    rSerializer.load("MinStress",mMinStress);
    rSerializer.load("MaxDetected",mMaxDetected);
    rSerializer.load("MinDetected",mMinDetected);
    rSerializer.load("PreviousStresses",mPreviousStresses);
    rSerializer.load("FatigueReductionFactor",mFatigueReductionFactor);
    rSerializer.load("ReversionFactorRelativeError",mReversionFactorRelativeError);
    rSerializer.load("MaxStressRelativeError",mMaxStressRelativeError);
    rSerializer.load("NumberOfCyclesGlobal",mNumberOfCyclesGlobal);
    rSerializer.load("NumberOfCyclesLocal",mNumberOfCyclesLocal);
    rSerializer.load("FatigueReductionParameter",mFatigueReductionParameter);
    rSerializer.load("PreviousMaxStress",mPreviousMaxStress);
    rSerializer.load("PreviousMinStress",mPreviousMinStress);
    rSerializer.load("WohlerStress",mWohlerStress);
    rSerializer.load("ThresholdStress",mThresholdStress);
    rSerializer.load("CyclesToFailure",mCyclesToFailure);
    rSerializer.load("NewCycleIndicator",mNewCycleIndicator);
}

}; // class
///@}

///@} addtogroup block

} // namespace Kratos
