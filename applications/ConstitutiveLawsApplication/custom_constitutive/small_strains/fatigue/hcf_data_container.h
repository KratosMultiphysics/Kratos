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
//
//

# pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/auxiliary_files/cl_integrators/high_cycle_fatigue_law_integrator.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

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
 * @author Alireza Taherzadeh-Fard
 */
class HCFDataContainer
{

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
        mCyclesToFailure = rFatigueVariables.CyclesToFailure;
        mReversionFactorRelativeError = rFatigueVariables.ReversionFactorRelativeError;
        mMaxStressRelativeError = rFatigueVariables.MaxStressRelativeError;
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

}; // class
///@}

///@} addtogroup block

} // namespace Kratos
