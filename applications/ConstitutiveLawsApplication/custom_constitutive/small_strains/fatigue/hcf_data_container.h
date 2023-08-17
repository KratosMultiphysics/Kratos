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
	HCFDataContainer()
    {};

    // Defining fatigue methods

    void CalculateSminAndSmax(const double CurrentStress,
                            double& rMaximumStress,
                            double& rMinimumStress,
                            const Vector& PreviousStresses,
                            bool& rMaxIndicator,
                            bool& rMinIndicator)
    {
        HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(CurrentStress,
                                                                            rMaximumStress,
                                                                            rMinimumStress,
                                                                            PreviousStresses,
                                                                            rMaxIndicator,
                                                                            rMinIndicator);
    }

    double CalculateTensionOrCompressionIdentifier(const Vector& rStressVector)
    {
        return HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(rStressVector);
    }

    double CalculateReversionFactor(const double MaxStress, const double MinStress)
    {
        return HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(MaxStress, MinStress);
    }

    void CalculateFatigueParameters(const double MaxStress,
                                    double ReversionFactor,
                                    const Properties& rMaterialParameters,
                                    double& rB0,
                                    double& rSth,
                                    double& rAlphat,
                                    double& rN_f)
    {
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(MaxStress,
                                                                    ReversionFactor,
                                                                    rMaterialParameters,
                                                                    rB0,
                                                                    rSth,
                                                                    rAlphat,
                                                                    rN_f);
    }

    void CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters,
                                                        const double MaxStress,
                                                        unsigned int LocalNumberOfCycles,
                                                        unsigned int GlobalNumberOfCycles,
                                                        const double B0,
                                                        const double Sth,
                                                        const double Alphat,
                                                        double& rFatigueReductionFactor,
                                                        double& rWohlerStress)
    {
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(rMaterialParameters,
                                                                                        MaxStress,
                                                                                        LocalNumberOfCycles,
                                                                                        GlobalNumberOfCycles,
                                                                                        B0,
                                                                                        Sth,
                                                                                        Alphat,
                                                                                        rFatigueReductionFactor,
                                                                                        rWohlerStress);
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
