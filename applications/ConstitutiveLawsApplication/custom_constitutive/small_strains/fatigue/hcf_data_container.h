//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

# pragma once

// System includes
#include <atomic>

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
 * @ingroup StructuralMechanicsApplication
 * @brief Define the initial state of the material in terms of initial stress/strain/F
 * @details Storages the information regarding initial stresses/strains/F
 * @author Alejandro Cornejo
 */
class KRATOS_API(KRATOS_CORE) HCFDataContainer
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
        HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(rStressVector);
    }

    double CalculateReversionFactor(const double MaxStress, const double MinStress)
    {
        HighCycleFatigueLawIntegrator<6>::CalculateReversionFactor(MaxStress, MinStress);
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
    double mMaxStress;
    double mMinStress;
    bool mMaxDetected;
    bool mMinDetected;
    double mFatigueReductionFactor;
    double mReversionFactorRelativeError;
    double mMaxStressRelativeError;
    unsigned int mNumberOfCyclesGlobal;
    unsigned int mNumberOfCyclesLocal;
    double mFatigueReductionParameter;
    double mPreviousMaxStress;
    double mPreviousMinStress;
    double mWohlerStress;
    double mThresholdStress;
    double mCyclesToFailure;
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
