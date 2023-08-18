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
        double max_stress = 0.0;
        double min_stress = 0.0;
        bool max_indicator = false;
        bool min_indicator = false;
        Vector previous_stresses = ZeroVector(2);
        double fatigue_reduction_factor = 1.0;
        double reversion_factor_relative_error = 0.0;
        double max_stress_relative_error = 0.0;
        unsigned int global_number_of_cycles = 1;
        unsigned int local_number_of_cycles = 1;
        double B0 = 0.0;
        double previous_max_stress = 0.0;
        double previous_min_stress = 0.0;
        double wohler_stress = 1.0;
        double s_th = 0.0;
        double cycles_to_failure = 0.0;
        bool new_cycle = false;
        double alphat = 0.0;
        double previous_reversion_factor = 0.0;
        double reversion_factor = 0.0;
    };

	HCFDataContainer()
    {};

    // Defining fatigue methods

    void CalculateSminAndSmax(const double CurrentStress,
                            HCFDataContainer::FatigueVariables &rFatigueVariables
                            )
    {
        HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(CurrentStress,
                                                                            rFatigueVariables.max_stress,
                                                                            rFatigueVariables.min_stress,
                                                                            rFatigueVariables.previous_stresses,
                                                                            rFatigueVariables.max_indicator,
                                                                            rFatigueVariables.min_indicator);
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
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueParameters(rFatigueVariables.max_stress,
                                                                    rFatigueVariables.reversion_factor,
                                                                    rMaterialParameters,
                                                                    rFatigueVariables.B0,
                                                                    rFatigueVariables.s_th,
                                                                    rFatigueVariables.alphat,
                                                                    rFatigueVariables.cycles_to_failure);
    }

    void CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables)
    {
        HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactorAndWohlerStress(rMaterialParameters,
                                                                                        rFatigueVariables.max_stress,
                                                                                        rFatigueVariables.local_number_of_cycles,
                                                                                        rFatigueVariables.global_number_of_cycles,
                                                                                        rFatigueVariables.B0,
                                                                                        rFatigueVariables.s_th,
                                                                                        rFatigueVariables.alphat,
                                                                                        rFatigueVariables.fatigue_reduction_factor,
                                                                                        rFatigueVariables.wohler_stress);
    }

    /**
     * @brief This method initializes all the values
     * in the FatigueVariables
     */
    void InitializeFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables)
    {
        rFatigueVariables.max_stress = mMaxStress;
        rFatigueVariables.min_stress = mMinStress;
        rFatigueVariables.max_indicator = mMaxDetected;
        rFatigueVariables.min_indicator = mMinDetected;
        rFatigueVariables.previous_stresses = mPreviousStresses;
        rFatigueVariables.fatigue_reduction_factor = mFatigueReductionFactor;
        rFatigueVariables.reversion_factor_relative_error = mReversionFactorRelativeError;
        rFatigueVariables.max_stress_relative_error = mMaxStressRelativeError;
        rFatigueVariables.global_number_of_cycles = mNumberOfCyclesGlobal;
        rFatigueVariables.local_number_of_cycles = mNumberOfCyclesLocal;
        rFatigueVariables.B0 = mFatigueReductionParameter;
        rFatigueVariables.previous_max_stress = mPreviousMaxStress;
        rFatigueVariables.previous_min_stress = mPreviousMinStress;
        rFatigueVariables.wohler_stress = mWohlerStress;
        rFatigueVariables.new_cycle = false;
        rFatigueVariables.s_th = mThresholdStress;
        rFatigueVariables.cycles_to_failure = mCyclesToFailure;
    }

    /**
     * @brief This method updates all member variables
     */
    void UpdateFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables)
    {
        mMaxStress = rFatigueVariables.max_stress;
        mMinStress = rFatigueVariables.min_stress;
        mMaxDetected = rFatigueVariables.max_indicator;
        mMinDetected = rFatigueVariables.min_indicator;
        mNumberOfCyclesGlobal = rFatigueVariables.global_number_of_cycles;
        mNumberOfCyclesLocal = rFatigueVariables.local_number_of_cycles;
        mNewCycleIndicator = rFatigueVariables.new_cycle;
        mFatigueReductionParameter = rFatigueVariables.B0;
        mPreviousMaxStress = rFatigueVariables.previous_max_stress;
        mPreviousMinStress = rFatigueVariables.previous_min_stress;
        mFatigueReductionFactor = rFatigueVariables.fatigue_reduction_factor;
        mWohlerStress = rFatigueVariables.wohler_stress;
        mThresholdStress = rFatigueVariables.s_th;
        mCyclesToFailure = rFatigueVariables.cycles_to_failure;
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
