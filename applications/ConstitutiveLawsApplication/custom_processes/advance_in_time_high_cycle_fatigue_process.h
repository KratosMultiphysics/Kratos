// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sergio Jimenez Reyes
//

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class AdvanceInTimeHighCycleFatigueProcess
 * @ingroup ConstitutiveLawsApplication
 * @brief This class determines the advance in time to be performed for a regular cyclic load for the high cycle fatigue CL
 * @author Sergio Jimenez
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION)AdvanceInTimeHighCycleFatigueProcess : public Process
{
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name  Enum's
    ///@{

protected:


public:
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AdvanceInTimeHighCycleFatigueProcess);

    /**
     * @brief This is the default constructor (double)
     * @param rModelPart The model part to be used
     * @param ThisParameters The input parameters
     */
	AdvanceInTimeHighCycleFatigueProcess(ModelPart& rModelPart, Parameters ThisParameters);

    // Destructor
    ~AdvanceInTimeHighCycleFatigueProcess() override = default;

    void Execute() override;


    /**
     * @brief This method checks which kind of load is being applied in the model in the current time
     */
    void CyclicLoad();

    /**
     * @brief This method computes the cycle time period per integration point
     * @param rCycleFound Bool variable indicating that a cycle has overcome at some integration point
     */
    void CyclePeriodPerIntegrationPoint(bool& rCycleFound);

    /**
     * @brief This method stablishes if stable conditions have been reached for initiating the advance strategy
     * @param rAdvancingStrategy Bool variable indicating weather advancing strategy will start or not
     * @param DamageIndicator Bool variable indicating that damage has iniciated at some point
     */
    void StableConditionForAdvancingStrategy(bool& rAdvancingStrategy, bool DamageIndicator);

    /**
     * @brief This method computes the time increment to be applied as an output of the advancing strategy
     * @param rIncrement Double variable corresponding to the time increment to apply
     */
    void TimeIncrement(double& rIncrement);

    /**
     * @brief This method properly applies the time increment in terms of cycle increment to all the integration points of the model
     * @param Increment Time increment to apply along the model
     */
    void TimeAndCyclesUpdate(const double Increment);

protected:
    // Member Variables
    ModelPart& mrModelPart;                     // The model part to compute
    Parameters mThisParameters;

}; // Class AdvanceInTimeHighCycleFatigueProcess

} // namespace Kratos