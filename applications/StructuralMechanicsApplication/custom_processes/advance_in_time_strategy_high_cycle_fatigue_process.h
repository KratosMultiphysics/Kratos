//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sergio Jimenez Reyes
//

#if !defined(KRATOS_ADVANCE_IN_TIME_STRATEGY_HIGH_CYCLE_FATIGUE_PROCESS)
#define KRATOS_ADVANCE_IN_TIME_STRATEGY_HIGH_CYCLE_FATIGUE_PROCESS

#include "processes/process.h"
#include "includes/model_part.h"
#include "structural_mechanics_application_variables.h"

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
 * @class AdvanceInTimeStrategyHighCycleFatigueProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This class determines the advance in time to be performed for a regular cyclic load for the high cycle fatigue CL
 * @author Sergio Jimenez
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)AdvanceInTimeStrategyHighCycleFatigueProcess : public Process
{
    ///@name Type Definitions
    ///@{
    typedef ModelPart::ElementsContainerType ElementsArrayType;

    ///@}
    ///@name  Enum's
    ///@{

protected:
    

public:
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AdvanceInTimeStrategyHighCycleFatigueProcess);

    /**
     * @brief This is the default constructor (double)
     * @param rModelPart The model part to be used
     * @param ThisParameters The input parameters
     */
	AdvanceInTimeStrategyHighCycleFatigueProcess(ModelPart& rModelPart, Parameters ThisParameters);

    // Destructor
    ~AdvanceInTimeStrategyHighCycleFatigueProcess() override = default;

    void Execute() override;

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
    void TimeAndCyclesUpdate(double Increment);

protected:
    // Member Variables
    ModelPart& mrModelPart;                     // The model part to compute
    Parameters mThisParameters;

}; // Class AdvanceInTimeStrategyHighCycleFatigueProcess

} // namespace Kratos

#endif /* KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS defined */