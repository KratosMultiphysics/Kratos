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

    // /**
    //  * @brief This method Calculates and extrapolates the
    //  * free energy indicator over the nodes
    //  * @param pNodeNormalizedFreeEnergyVector the indicator container
    //  */
    // void NormalizedFreeEnergyExtrapolation(NodeNormalizedFreeEnergy *pNodeNormalizedFreeEnergyVector);
    
    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable);

    /**
     * @brief This method computes the free energy 
     * @param rStrainVector The strain vector
     * @param rStressVector The stress vector
     * @param Damage The damage variable
     * @param rMatProps The material properties
     * @param rGeometry The geometry of the element
     */
    void CyclePeriodPerIntegrationPoint(bool& rCycleFound);

    /**
     * @brief This method computes the free energy 
     * @param rStrainVector The strain vector
     * @param rStressVector The stress vector
     * @param Damage The damage variable
     * @param rMatProps The material properties
     * @param rGeometry The geometry of the element
     */
    void StableConditionForAdvancingStrategy(bool& rAdvancingStrategy, bool DamageIndicator);

    /**
     * @brief This method computes the free energy 
     * @param rStrainVector The strain vector
     * @param rStressVector The stress vector
     * @param Damage The damage variable
     * @param rMatProps The material properties
     * @param rGeometry The geometry of the element
     */
    void TimeIncrement(double& rIncrement);

    /**
     * @brief This method computes the free energy 
     * @param rStrainVector The strain vector
     * @param rStressVector The stress vector
     * @param Damage The damage variable
     * @param rMatProps The material properties
     * @param rGeometry The geometry of the element
     */
    void TimeAndCyclesUpdate(double Increment);

protected:
    // Member Variables
    ModelPart& mrModelPart;                     // The model part to compute
    Parameters mThisParameters;

}; // Class AdvanceInTimeStrategyHighCycleFatigueProcess

} // namespace Kratos

#endif /* KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS defined */