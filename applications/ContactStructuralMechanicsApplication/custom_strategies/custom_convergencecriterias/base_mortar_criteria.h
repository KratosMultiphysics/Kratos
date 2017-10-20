// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BASE_MORTAR_CRITERIA_H)
#define  KRATOS_BASE_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_processes/aalm_adapt_penalty_value_process.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
///@{
    
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

/** @brief Custom convergence criteria for the mortar condition 
 */
template<class TSparseSpace, class TDenseSpace>
class BaseMortarConvergenceCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( BaseMortarConvergenceCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace                              SparseSpaceType;

    typedef typename BaseType::TDataType                    TDataType;

    typedef typename BaseType::DofsArrayType            DofsArrayType;

    typedef typename BaseType::TSystemMatrixType    TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType    TSystemVectorType;
    
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    
    typedef ModelPart::NodesContainerType              NodesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    BaseMortarConvergenceCriteria()
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
    }

    ///Copy constructor 
    BaseMortarConvergenceCriteria( BaseMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    ~BaseMortarConvergenceCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        // We update the normals if necessary
        if (rModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            // Update normal of the conditions
            ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); 
        }
        
        // We recalculate the penalty parameter
        if (rModelPart.GetProcessInfo()[ADAPT_PENALTY] == true)
        {
            AALMAdaptPenaltyValueProcess aalm_adaptation_of_penalty = AALMAdaptPenaltyValueProcess( rModelPart );
            aalm_adaptation_of_penalty.Execute();
        }
        
        return true;
    }
    
    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */

    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        // Set to zero the weighted gap
        ResetWeightedGap(rModelPart);
        
        ConditionsArrayType& conditions_array = rModelPart.GetSubModelPart("Contact").Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        #pragma omp parallel for 
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            it_cond->AddExplicitContribution(rModelPart.GetProcessInfo());
        }
        
        return true;
    }
    
    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_ERROR << "WARNING:: YOUR ARE CALLING THE BASE MORTAR CRITERIA" << std::endl;
    }
    
    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
    }
    
    /**
     * This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override 
    {
    }

    /**
     * This function initializes the non linear iteration
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void InitializeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
    }
    
    /**
     * This function finalizes the non linear iteration
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
    }
    
    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    /**
     * This method resets the weighted gap in the nodes of the problem
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     */
    
    virtual void ResetWeightedGap(ModelPart& rModelPart)
    {       
        NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());

        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            it_node->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
        }
    }
    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class BaseMortarConvergenceCriteria 

///@name Explicit Specializations
///@{

}  // namespace Kratos 

#endif /* KRATOS_BASE_MORTAR_CRITERIA_H  defined */

