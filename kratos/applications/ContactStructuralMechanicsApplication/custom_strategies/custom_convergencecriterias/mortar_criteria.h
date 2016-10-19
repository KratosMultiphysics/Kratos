// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_MORTAR_CRITERIA_H)
#define  KRATOS_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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

/** @brief Custom convergence criteria for the mortar condition 
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class MortarConvergenceCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( MortarConvergenceCriteria );

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
    MortarConvergenceCriteria()
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mInitialPreviousState = false;
    }

    ///Copy constructor 
    MortarConvergenceCriteria( MortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
      ,mInitialPreviousState(rOther.mInitialPreviousState)
    {
    }

    /// Destructor
    virtual ~MortarConvergenceCriteria() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * Criterias that need to be called after getting the solution 
     * @param 
     * @return 
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        // TODO: Add the criteria of AUXILIAR_SLIP
        if (mInitialPreviousState == false)
        {
            NodesArrayType& pNode  = rModelPart.GetSubModelPart("Contact").Nodes();
            NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
            NodesArrayType::iterator it_node_end   = pNode.ptr_end();
            
            for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
            {
                if (node_it->Is(SLAVE))
                {
                    if (node_it->Is(ACTIVE))
                    {
                        node_it->GetValue(AUXILIAR_ACTIVE) = true;
                    }
                    else
                    {
                        node_it->GetValue(AUXILIAR_ACTIVE) = false;
                    }
                    
                    if (node_it->Is(SLIP))
                    {
                        node_it->GetValue(AUXILIAR_SLIP) = true;
                    }
                    else
                    {
                        node_it->GetValue(AUXILIAR_SLIP) = false;
                    }
                }
            }
            
            mInitialPreviousState = true;
            return false;
        }
        
        bool is_converged        = true;
        bool active_is_converged = true;
        bool slip_is_converged   = true;
        
        NodesArrayType& pNode  = rModelPart.GetSubModelPart("Contact").Nodes();
        NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
        NodesArrayType::iterator it_node_end   = pNode.ptr_end();
        
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            if (node_it->Is(SLAVE))
            {
                // NORMAL DIRECTION
                bool aux_bool_normal;
                if (node_it->Is(ACTIVE))
                {
                    aux_bool_normal = true;
                }
                else
                {
                    aux_bool_normal = false;
                }
                if (node_it->GetValue(AUXILIAR_ACTIVE) != aux_bool_normal)
                {                            
                    node_it->GetValue(AUXILIAR_ACTIVE) = aux_bool_normal;
                    active_is_converged = false;
                }
                
                // TANGENT DIRECTION
                bool aux_bool_tangent;
                if (node_it->Is(SLIP))
                {
                    aux_bool_tangent = true;
                }
                else
                {
                    aux_bool_tangent = false;
                }
                if (node_it->GetValue(AUXILIAR_SLIP) != aux_bool_tangent)
                {                            
                    node_it->GetValue(AUXILIAR_SLIP) = aux_bool_tangent;
                    slip_is_converged = false;
                }
                
                is_converged = active_is_converged && slip_is_converged;
            }
        }
        
        if (this->GetEchoLevel() > 0)
        {
            if (is_converged == true)
            {
                std::cout << "Convergence is achieved in ACTIVE/INACTIVE and STICK/SLIP mortar nodes check" << std::endl;
            }
            else if ((not active_is_converged && slip_is_converged ) == true)
            {
                std::cout << "Convergence is not achieved in ACTIVE/INACTIVE mortar nodes check. RECALCULATING...." << std::endl;
            }
            else if ((active_is_converged && not slip_is_converged ) == true)
            {
                std::cout << "Convergence is not achieved in STICK/SLIP mortar nodes check. RECALCULATING...." << std::endl;
            }
            else
            {
                std::cout << "Convergence is not achieved in ACTIVE/INACTIVE and STICK/SLIP mortar nodes check. RECALCULATING...." << std::endl;
            }
        }
        
        return is_converged;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This function 
     * @param 
     * @return 
     */ 
    void Initialize(ModelPart& rModelPart)
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This function 
     * @param 
     * @return 
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        mInitialPreviousState = false;        
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This function 
     * @param 
     * @return 
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) {}

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

    bool       mInitialPreviousState;
    
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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class ClassName 


}  // namespace Kratos 

#endif /* KRATOS_MORTAR_CRITERIA_H  defined */

