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
//       ,mPreviousState(rOther.mPreviousState)
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
        if (mInitialPreviousState == false)
        {
            NodesArrayType& pNode  = rModelPart.Nodes();
            NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
            NodesArrayType::iterator it_node_end   = pNode.ptr_end();
            
            for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
            {
                if (node_it->Is(INTERFACE))
                {
                    if (node_it->Is(ACTIVE))
                    {
                        node_it->GetValue(AUXILIAR_BOOLEAN) = true;
                    }
                    else
                    {
                        node_it->GetValue(AUXILIAR_BOOLEAN) = false;
                    }
//                     mPreviousState.push_back(aux_bool);   
                }
            }
            
            mInitialPreviousState = true;
            return false;
        }
        
        bool is_converged = true;
        
        NodesArrayType& pNode  = rModelPart.Nodes();
        NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
        NodesArrayType::iterator it_node_end   = pNode.ptr_end();
        
//         unsigned int aux_index = 0;
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            if (node_it->Is(INTERFACE))
            {
                bool aux_bool;
                if (node_it->Is(ACTIVE))
                {
                    aux_bool = true;
                }
                else
                {
                    aux_bool = false;
                }
                if (node_it->GetValue(AUXILIAR_BOOLEAN) != aux_bool)
                {                            
                    node_it->GetValue(AUXILIAR_BOOLEAN) = aux_bool;
                    is_converged = false;
                }
//                 aux_index += 1;
            }
        }
        
        if (this->GetEchoLevel() > 0)
        {
            if (is_converged == true)
            {
                std::cout << "Convergence is achieved in ACTIVE/INACTIVE mortar nodes check" << std::endl;
            }
            else
            {
                std::cout << "Convergence is not achieved in ACTIVE/INACTIVE mortar nodes check. RECALCULATING...." << std::endl;
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
//     std::vector<bool> mPreviousState;
    
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

