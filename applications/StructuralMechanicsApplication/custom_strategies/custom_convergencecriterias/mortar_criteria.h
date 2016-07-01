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
      ,mPreviousState(rOther.mPreviousState)
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
            ConditionsArrayType& pCond  = rModelPart.Conditions();
            ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
            ConditionsArrayType::iterator it_end   = pCond.ptr_end();
            
            for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
            {
                if (cond_it->Is(ACTIVE))
                {
                    std::vector<contact_container> *& ContactContainers = cond_it->GetValue(CONTACT_CONTAINERS);
                    for (unsigned int i = 0; i < ContactContainers->size(); i++)
                    {
                        for (unsigned int j = 0; j < (*ContactContainers)[i].active_nodes_slave.size(); j++)
                        {
                            mPreviousState.push_back((*ContactContainers)[i].active_nodes_slave[j]);   
                        }
                    }
                }
            }
            
            mInitialPreviousState = true;
            return false;
        }
        
        bool is_converged = true;
        
        ConditionsArrayType& pCond  = rModelPart.Conditions();
        ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_end   = pCond.ptr_end();
        
        unsigned int aux_index = 0;
        for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
        {
            if (cond_it->Is(ACTIVE))
            {
                
                if (this->GetEchoLevel() > 1)
                {
                    std::cout << std::endl;
                    std::cout << "The current state of active/inactive nodes for condition " << cond_it->Id() << " is:" << std::endl;
                    for (unsigned int i = 0; i < mPreviousState.size(); i++)
                    {
                        std::cout << mPreviousState[i]<< " ";
                    }
                    std::cout << std::endl;
                }
                
                std::vector<contact_container> *& ContactContainers = cond_it->GetValue(CONTACT_CONTAINERS);
                for (unsigned int i = 0; i < ContactContainers->size(); i++)
                {
                    for (unsigned int j = 0; j < (*ContactContainers)[i].active_nodes_slave.size(); j++)
                    {
                        if (mPreviousState[aux_index] != (*ContactContainers)[i].active_nodes_slave[j])
                        {                            
                            mPreviousState[aux_index] = (*ContactContainers)[i].active_nodes_slave[j];
                            is_converged = false;
                        }
                        aux_index += 1;
                    }
                }
                
                if (this->GetEchoLevel() > 1)
                {
                    std::cout << "The new state of active/inactive nodes for condition " << cond_it->Id() << " is:" << std::endl;
                    for (unsigned int i = 0; i < mPreviousState.size(); i++)
                    {
                        std::cout << mPreviousState[i]<< " ";
                    }
                    std::cout << std::endl;
                }
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
    std::vector<bool> mPreviousState;
    
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

