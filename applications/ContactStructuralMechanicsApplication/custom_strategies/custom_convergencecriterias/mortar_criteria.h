// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
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
template<class TSparseSpace, class TDenseSpace>
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
        // First we check if we are in the frictional case
        const bool FrictionalCase = rModelPart.Is(SLIP);
            
        if (mInitialPreviousState == false)
        {            
            // We iterate over the nodes
            NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
            
            auto numNodes = pNode.end() - pNode.begin();
            
            #pragma omp parallel for
            for(unsigned int i = 0; i < numNodes; i++) 
            {
                auto itNode = pNode.begin() + i;
 
                itNode->SetValue(AUXILIAR_ACTIVE, itNode->Is(ACTIVE));

                if (FrictionalCase == true)
                {
                    itNode->SetValue(AUXILIAR_SLIP, itNode->Is(SLIP));
                }

            }
            
            mInitialPreviousState = true;
            return false;
        }
        
        bool IsConverged = true;
        
        // We iterate again over the nodes
        NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            // NORMAL DIRECTION
            if (itNode->GetValue(AUXILIAR_ACTIVE) != itNode->Is(ACTIVE))
            {                            
                itNode->SetValue(AUXILIAR_ACTIVE, itNode->Is(ACTIVE));
                #pragma omp critical
                IsConverged = false;
            }
            
            if (FrictionalCase == true)
            {
                // TANGENT DIRECTION
                if (itNode->GetValue(AUXILIAR_SLIP) != itNode->Is(SLIP))
                {                            
                    itNode->SetValue(AUXILIAR_SLIP, itNode->Is(SLIP));
                    #pragma omp critical
                    IsConverged = false;
                }
            }
        }
        
        if (this->GetEchoLevel() > 0)
        {
            if (IsConverged == true)
            {
                if (FrictionalCase == true)
                {
                    std::cout << "Convergence is achieved in ACTIVE/INACTIVE mortar nodes check" << std::endl;
                }
                else
                {
                    std::cout << "Convergence is achieved in ACTIVE/INACTIVE and STICK/SLIP mortar nodes check" << std::endl;
                }
            }
            else
            {
                if (FrictionalCase == true)
                {
                    std::cout << "Convergence is not achieved in ACTIVE/INACTIVE and STICK/SLIP mortar nodes check. RECALCULATING...." << std::endl;
                }
                else
                {
                    std::cout << "Convergence is not achieved in ACTIVE/INACTIVE mortar nodes check. RECALCULATING...." << std::endl;
                }
            }
        }
        
        return IsConverged;
    }
    
    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart)
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }
    
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
    ) 
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

    bool mInitialPreviousState;
    
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

///@name Explicit Specializations
///@{

}  // namespace Kratos 

#endif /* KRATOS_MORTAR_CRITERIA_H  defined */

