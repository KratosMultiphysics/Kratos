// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Mohamed Khalil
//                   Vicente Mataix Ferrándiz
//

#if !defined(RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_CONTACT_SCHEME)
#define RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_CONTACT_SCHEME

/* System Includes */

/* External Includes */

/* Project includes */
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/point.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/contact_utilities.h"
#include "includes/condition.h"

namespace Kratos {

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
    
/** \brief  Short class definition.
This class provides the implementation of the condition contribution methods of the solution strategy
for a contact problem employing the standard lagrange multipliers constraining formulation 
*/

template<class TSparseSpace, class TDenseSpace >
class ResidualBasedIncrementalUpdateStaticContactScheme :
    public ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace >
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalUpdateStaticContactScheme );

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    
    typedef typename BaseType::ElementsArrayType   ElementsArrayType;      
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    
    typedef Point<3> PointType;
    
    
    ResidualBasedIncrementalUpdateStaticContactScheme()
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>( )
    {
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY;
        
        bool condition_is_active = true;
        if( (rCurrentCondition)->IsDefined(ACTIVE) == true)
        {
            condition_is_active = (rCurrentCondition)->Is(ACTIVE);
        }
        
        if (condition_is_active == true)
        {
            ( rCurrentCondition )->InitializeNonLinearIteration( CurrentProcessInfo );
            BaseType::Condition_CalculateSystemContributions( rCurrentCondition, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );
        }
        
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY;

        bool condition_is_active = true;
        if( (rCurrentCondition)->IsDefined(ACTIVE) == true)
        {
            condition_is_active = (rCurrentCondition)->Is(ACTIVE);
        }
        
        if (condition_is_active == true)
        {
            ( rCurrentCondition )->InitializeNonLinearIteration( CurrentProcessInfo );
            BaseType::Condition_Calculate_RHS_Contribution( rCurrentCondition, RHS_Contribution, EquationId, CurrentProcessInfo );
        }
        
        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /** 
     * Function that returns the list of Degrees of freedom to be assembled in the system for a Given Element
     */
    
    void GetConditionDofList(
        Condition::Pointer rCurrentCondition,
        Element::DofsVectorType& ConditionDofList,
        ProcessInfo& CurrentProcessInfo)
    {
        bool condition_is_active = true;
        if( (rCurrentCondition)->IsDefined(ACTIVE) == true)
        {
            condition_is_active = (rCurrentCondition)->Is(ACTIVE);
        }
        
        if (condition_is_active == true)
        {
            rCurrentCondition->GetDofList(ConditionDofList, CurrentProcessInfo);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    virtual void InitializeConditions(
        ModelPart& rModelPart)
    {
        KRATOS_TRY

        if( this->mElementsAreInitialized==false )
            KRATOS_THROW_ERROR(std::logic_error, "Before initilizing Conditions, initialize Elements FIRST","")

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            typename ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin() + ConditionPartition[k];
            typename ConditionsArrayType::iterator CondEnd = rModelPart.Conditions().begin() + ConditionPartition[k + 1];

            for (typename ConditionsArrayType::iterator itCond = CondBegin; itCond != CondEnd; itCond++)
            {
                bool condition_is_active = true;
                if( (itCond)->IsDefined(ACTIVE) == true)
                {
                    condition_is_active = (itCond)->Is(ACTIVE);
                }
                
                if ( condition_is_active == true )
                {
                    itCond->Initialize(); //function to initialize the condition
                }

            }

        }

        this->mConditionsAreInitialized = true;
        KRATOS_CATCH("")
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    virtual void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        //initialize solution step for all of the elements
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for ( typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            (it) -> InitializeSolutionStep(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for ( typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it )
        {
            bool condition_is_active = true;
            if( (it)->IsDefined(ACTIVE) == true)
            {
                condition_is_active = (it)->Is(ACTIVE);
            }
            
            if ( condition_is_active == true )
            {
                (it) -> InitializeSolutionStep(CurrentProcessInfo);
            }

        }
        KRATOS_CATCH("")
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void InitializeNonLinIteration(
      ModelPart& rModelPart,
      TSystemMatrixType& A,
      TSystemVectorType& Dx,
      TSystemVectorType& b
    )
    {
        KRATOS_TRY;
        
        // It resets the weighted gap and slip 
        ContactUtilities::ResetWeightedGapSlip(rModelPart);
        
        // Initializes the non-linear iteration for all the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            typename ElementsArrayType::iterator ElementsEnd   = rElements.begin() + ElementPartition[k + 1];

            for (typename ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                itElem->InitializeNonLinearIteration(CurrentProcessInfo);
            }
        }

        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart ); 
        
        // Initializes the non-linear iteration for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();
        
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            typename ConditionsArrayType::iterator ConditionsEnd   = rConditions.begin() + ConditionPartition[k + 1];

            for (typename ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                bool condition_is_active = true;
                if( (itCond)->IsDefined(ACTIVE) == true)
                {
                    condition_is_active = (itCond)->Is(ACTIVE);
                }
                
                if ( condition_is_active == true )
                {
                    itCond->InitializeNonLinearIteration(CurrentProcessInfo);
                }
            }
        }
        
        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * It initializes a non-linear iteration (for an individual condition)
     * @param rCurrentConditiont: The condition to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(
      Condition::Pointer rCurrentCondition,
      ProcessInfo& CurrentProcessInfo
    )
    {
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param rCurrentConditiont: The element to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(
      Element::Pointer rCurrentElement,
      ProcessInfo& CurrentProcessInfo
    )
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY;
        
        // Finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            typename ElementsArrayType::iterator ElementsEnd   = rElements.begin() + ElementPartition[k + 1];

            for (typename ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                itElem->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }

        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart ); 
        
        // Finalizes solution step for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            typename ConditionsArrayType::iterator ConditionsEnd   = rConditions.begin() + ConditionPartition[k + 1];

            for (typename ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                bool condition_is_active = true;
                if( (itCond)->IsDefined(ACTIVE)== true)
                {
                    condition_is_active = (itCond)->Is(ACTIVE);
                }
                
                if ( condition_is_active == true )
                {
                    itCond->FinalizeSolutionStep(CurrentProcessInfo);
                }
            }
        }
        
        ContactUtilities::ReComputeActiveInactive( rModelPart );  
        
        KRATOS_CATCH("");
    }
    
    /**
     * Function to be called when it is needed to finalize an iteration. It is designed to be called at the end of each non linear iteration
     */
    
    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        )
    {
        KRATOS_TRY;
        
        // Finalizes non linear iteration for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            typename ElementsArrayType::iterator ElementsEnd   = rElements.begin() + ElementPartition[k + 1];

            for (typename ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                itElem->FinalizeNonLinearIteration(CurrentProcessInfo);
            }
        }

        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart ); 

        // Finalizes non linear iteration for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            typename ConditionsArrayType::iterator ConditionsEnd   = rConditions.begin() + ConditionPartition[k + 1];

            for (typename ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                bool condition_is_active = true;
                if( (itCond)->IsDefined(ACTIVE)== true)
                {
                    condition_is_active = (itCond)->Is(ACTIVE);
                }
                
                if ( condition_is_active == true )
                {
                    itCond->FinalizeNonLinearIteration(CurrentProcessInfo);
                }
            }
        }
                
        ContactUtilities::ReComputeActiveInactive(rModelPart); 
        
        KRATOS_CATCH("");
    }
};

}  // namespace Kratos

#endif /* RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_CONTACT_SCHEME */
