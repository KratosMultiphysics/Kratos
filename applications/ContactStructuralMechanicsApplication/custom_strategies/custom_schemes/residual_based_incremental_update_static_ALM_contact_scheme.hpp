// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    icente Mataix Ferrándiz
//

#if !defined(RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_ALM_CONTACT_SCHEME)
#define RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_ALM_CONTACT_SCHEME

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
class ResidualBasedIncrementalUpdateStaticALMContactScheme :
    public ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace >
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalUpdateStaticALMContactScheme );

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
    
    
    ResidualBasedIncrementalUpdateStaticALMContactScheme()
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>( )
    {
    }
    
    /**
     * It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;
        
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        // Initializes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        const int nelem = static_cast<int>(rModelPart.Elements().size());
        typename ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin();

        #pragma omp parallel for
        for(int i = 0;  i < nelem; i++)
        {
            typename ElementsArrayType::iterator itElem = ElemBegin + i;
            
            itElem->InitializeSolutionStep(CurrentProcessInfo);
        }
        
        // Initializes solution step for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        const int ncond = static_cast<int>(rModelPart.Conditions().size());
        typename ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin();

        #pragma omp parallel for
        for(int i = 0;  i < ncond; i++)
        {
            typename ConditionsArrayType::iterator itCond = CondBegin + i;
            
            bool condition_is_active = true;
            if( (itCond)->IsDefined(ACTIVE)== true)
            {
                condition_is_active = (itCond)->Is(ACTIVE);
            }
            
            if ( condition_is_active == true )
            {
                itCond->InitializeSolutionStep(CurrentProcessInfo);
            }
        }
        
        KRATOS_CATCH("");
    }

    /**
     * It initializes a non-linear iteration (for the element)
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */
    
    void InitializeNonLinIteration(
      ModelPart& rModelPart,
      TSystemMatrixType& A,
      TSystemVectorType& Dx,
      TSystemVectorType& b
      ) override
    {
        KRATOS_TRY;
        
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        
        // Initializes the non-linear iteration for all the elements
        ElementsArrayType& rElements = rModelPart.Elements();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        const int nelem = static_cast<int>(rModelPart.Elements().size());
        typename ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin();

        #pragma omp parallel for
        for(int i = 0;  i < nelem; i++)
        {
            typename ElementsArrayType::iterator itElem = ElemBegin + i;
            
            itElem->InitializeNonLinearIteration(CurrentProcessInfo);
        }

        if (rModelPart.Is(INTERACTION) == true)
        {
            // Update normal of the conditions
            ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); 
        }
        
        // Reset the weighted variables
        ContactUtilities::ResetWeightedALMFrictionlessValues( rModelPart ); // FIXME: Just for frictionless case
        
        // Initializes the non-linear iteration for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();
        
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        const int ncond = static_cast<int>(rModelPart.Conditions().size());
        typename ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin();

        #pragma omp parallel for
        for(int i = 0;  i < ncond; i++)
        {
            typename ConditionsArrayType::iterator itCond = CondBegin + i;
            
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
        
        KRATOS_CATCH("");
    }
    
    /**
     * Function called once at the end of a solution step, after convergence is reached if
     * an iterative process is needed
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;
        
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        
        // Finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        const int nelem = static_cast<int>(rModelPart.Elements().size());
        typename ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin();

        #pragma omp parallel for
        for(int i = 0;  i < nelem; i++)
        {
            typename ElementsArrayType::iterator itElem = ElemBegin + i;
            
            itElem->FinalizeSolutionStep(CurrentProcessInfo);
        }

        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); 
        
        // Finalizes solution step for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        const int ncond = static_cast<int>(rModelPart.Conditions().size());
        typename ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin();

        #pragma omp parallel for
        for(int i = 0;  i < ncond; i++)
        {
            typename ConditionsArrayType::iterator itCond = CondBegin + i;
            
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
        
//         // It resets the visited flag
//         ContactUtilities::ResetVisited(rModelPart);
        
//         // It recomputes the active/inactive pair
//         ContactUtilities::ReComputeActiveInactiveALMFrictionless( rModelPart );  // FIXME: Just for frictionless case
        
        KRATOS_CATCH("");
    }
    
    /**
     * Function to be called when it is needed to finalize an iteration. It is designed to be called at the end of each non linear iteration
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */
    
    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;
        
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        
        // Finalizes non linear iteration for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        const int nelem = static_cast<int>(rModelPart.Elements().size());
        typename ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin();

        #pragma omp parallel for
        for(int i = 0;  i < nelem; i++)
        {
            typename ElementsArrayType::iterator itElem = ElemBegin + i;
            
            itElem->FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        if (rModelPart.Is(INTERACTION) == true)
        {
            // Update normal of the conditions
            ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); 
        }

        // Finalizes non linear iteration for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        const int ncond = static_cast<int>(rModelPart.Conditions().size());
        typename ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin();

        #pragma omp parallel for
        for(int i = 0;  i < ncond; i++)
        {
            typename ConditionsArrayType::iterator itCond = CondBegin + i;
            
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
            
        // It resets the visited flag
        ContactUtilities::ResetVisited(rModelPart);
        
        // It recomputes the active/inactive pair
        ContactUtilities::ReComputeActiveInactiveALMFrictionless( rModelPart );  // FIXME: Just for frictionless case
        
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
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
    ///@{

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
}; /* Class ResidualBasedIncrementalUpdateStaticALMContactScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

}  // namespace Kratos

#endif /* RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_ALM_CONTACT_SCHEME */
