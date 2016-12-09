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

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_CONTACT_SCHEME )
#define  KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_CONTACT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/point.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
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
class ResidualBasedBossakDisplacementContactScheme :
    public ResidualBasedBossakDisplacementScheme< TSparseSpace, TDenseSpace >
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBossakDisplacementContactScheme );

    typedef ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType                                   TDataType;

    typedef typename BaseType::DofsArrayType                           DofsArrayType;

    typedef typename Element::DofsVectorType                          DofsVectorType;

    typedef typename BaseType::TSystemMatrixType                   TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                   TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType           LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType           LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                             NodesArrayType;

    typedef ModelPart::ElementsContainerType                       ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                   ConditionsArrayType;

    typedef typename BaseType::Pointer                               BaseTypePointer;
    
    ResidualBasedBossakDisplacementContactScheme(double rAlpham = 0.0)
        : ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>(rAlpham)
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
    )
    {
        KRATOS_TRY;

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
                
        // Initialize solution step for all of the elements
        ElementsArrayType& pElements = rModelPart.Elements();
        for ( typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            it->InitializeSolutionStep(CurrentProcessInfo);
        }
        
        // Initialize solution step for all of the conditions
        ConditionsArrayType& pConditions = rModelPart.Conditions();
        for ( typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it )
        {
            bool condition_is_active = true;
            if( it->IsDefined(ACTIVE) == true)
            {
                condition_is_active = it->Is(ACTIVE);
            }
            
            if ( condition_is_active == true )
            {
                it->InitializeSolutionStep(CurrentProcessInfo);
            }

        }

        double DeltaTime = CurrentProcessInfo[DELTA_TIME];

        double beta = 0.25;
        if (CurrentProcessInfo.Has(NEWMARK_BETA))
        {
            beta = CurrentProcessInfo[NEWMARK_BETA];
        }
        double gamma = 0.5;
        if (CurrentProcessInfo.Has(NEWMARK_GAMMA))
        {
            gamma = CurrentProcessInfo[NEWMARK_GAMMA];
        }

        BaseType::CalculateNewmarkCoefficients(beta, gamma);

        if (DeltaTime < 1.0e-24)
        {
            KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Scheme DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part ";
        }

        // Initializing Newmark constants
        BaseType::mNewmark.c0 = ( 1.0 / (BaseType::mNewmark.beta * DeltaTime * DeltaTime) );
        BaseType::mNewmark.c1 = ( BaseType::mNewmark.gamma / (BaseType::mNewmark.beta * DeltaTime) );
        BaseType::mNewmark.c2 = ( 1.0 / (BaseType::mNewmark.beta * DeltaTime) );
        BaseType::mNewmark.c3 = ( 0.5 / (BaseType::mNewmark.beta) - 1.0 );
        BaseType::mNewmark.c4 = ( (BaseType::mNewmark.gamma / BaseType::mNewmark.beta) - 1.0  );
        BaseType::mNewmark.c5 = ( DeltaTime * 0.5 * ( ( BaseType::mNewmark.gamma / BaseType::mNewmark.beta ) - 2.0 ) );

        KRATOS_CATCH( "" );
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
    )
    {
        KRATOS_TRY;

        // It resets the weighted gap and slip 
        ContactUtilities::ResetVisited(rModelPart);
        
        // Initializes the non-linear iteration for all the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

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

        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); 
        
        // Reset the weighted variables
        ContactUtilities::ResetWeightedValues( rModelPart ); 
        
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
        )
    {
        KRATOS_TRY;
            
        // Finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

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
       
//         ContactUtilities::ReComputeActiveInactive( rModelPart ); 
        
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
        )
    {
        KRATOS_TRY;
        
        // Finalizes non linear iteration for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

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
                
        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); 
        
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
        
        ContactUtilities::ReComputeActiveInactive( rModelPart); 
        
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
}; /* Class ResidualBasedBossakDisplacementContactScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_CONTACT_SCHEME defined */
