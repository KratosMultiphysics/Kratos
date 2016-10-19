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

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    
    typedef typename BaseType::ElementsArrayType   ElementsArrayType;      
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    
    typedef Point<3> PointType;
    
    
    ResidualBasedBossakDisplacementContactScheme(double rAlpham = 0.0)
        : ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>(rAlpham)
    {
        // For pure Newmark Scheme
        mAlpha.f = 0.0;
        mAlpha.m = rAlpham;

        mNewmark.beta= (1.0 + mAlpha.f - mAlpha.m) * (1.0 + mAlpha.f - mAlpha.m) * 0.25;
        mNewmark.gamma= 0.5 + mAlpha.f - mAlpha.m;

        // Allocate auxiliary memory
        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(NumThreads);
        mMatrix.D.resize(NumThreads);

        mVector.v.resize(NumThreads);
        mVector.a.resize(NumThreads);
        mVector.ap.resize(NumThreads);
    }
    
    /**
     * This function is designed to be called in the builder and solver to introduce
     * @param rCurrentCondition: The condition to compute
     * @param LHS_Contribution: The LHS matrix contribution
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */
    
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
            const int thread = OpenMPUtils::ThisThread();
                    
            // Basic operations for the condition considered
            (rCurrentCondition)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

            (rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);

            (rCurrentCondition)->CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentCondition)->CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

            AddDynamicsToLHS  (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

            AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);
        }
        
        KRATOS_CATCH("");
    }
    
    /**
     * This function is designed to calculate just the RHS contribution
     * @param rCurrentElemen: The element to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */
    
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
            const int thread = OpenMPUtils::ThisThread();
            
            // Basic operations for the condition considered
            (rCurrentCondition)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

            (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);

            (rCurrentCondition)->CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentCondition)->CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

            // Adding the dynamic contributions (static is already included)
            AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);
        }
        
        KRATOS_CATCH("");
    }
    
    /** 
     * Function that returns the list of Degrees of freedom to be assembled in the system for a Given condition
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
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
        (rCurrentCondition)->InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param rCurrentElement: The element to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(
      Element::Pointer rCurrentElement,
      ProcessInfo& CurrentProcessInfo
    )
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
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

    struct GeneralAlphaMethod
    {
        double f;  // Alpha Hilbert
        double m;  // Alpha Bosssak
    };

    struct NewmarkMethod
    {
        double beta;
        double gamma;

        // System constants
        double c0;
        double c1;
        double c2;
        double c3;
        double c4;
        double c5;
        double c6;
    };

    struct  GeneralMatrices
    {
        std::vector< Matrix > M;     // First derivative matrix  (usually mass matrix)
        std::vector< Matrix > D;     // Second derivative matrix (usually damping matrix)
    };

    struct GeneralVectors
    {
        std::vector< Vector > v;    // Velocity
        std::vector< Vector > a;    // Acceleration
        std::vector< Vector > ap;   // Previous acceleration
    };

    GeneralAlphaMethod  mAlpha;
    NewmarkMethod       mNewmark;

    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Updating first time Derivative
     * @param CurrentVelocity: The current velocity
     * @param DeltaDisplacement: The increment of displacement
     * @param PreviousVelocity: The previous velocity
     * @param PreviousAcceleration: The previous acceleration
     */

    inline void UpdateVelocity(
        array_1d<double, 3 > & CurrentVelocity,
        const array_1d<double, 3 > & DeltaDisplacement,
        const array_1d<double, 3 > & PreviousVelocity,
        const array_1d<double, 3 > & PreviousAcceleration
    )
    {
        noalias(CurrentVelocity) =  (mNewmark.c1 * DeltaDisplacement - mNewmark.c4 * PreviousVelocity
                                     - mNewmark.c5 * PreviousAcceleration);
    }

    /**
     * Updating second time Derivative
     * @param CurrentVelocity: The current velocity
     * @param DeltaDisplacement: The increment of displacement
     * @param PreviousVelocity: The previous velocity
     * @param PreviousAcceleration: The previous acceleration
     */

    inline void UpdateAcceleration(
        array_1d<double, 3 > & CurrentAcceleration,
        const array_1d<double, 3 > & DeltaDisplacement,
        const array_1d<double, 3 > & PreviousVelocity,
        const array_1d<double, 3 > & PreviousAcceleration
    )
    {
        noalias(CurrentAcceleration) =  (mNewmark.c0 * DeltaDisplacement - mNewmark.c2 * PreviousVelocity
                                         -  mNewmark.c3 * PreviousAcceleration);
    }

    /**
     * It adds the dynamic LHS contribution of the elements: M*c0 + D*c1 + K
     * @param LHS_Contribution: The dynamic contribution for the LHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param CurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        )
    {

        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += M * (1.0 - mAlpha.m) * mNewmark.c0;

            // std::cout<<" Mass Matrix "<<M<<" coeficient "<<(1-mAlpha.m)*mNewmark.c0<<std::endl;
        }

        // Adding  damping contribution
        if (D.size1() != 0) // if D matrix declared
        {
            noalias(LHS_Contribution) += D * (1.0 - mAlpha.f) * mNewmark.c1;

        }
    }

    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param CurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m);

            rCurrentElement->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread];

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
            //KRATOS_WATCH( prod(M, mVector.a[thread] ) )

        }

        // Adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
        }
    }

    /**
     * It adds the dynamic RHS contribution of the condition: b - M*a - D*v
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param CurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m);

            rCurrentCondition->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread];

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        // Adding damping contribution
        // Damping contribution
        if (D.size1() != 0)
        {
            rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v [thread]);
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
