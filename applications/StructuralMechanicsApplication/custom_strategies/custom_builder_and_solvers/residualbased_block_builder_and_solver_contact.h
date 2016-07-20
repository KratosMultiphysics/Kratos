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

#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_CONTACT )
#define  KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_CONTACT

/* System includes */

#include "utilities/openmp_utils.h"


#include <unordered_set>
/* External includes */
#include "boost/smart_ptr.hpp"

#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "custom_utilities/contact_utilities.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"

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
    
/** \brief 
 * 
*/

template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolverContact
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverContact);


    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ResidualBasedBlockBuilderAndSolverContact(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
    }

    /// Destructor.
    virtual ~ResidualBasedBlockBuilderAndSolverContact()
    {
    }

    ///@}
    ///@name Operators
    ///@{
    
    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY;

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
        {
            norm_b = TSparseSpace::TwoNorm(b);
        }
        else
        {
            norm_b = 0.00;
        }

        if (norm_b != 0.00)
        {
            // Do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
        {
            TSparseSpace::SetToZero(Dx);
        }

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY;

        BuildRHS(pScheme, rModelPart, b);
        SystemSolve(A, Dx, b);

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY;
        
        BuildRHSNoDirichlet(pScheme,rModelPart,b);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k < ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            const std::size_t i = dof_iterator->EquationId();

            if (dof_iterator->IsFixed() == true)
            {
                b[i] = 0.0f;       
            }
        }

        KRATOS_CATCH("");

    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        )
    {
        TSparseSpace::SetToZero(b);

        // Refresh RHS to have the correct reactions
        BuildRHSNoDirichlet(pScheme, rModelPart, b);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: DOFS are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k < ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;

            if (dof_iterator->IsFixed())
            {
                const int i = (dof_iterator)->EquationId();
                (dof_iterator)->GetSolutionStepReactionValue() = -b[i];
            }
        }
    }
    
    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
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

    void AssembleRHS(
        TSystemVectorType& b,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
        )
    {
        unsigned int local_size = RHS_Contribution.size();


        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            // ASSEMBLING THE SYSTEM VECTOR
            double& b_value = b[i_global];
            const double& rhs_value = RHS_Contribution[i_local];

            #pragma omp atomic
            b_value += rhs_value;
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
    
    
    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY;

        // Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& ConditionsArray = rModelPart.Conditions();

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        // Assemble all elements
        // For (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)

        unsigned int nelements = static_cast<int>(pElements.size());
        #pragma omp parallel for firstprivate(nelements, RHS_Contribution, EquationId)
        for(unsigned int i = 0; i < nelements; i++)
        {
            typename ElementsArrayType::iterator it = pElements.begin() + i;
            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool element_is_active = true;
            if( (it)->IsDefined(ACTIVE) )
            {
                element_is_active = (it)->Is(ACTIVE);
            }

            if(element_is_active)
            {
                // Calculate elemental Right Hand Side Contribution
                pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                // Assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId);
            }
        }
        
        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart ); 

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // Assemble all conditions
        // For (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        unsigned int nconditions = static_cast<int>(ConditionsArray.size());
        #pragma omp parallel for firstprivate(nconditions, RHS_Contribution, EquationId)
        for (unsigned int i = 0; i < nconditions; i++)
        {
            auto it = ConditionsArray.begin() + i;
            // Detect if the element is active or not. If the user did not make any choice the element
            // It is active by default
            bool condition_is_active = true;
            if( (it)->IsDefined(ACTIVE) )
            {
                condition_is_active = (it)->Is(ACTIVE);
            }

            if(condition_is_active)
            {
                // Calculate elemental contribution
                pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                // Assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId);
            }
        }

        KRATOS_CATCH("");

    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; /* Class ResidualBasedBlockBuilderAndSolverContact */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_CONTACT  defined */
