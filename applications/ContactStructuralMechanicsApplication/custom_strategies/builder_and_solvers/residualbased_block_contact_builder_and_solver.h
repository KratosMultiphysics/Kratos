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

#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_CONTACT_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_BLOCK_CONTACT_BUILDER_AND_SOLVER

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

VERSION MODIFIED TO FIT THE NEEDS OF THE CONTACT PROBLEM

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}

 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockContactBuilderAndSolver
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockContactBuilderAndSolver);

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

    /** Constructor.
     */
    ResidualBasedBlockContactBuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
//         std::cout << "Using the standard builder and solver " << std::endl;
    }

    /** Destructor.
     */
    ~ResidualBasedBlockContactBuilderAndSolver() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * Function to perform the building and solving phase at the same time.
     * It is ideally the fastest and safer function to use when it is possible to solve just after building
     */    
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        this->Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "Before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "Unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        this->SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "System solve time: " << stop_solve - start_solve << std::endl;
        }

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "After the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "Unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
     * Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        this->BuildRHS(pScheme, rModelPart, b);
        this->SystemSolve(A, Dx, b);

        KRATOS_CATCH("")
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
    
    ///@}

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedBlockContactBuilderAndSolver */

///@}
///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_CONTACT_BUILDER_AND_SOLVER  defined */
