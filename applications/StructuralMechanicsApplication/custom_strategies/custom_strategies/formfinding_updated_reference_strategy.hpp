// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Anna Rehr
//

#if !defined(KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY )
#define  KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY
// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "includes/gid_io.h"


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

/**
 * @class FormfindingUpdatedReferenceStrategy
 *
 * @ingroup StrucutralMechanicsApplication
 *
 * @brief inherited class from ResidualBasedNewtonRaphsonStrategy for formfinding
 *
 * @details additions for formfinding: update the reference configuration for each element, initialize the elements for formfinding,
 * adaption line search for formfinding, print formfinding output (different nonlinear iterations)
 *
 * @author Anna Rehr
 */

    template<class TSparseSpace,
    class TDenseSpace, // = DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
    class FormfindingUpdatedReferenceStrategy
        : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        ///@name Type Definitions
        ///@{
        typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

        // Counted pointer of ClassName
        KRATOS_CLASS_POINTER_DEFINITION(FormfindingUpdatedReferenceStrategy);

        typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
        typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
        typedef typename BaseType::TSchemeType TSchemeType;
        typedef GidIO<> IterationIOType;
        typedef IterationIOType::Pointer IterationIOPointerType;
        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        ///@}
        ///@name Life Cycle

        ///@{

        /**
        * Constructor.
        */

        FormfindingUpdatedReferenceStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            typename TLinearSolver::Pointer pNewLinearSolver,
            typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
            int MaxIterations = 30,
            bool CalculateReactions = false,
            bool ReformDofSetAtEachStep = false,
            bool MoveMeshFlag = false,
            bool PrintIterations = false,
            bool IncludeLineSearch = false
            )
            : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewLinearSolver,
                pNewConvergenceCriteria,
                MaxIterations,
                CalculateReactions,
                ReformDofSetAtEachStep,
                MoveMeshFlag),
                mPrintIterations(PrintIterations),
                mIncludeLineSearch(IncludeLineSearch)
        {
        }

        // constructor with Builder and Solver
        FormfindingUpdatedReferenceStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            typename TLinearSolver::Pointer pNewLinearSolver,
            typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
            typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
            int MaxIterations = 30,
            bool CalculateReactions = false,
            bool ReformDofSetAtEachStep = false,
            bool MoveMeshFlag = false,
            bool PrintIterations = false,
            bool IncludeLineSearch = false
            )
            : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                model_part,
                pScheme,
                pNewLinearSolver,
                pNewConvergenceCriteria,
                pNewBuilderAndSolver,
                MaxIterations,
                CalculateReactions,
                ReformDofSetAtEachStep,
                MoveMeshFlag),
            mPrintIterations(PrintIterations),
            mIncludeLineSearch(IncludeLineSearch)
        {
        }

        /**
        * Destructor.
        */

        ~FormfindingUpdatedReferenceStrategy() override
        {
        }

        /**
        * InitializeSolutionStep. 
        * In addition to the base class initialization, the elements are initialized for formfinding
        * in each step.
        */
        void InitializeSolutionStep() override
        {
            for (auto& element : BaseType::GetModelPart().Elements())
                element.Initialize();
            for (auto& condition : BaseType::GetModelPart().Conditions())
                condition.Initialize();

            BaseType::InitializeSolutionStep();
        }

          ///@}
    protected:

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{

        bool mPrintIterations;
        bool mIncludeLineSearch;
        //IterationIOPointerType mpIterationIO;


        ///@}
        ///@name Private Operators
        ///@{

        /**
        * Copy constructor.
        */
        FormfindingUpdatedReferenceStrategy(const FormfindingUpdatedReferenceStrategy& Other)
        {
        };

        ///@}
    }; /* Class FormfindingUpdatedReferenceStrategy */

} /* namespace Kratos. */

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY defined */
