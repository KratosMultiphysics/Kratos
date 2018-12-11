//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_HEXAHEDRA_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_HEXAHEDRA_NEWTON_RAPHSON_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Extrapolation process
#include "custom_processes/integration_values_extrapolation_to_nodes_process.h"


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
 * @class HexahedraNewtonRaphsonStrategy
 * @ingroup KratosCore
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 * @author Alejandro Cornejo
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class HexahedraNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(HexahedraNewtonRaphsonStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TDataType TDataType;
    typedef TSparseSpace SparseSpaceType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    HexahedraNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY;

        BaseType::mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        BaseType::SetMaxIterationNumber(MaxIterations);
        BaseType::mCalculateReactionsFlag = CalculateReactions;

        BaseType::mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        BaseType::mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        BaseType::mpScheme = pScheme;

        // Saving the linear solver
        BaseType::mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        BaseType::mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(BaseType::mpLinearSolver));

        // Set flags to start correcty the calculations
        BaseType::mSolutionStepIsInitialized = false;
        BaseType::mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        BaseType::GetBuilderAndSolver()->SetCalculateReactionsFlag(BaseType::mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        BaseType::GetBuilderAndSolver()->SetReshapeMatrixFlag(BaseType::mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        BaseType::SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        BaseType::mpA = TSparseSpace::CreateEmptyMatrixPointer();
        BaseType::mpDx = TSparseSpace::CreateEmptyVectorPointer();
        BaseType::mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    HexahedraNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY

        BaseType::mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        BaseType::SetMaxIterationNumber(MaxIterations);
        BaseType::mCalculateReactionsFlag = CalculateReactions;

        BaseType::mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        BaseType::mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        BaseType::mpScheme = pScheme;

        // Saving the linear solver
        BaseType::mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        BaseType::mpBuilderAndSolver = pNewBuilderAndSolver;

        // Set flags to start correcty the calculations
        BaseType::mSolutionStepIsInitialized = false;
        BaseType::mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        BaseType::GetBuilderAndSolver()->SetCalculateReactionsFlag(BaseType::mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        BaseType::GetBuilderAndSolver()->SetReshapeMatrixFlag(BaseType::mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        BaseType::SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        BaseType::mpA = TSparseSpace::CreateEmptyMatrixPointer();
        BaseType::mpDx = TSparseSpace::CreateEmptyVectorPointer();
        BaseType::mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~HexahedraNewtonRaphsonStrategy() override
    {
        BaseType::Clear();
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();

        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // We extrapolate the stress to the nodes to perform the smoothing
        ModelPart& rModelPart = BaseType::GetModelPart();
        Parameters extrapolation_parameters(R"(
        {
        	"list_of_variables": ["STRESS_VECTOR","STRAIN_VECTOR"]
        })");
        auto extrapolation_process = IntegrationValuesExtrapolationToNodesProcess(rModelPart, extrapolation_parameters);
        extrapolation_process.Execute();

        is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        //function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);
            p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        } else {
            TSparseSpace::SetToZero(rDx); //Dx=0.00;
            TSparseSpace::SetToZero(rb);
            p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }

        // Debugging info
        BaseType::EchoInfo(iteration_number);

        // Updating the results stored in the database
        BaseType::UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        if (is_converged == true) {
            //initialisation of the convergence criteria
            BaseType::mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                TSparseSpace::SetToZero(rb);
                p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);
            }
            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);
        }

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false && iteration_number++ < BaseType::mMaxIterationNumber) {
            //setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);
            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);
            extrapolation_process.Execute();

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0) {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false) {
                    if (BaseType::GetKeepSystemConstantDuringIterations() == false) {
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);
                        p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                    } else {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);
                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                    }
                } else {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);
                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                }
            } else {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            BaseType::EchoInfo(iteration_number);

            // Updating the results stored in the database
            BaseType::UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

            residual_is_updated = false;

            if (is_converged == true) {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);
                    residual_is_updated = true;
                }

                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            BaseType::MaxIterationsExceeded();

        //calculate reactions if required
        if (BaseType::mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);

        return is_converged;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
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

  protected:
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

    /**
     * Copy constructor.
     */

    HexahedraNewtonRaphsonStrategy(const HexahedraNewtonRaphsonStrategy &Other){};

    ///@}

}; /* Class HexahedraNewtonRaphsonStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY  defined */
