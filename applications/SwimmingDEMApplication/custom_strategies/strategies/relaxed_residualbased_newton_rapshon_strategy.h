//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

#if !defined(KRATOS_RELAXED_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_RELAXED_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "swimming_dem_application_variables.h"

//default builder and solver
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

/**
 * @class RelaxedResidualBasedNewtonRaphsonStrategy
 * @ingroup KratosCore
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm including a relaxation term
 * @author Joaquin Gonzalez-Usua
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class RelaxedResidualBasedNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(RelaxedResidualBasedNewtonRaphsonStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseStrategyType;

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef RelaxedResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

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

    // /**
    //  * @brief Default constructor
    //  */
    explicit RelaxedResidualBasedNewtonRaphsonStrategy() : BaseType()
    {
    }

    // /**
    //  * @brief Default constructor. (with parameters)
    //  * @param rModelPart The model part of the problem
    //  * @param ThisParameters The configuration parameters
    //  */
    explicit RelaxedResidualBasedNewtonRaphsonStrategy(ModelPart& rModelPart)
        : RelaxedResidualBasedNewtonRaphsonStrategy(rModelPart, RelaxedResidualBasedNewtonRaphsonStrategy::GetDefaultParameters())
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit RelaxedResidualBasedNewtonRaphsonStrategy(ModelPart& rModelPart, Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters),
            mSolutionStepIsInitialized(false),
            mInitializeWasPerformed(false),
            mKeepSystemConstantDuringIterations(false)
    {
    }

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
    explicit RelaxedResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
            mpScheme(pScheme),
            mpConvergenceCriteria(pNewConvergenceCriteria),
            mReformDofSetAtEachStep(ReformDofSetAtEachStep),
            mCalculateReactionsFlag(CalculateReactions),
            mSolutionStepIsInitialized(false),
            mMaxIterationNumber(MaxIterations),
            mInitializeWasPerformed(false),
            mKeepSystemConstantDuringIterations(false)
    {
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit RelaxedResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
        mpScheme(pScheme),
        mpBuilderAndSolver(pNewBuilderAndSolver),
        mpConvergenceCriteria(pNewConvergenceCriteria),
        mReformDofSetAtEachStep(ReformDofSetAtEachStep),
        mCalculateReactionsFlag(CalculateReactions),
        mSolutionStepIsInitialized(false),
        mMaxIterationNumber(MaxIterations),
        mInitializeWasPerformed(false),
        mKeepSystemConstantDuringIterations(false)
    {
        KRATOS_TRY
        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // Tells to the builder and solver if the reactions have to be Calculated or not
        p_builder_and_solver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        p_builder_and_solver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();
        KRATOS_CATCH("")
    }


    /**
     * Constructor with Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param Settings Settings used in the strategy
     */
    RelaxedResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters Settings)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, Settings),
            mpScheme(pScheme),
            mpConvergenceCriteria(pNewConvergenceCriteria),
            mSolutionStepIsInitialized(false),
            mInitializeWasPerformed(false),
            mKeepSystemConstantDuringIterations(false)
    {
    }

    /**
     * @brief Constructor specifying the builder and solver and using Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param Settings Settings used in the strategy
     */
    RelaxedResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters Settings)
        : BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, Settings),
            mpScheme(pScheme),
            mpBuilderAndSolver(pNewBuilderAndSolver),
            mpConvergenceCriteria(pNewConvergenceCriteria),
            mSolutionStepIsInitialized(false),
            mInitializeWasPerformed(false),
            mKeepSystemConstantDuringIterations(false)
    {
    }

    /**
    * @brief Destructor.
    * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
    */
    ~RelaxedResidualBasedNewtonRaphsonStrategy() override
    {
        // If the linear solver has not been deallocated, clean it before
        // deallocating mpA. This prevents a memory error with the ML
        // solver (which holds a reference to it).
        // NOTE: The linear solver is hold by the B&S
        auto p_builder_and_solver = this->GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->Clear();
        }

        // Deallocating system vectors to avoid errors in MPI. Clear calls
        // TrilinosSpace::Clear for the vectors, which preserves the Map of
        // current vectors, performing MPI calls in the process. Due to the
        // way Python garbage collection works, this may happen after
        // MPI_Finalize has already been called and is an error. Resetting
        // the pointers here prevents Clear from operating with the
        // (now deallocated) vectors.
        mpA.reset();
        mpDx.reset();
        mpb.reset();

        Clear();
    }
    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        BaseType::SetScheme(pScheme);
    };

    // /**
    //  * @brief Get method for the time scheme
    //  * @return mpScheme: The pointer to the time scheme considered
    //  */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

//*********************************************************************************
/**OPERATIONS ACCESSIBLE FROM THE INPUT: **/
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    /**
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    /**
     * @brief This method sets the flag mFullUpdateFlag
     * @param UseOldStiffnessInFirstIterationFlag The flag that tells if
     */
    void SetUseOldStiffnessInFirstIterationFlag(bool UseOldStiffnessInFirstIterationFlag)
    {
        mUseOldStiffnessInFirstIteration = UseOldStiffnessInFirstIterationFlag;
    }



    bool GetUseOldStiffnessInFirstIterationFlag()
    {
        return mUseOldStiffnessInFirstIteration;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    /**
     * @brief This method sets the flag mMaxIterationNumber
     * @param MaxIterationNumber This is the maximum number of on linear iterations
     */
    void SetMaxIterationNumber(unsigned int MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }

    /**
     * @brief This method gets the flag mMaxIterationNumber
     * @return mMaxIterationNumber: This is the maximum number of on linear iterations
     */
    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
    }
    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }
    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        BaseType::Clear();
    }
    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (!mSolutionStepIsInitialized) {
            // Pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
            ModelPart& r_model_part = BaseType::GetModelPart();

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
            {
                //setting up the list of the DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
                KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup Dofs Time: "
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                //shaping correctly the system
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(r_model_part);
                KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup System Time: "
                    << setup_system_time.ElapsedSeconds() << std::endl;

                //setting up the Vectors involved to the correct size
                BuiltinTimer system_matrix_resize_time;
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb,
                                                                 r_model_part);
                KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Matrix Resize Time: "
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;
            }

            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Construction Time: "
                << system_construction_time.ElapsedSeconds() << std::endl;

            TSystemMatrixType& rA  = *mpA;
            TSystemVectorType& rDx = *mpDx;
            TSystemVectorType& rb  = *mpb;

            // Initial operations ... things that are constant over the Solution Step
            p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);

            // Initial operations ... things that are constant over the Solution Step
            p_scheme->InitializeSolutionStep(r_model_part, rA, rDx, rb);

            // Initialisation of the convergence criteria
            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(rb);
                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            mpConvergenceCriteria->InitializeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                TSparseSpace::SetToZero(rb);

            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& r_model_part = BaseType::GetModelPart();

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        p_scheme->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        p_builder_and_solver->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        //Cleaning memory after the solution
        p_scheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            this->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;
        double alpha = r_model_part.GetProcessInfo()[RELAXATION_ALPHA];
        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            if (mUseOldStiffnessInFirstIteration){
                p_builder_and_solver->BuildAndSolveLinearizedOnPreviousIteration(p_scheme, r_model_part, rA, rDx, rb,BaseType::MoveMeshFlag());
            } else {
                p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
            }
        } else {
            TSparseSpace::SetToZero(rDx);  // Dx = 0.00;
            TSparseSpace::SetToZero(rb);
            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }
        // Debugging info
        EchoInfo(iteration_number);
        TSparseSpace::InplaceMult(rDx, alpha);
        // Updating the results stored in the database
        UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());
        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (is_converged) {
            if (mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
        }

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false &&
               iteration_number++ < mMaxIterationNumber)
        {
            //setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (GetKeepSystemConstantDuringIterations() == false)
                    {
                        //A = 0.00;
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                        TSparseSpace::InplaceMult(rDx, alpha);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
            }
            else
            {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            EchoInfo(iteration_number);

            // Updating the results stored in the database
            UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            residual_is_updated = false;

            if (is_converged == true)
            {
                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber) {
            MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("RelaxedResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << mMaxIterationNumber << " iterations" << std::endl;
        }

        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (residual_is_updated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(mb);
            //    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, mb);
        }

        //calculate reactions if required
        if (mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        return is_converged;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "relaxed_newton_raphson_strategy";
    }

    ///@}
    ///@name Operators

    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access

    ///@{
            /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    TSystemMatrixType &GetSystemMatrix() override
    {
        TSystemMatrixType &mA = *mpA;

        return mA;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    TSystemVectorType& GetSystemVector() override
    {
        TSystemVectorType& mb = *mpb;

        return mb;
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    TSystemVectorType& GetSolutionVector() override
    {
        TSystemVectorType& mDx = *mpDx;

        return mDx;
    }

    /**
     * @brief Set method for the flag mKeepSystemConstantDuringIterations
     * @param Value If we consider constant the system of equations during the iterations
     */
    void SetKeepSystemConstantDuringIterations(bool Value)
    {
        mKeepSystemConstantDuringIterations = Value;
    }

    /**
     * @brief Get method for the flag mKeepSystemConstantDuringIterations
     * @return True if we consider constant the system of equations during the iterations, false otherwise
     */
    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "RelaxedResidualBasedNewtonRaphsonStrategy";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

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

    typename TSchemeType::Pointer mpScheme = nullptr; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver = nullptr; /// The pointer to the builder and solver employed
    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria = nullptr; /// The pointer to the convergence criteria employed

    TSystemVectorPointerType mpDx; /// The increment in the solution
    TSystemVectorPointerType mpb; /// The RHS vector of the system of equations
    TSystemMatrixPointerType mpA; /// The LHS matrix of the system of equations

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    /**
     * @brief Flag telling if it is needed or not to compute the reactions
     * @details default = true
     */
    bool mCalculateReactionsFlag;

    /**
     * @brief Flag telling if a full update of the database will be performed at the first iteration
     * @details default = false
     */
    bool mUseOldStiffnessInFirstIteration = false;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    unsigned int mMaxIterationNumber; /// The maximum number of iterations, 30 by default

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mKeepSystemConstantDuringIterations; // Flag to allow keeping system matrix constant during iterations

    virtual void UpdateDatabase(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh) override
    {
        BaseType::UpdateDatabase(rA,rDx,rb,MoveMesh);
    }


    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber) override
    {
        BaseType::EchoInfo(IterationNumber);
    }

    virtual void MaxIterationsExceeded() override
    {
        KRATOS_INFO_IF("RelaxedResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
            << "ATTENTION: max iterations ( " << mMaxIterationNumber
            << " ) exceeded!" << std::endl;
    }

        /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
    }

    void WriteDofInfo(std::string FileName, const TSystemVectorType& rDX)
    {
        BaseType::WriteDofInfo(FileName, rDX);
    }
    /*
    * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
    */

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

    /*
     * Copy constructor.
     */

    RelaxedResidualBasedNewtonRaphsonStrategy(const RelaxedResidualBasedNewtonRaphsonStrategy &Other){};

    ///@}

}; /* Class RelaxedResidualBasedNewtonRaphsonStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_RELAXED_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY  defined */