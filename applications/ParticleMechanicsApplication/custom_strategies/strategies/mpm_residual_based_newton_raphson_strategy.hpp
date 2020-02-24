//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


#if !defined(KRATOS_MPM_RESIDUAL_BASED_NEWTON_RAPHSON_STRATEGY )
#define      KRATOS_MPM_RESIDUAL_BASED_NEWTON_RAPHSON_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"

#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

// Application includes
#include "particle_mechanics_application_variables.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/// Short class definition.

/**
 * @class MPMResidualBasedNewtonRaphsonStrategy
 * @ingroup KratosParticle
 * @brief Newton Raphson strategy suited for MPM simulations
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is achieved) using a Newton Raphson algorithm
 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class MPMResidualBasedNewtonRaphsonStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(MPMResidualBasedNewtonRaphsonStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructors.
     */
    MPMResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
    }

    MPMResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY;

        mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        SetMaxIterationNumber(MaxIterations);

        mCalculateReactionsFlag = CalculateReactions;

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (mpLinearSolver)
                             );

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;

        mInitializeWasPerformed = false;

        mFinalizeSolutionStep = true;

        // Tells to the Builder And Solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        KRATOS_CATCH( "" );
    }

    MPMResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY;

        mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = pNewBuilderAndSolver;

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;

        mInitializeWasPerformed = false;

        mFinalizeSolutionStep = true;

        // Tells to the Builder And Solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        KRATOS_CATCH( "" );
    }

    /** Destructor.
     */
    virtual ~MPMResidualBasedNewtonRaphsonStrategy()
    {
    }

    /** Destructor.
     */

    // Set and Get Scheme ... containing Builder, Update and other

    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        mReformDofSetAtEachStep = flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    void SetMaxIterationNumber(unsigned int MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }

    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    void SetFinalizeSolutionStepFlag(bool FinalizeSolutionStepFlag = true)
    {
        mFinalizeSolutionStep = FinalizeSolutionStepFlag;
    }

    bool GetFinalizeSolutionStepFlag()
    {
        return mFinalizeSolutionStep;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     */
    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
        mpConvergenceCriteria->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY;
        // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        // If the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
            Initialize();

        // Initialize solution step
        if (mSolutionStepIsInitialized == false)
            InitializeSolutionStep();

        DofsArrayType& r_dof_set = GetBuilderAndSolver()->GetDofSet();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        GetScheme()->Predict(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

        // Move the mesh if needed
        if (this->MoveMeshFlag() == true) BaseType::MoveMesh();

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        // if the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
        {
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy",this->GetEchoLevel() >1) << "Initializing solving strategy" << std::endl;
            KRATOS_ERROR_IF(mInitializeWasPerformed == true) << "Initialization was already performed " << mInitializeWasPerformed << std::endl;

            // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy",this->GetEchoLevel() >1) << "Initializing scheme" << std::endl;
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            // Initialize The Elements - OPERATIONS TO BE DONE ONCE
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy",this->GetEchoLevel() >1) << "Initializing elements" << std::endl;
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            // Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy",this->GetEchoLevel() >1) << "Initializing conditions" << std::endl;
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            // Initialisation of the convergence criteria
            typename TConvergenceCriteriaType::Pointer p_convergence_criteria = mpConvergenceCriteria;
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy",this->GetEchoLevel() >1) << "Initializing convergence criteria"<<std::endl;
            if (p_convergence_criteria->IsInitialized() == false)
                p_convergence_criteria->Initialize(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        // Set up the system, operation performed just once unless it is required to reform the dof set at each iteration
        if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
        {
            // Setting up the list of the DOFs to be solved
            p_builder_and_solver->SetUpDofSet(p_scheme, BaseType::GetModelPart());

            // Shaping correctly the system
            p_builder_and_solver->SetUpSystem(BaseType::GetModelPart());

        }

        // Prints informations about the current time
        if (this->GetEchoLevel() == 2 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb = *mpb;
        DofsArrayType& r_dof_set = p_builder_and_solver->GetDofSet();

        // Initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;

        p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);
        is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

        KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "PreCriteria:"
            << "\tIs_converged: " << is_converged  << "\tmRebuildLevel: " << BaseType::mRebuildLevel
            << "\tmStiffnessMatrixIsBuilt: " << BaseType::mStiffnessMatrixIsBuilt << std::endl;

        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "SetToZero the matrix and vectors of the system"<<std::endl;

            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build and Solve"<<std::endl;

            p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }
        else
        {
            TSparseSpace::SetToZero(rDx); // rDx=0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "BuildRHSAndSolve"<<std::endl;
        }


        if (this->GetEchoLevel() == 3) // If it is needed to print the debug info
        {
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "solution obtained = " << rDx << std::endl;
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "RHS  = " << rb << std::endl;
        }
        else if (this->GetEchoLevel() == 4) // Print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rb);
        }

        // Update results
        r_dof_set = p_builder_and_solver->GetDofSet();
        p_scheme->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

        if (is_converged == true)
        {
            // Initialisation of the convergence criteria
            r_dof_set = p_builder_and_solver->GetDofSet();
            mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(rb);
                p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        }

        KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Starting Nonlinear iteration"<<std::endl;

        // Iteration Loop
        while (is_converged == false &&
                iteration_number++<mMaxIterationNumber)
        {
            // Setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

            is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

            // Call the linear system solver to find the correction rDx. It is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false )
                {
                    KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Iteration Number: " << iteration_number <<std::endl;

                    if( GetKeepSystemConstantDuringIterations() == false)
                    {
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build and Solve"<<std::endl;
                        p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build RHS and Solve" <<std::endl;
                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build RHS and Solve" <<std::endl;
                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                }
            }
            else
            {
                KRATOS_WARNING("MPMNewtonRaphsonStrategy") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Updating the results stored in the database
            r_dof_set = p_builder_and_solver->GetDofSet();

            p_scheme->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

            p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

            // Move the mesh if needed
            if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

            // If converged
            if (is_converged == true)
            {
                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);

                }

                is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
            }
        }


        // Plot a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            if (this->GetEchoLevel() > 1)
                MaxIterationsExceeded();
        }

        return true;
    }

    /**
     * @brief This should be considered as a "post solution" convergence check which is useful for coupled analysis - the convergence criteria used is the one used inside the "solve" step
     */
    bool IsConverged() override
    {
        KRATOS_TRY;

        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb = *mpb;

        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            GetBuilderAndSolver()->BuildRHS(GetScheme(), BaseType::GetModelPart(), rb);
        }

        DofsArrayType& r_dof_set = GetBuilderAndSolver()->GetDofSet();

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        KRATOS_CATCH( "" );

    }

    /**
     * @brief This operations should be called before printing the results when non trivial results
     * (e.g. stresses)
     * Need to be calculated given the solution of the step
     * @details This operations should be called only when needed, before printing as it can involve a non
     * negligible cost
     */
    void CalculateOutputData() override
    {
        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb = *mpb;

        DofsArrayType& r_dof_set = GetBuilderAndSolver()->GetDofSet();
        GetScheme()->CalculateOutputData(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        SparseSpaceType::Clear(mpA);
        TSystemMatrixType& rA = *mpA;
        SparseSpaceType::Resize(rA, 0, 0);

        SparseSpaceType::Clear(mpDx);
        TSystemVectorType& rDx = *mpDx;
        SparseSpaceType::Resize(rDx, 0);

        SparseSpaceType::Clear(mpb);
        TSystemVectorType& rb = *mpb;
        SparseSpaceType::Resize(rb, 0);

        // Setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        KRATOS_CATCH( "" );
    }

    /*@} */
    /**@name Operators
     */
    /*@{ */

    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& rA = *mpA;
        return rA;
    }

    /**
     * @brief Set method for the flag mKeepSystemConstantDuringIterations
     * @param Value If we consider constant the system of equations during the iterations
     */
    void SetKeepSystemConstantDuringIterations(bool value)
    {
        mKeepSystemConstantDuringIterations = value;
    }

    /**
     * @brief Get method for the flag mKeepSystemConstantDuringIterations
     * @return True if we consider constant the system of equations during the iterations, false otherwise
     */
    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }

    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */

    /*@} */

private:

    /**@name Protected static Member Variables */
    /*@{ */

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */

    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

    /*@} */

protected:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */

    typename TLinearSolver::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employed
    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria; /// The pointer to the convergence criteria employed

    TSystemVectorPointerType mpDx; /// The incremement in the solution
    TSystemVectorPointerType mpb;  /// The RHS vector of the system of equations
    TSystemMatrixPointerType mpA;  /// The LHS matrix of the system of equations

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reform at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    // Flag telling if it is needed or not to compute the reactions
    bool mCalculateReactionsFlag;

    // Flag to set as initialized the solution step
    bool mSolutionStepIsInitialized;

    // The maximum number of iterations, 30 by default
    unsigned int mMaxIterationNumber;

    // Flag to set as initialized the strategy
    bool mInitializeWasPerformed;

    // Flag to allow keeping system matrix constant during iterations
    bool mKeepSystemConstantDuringIterations;

    // Flag to allow to not finalize the solution step, so the historical variables are not updated
    bool mFinalizeSolutionStep;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        // Initialize solution step
        if (mSolutionStepIsInitialized == false)
        {
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
            typename TSchemeType::Pointer p_scheme = GetScheme();

            // Setting up the Vectors involved to the correct size
            p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb, BaseType::GetModelPart());

            TSystemMatrixType& rA = *mpA;
            TSystemVectorType& rDx = *mpDx;
            TSystemVectorType& rb = *mpb;

            // Initial operations ... things that are constant over the Solution Step
            p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);

            // Initial operations ... things that are constant over the Solution Step
            p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);

            mSolutionStepIsInitialized = true;
        }

        KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Initialize Solution Step in strategy finished" <<std::endl;

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        typename TSchemeType::Pointer p_scheme = GetScheme();

        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb = *mpb;

        if (mCalculateReactionsFlag == true)
        {
            p_builder_and_solver->CalculateReactions(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }

        // Calling r_dof_set
        DofsArrayType& r_dof_set = p_builder_and_solver->GetDofSet();

        /* Finalization of the solution step,
         * operations to be done after achieving convergence, for example the
         * Final Residual Vector (rb) has to be saved in there to avoid error accumulation
         */
        if( mFinalizeSolutionStep )
        {
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Calling FinalizeSolutionStep" <<std::endl;

            p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);
            p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);
            mpConvergenceCriteria->FinalizeSolutionStep(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        }

        // Cleaning memory after the solution
        p_scheme->Clean();

        // Reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            SparseSpaceType::Clear(mpA);
            SparseSpaceType::Clear(mpDx);
            SparseSpaceType::Clear(mpb);

            this->Clear();
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief This method prints information after reach the max number of iterations
     */
    void MaxIterationsExceeded()
    {
        KRATOS_WARNING("MPMNewtonRaphsonStrategy") << "***************************************************" << std::endl;
        KRATOS_WARNING("MPMNewtonRaphsonStrategy") << "******* ATTENTION: max iterations exceeded ********" << std::endl;
        KRATOS_WARNING("MPMNewtonRaphsonStrategy") << "***************************************************" << std::endl;

    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY;

        BaseType::Check();
        GetBuilderAndSolver()->Check(BaseType::GetModelPart());
        GetScheme()->Check(BaseType::GetModelPart());
        mpConvergenceCriteria->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH( "" );
    }


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /** Copy constructor.
     */
    MPMResidualBasedNewtonRaphsonStrategy(const MPMResidualBasedNewtonRaphsonStrategy& Other)
    {
    };


    /*@} */

}; /* Class MPMResidualBasedNewtonRaphsonStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_MPM_RESIDUAL_BASED_NEWTON_RAPHSON_STRATEGY  defined */

