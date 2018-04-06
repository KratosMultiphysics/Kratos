//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/openmp_utils.h"

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
 * @class ResidualBasedNewtonRaphsonStrategy
 * @ingroup KratosCore
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 * @author Riccardo Rossi
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ResidualBasedNewtonRaphsonStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedNewtonRaphsonStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    
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
    ResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
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
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(mpLinearSolver));

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

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
    ResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY

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

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

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
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~ResidualBasedNewtonRaphsonStrategy() override
    {
        Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
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

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
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
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY
        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
        {
            Initialize();
        }

        //initialize solution step
        if (mSolutionStepIsInitialized == false)
        {
            InitializeSolutionStep();
        }

        TSystemMatrixType &A = *mpA;
        TSystemVectorType &Dx = *mpDx;
        TSystemVectorType &b = *mpb;

        GetScheme()->Predict(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), A, Dx, b);

        //move the mesh if needed
        if (this->MoveMeshFlag() == true)
        {
            BaseType::MoveMesh();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     */    
    void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            //pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (pScheme->SchemeIsInitialized() == false)
            {
                pScheme->Initialize(BaseType::GetModelPart());
            }

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (pScheme->ElementsAreInitialized() == false)
            {
                pScheme->InitializeElements(BaseType::GetModelPart());
            }

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (pScheme->ConditionsAreInitialized() == false)
            {
                pScheme->InitializeConditions(BaseType::GetModelPart());
            }

            //initialisation of the convergence criteria
            if (mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
            {
                mpConvergenceCriteria->Initialize(BaseType::GetModelPart());
            }

            //pointers needed in the solution
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
            const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
            {
                //setting up the list of the DOFs to be solved
                double setup_dofs_begintime = OpenMPUtils::GetCurrentTime();
                pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());
                if (this->GetEchoLevel() > 0 && rank == 0)
                {
                    double setup_dofs_endtime = OpenMPUtils::GetCurrentTime();
                    KRATOS_INFO("setup_dofs_time") << "setup_dofs_time : " << setup_dofs_endtime - setup_dofs_begintime << std::endl;
                }
                //shaping correctly the system
                double setup_system_begin = OpenMPUtils::GetCurrentTime();
                pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
                if (this->GetEchoLevel() > 0 && rank == 0)
                {
                    double setup_system_end = OpenMPUtils::GetCurrentTime();
                    KRATOS_INFO("setup_system_time") << rank << ": setup_system_time : " << setup_system_end - setup_system_begin << std::endl;
                }
            }

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief The problem of interest is solved.
     * @details This function calls sequentially: Initialize(), InitializeSolutionStep(), Predict(),
     * SolveSolutionStep() and FinalizeSolutionStep().
     * All those functions can otherwise be called separately.
     */    
    double Solve() override
    {
        Initialize();
        InitializeSolutionStep();
        Predict();
        SolveSolutionStep();
        FinalizeSolutionStep();
        return 0.00;
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpA != nullptr)
            SparseSpaceType::Clear(mpA);
        if (mpDx != nullptr)
            SparseSpaceType::Clear(mpDx);
        if (mpb != nullptr)
            SparseSpaceType::Clear(mpb);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        KRATOS_CATCH("");
    }

    /**
     * @brief This should be considered as a "post solution" convergence check which is useful for coupled analysis - the convergence criteria used is the one used inside the "solve" step
     */    
    bool IsConverged() override
    {
        KRATOS_TRY;

        TSystemMatrixType &A = *mpA;
        TSystemVectorType &Dx = *mpDx;
        TSystemVectorType &b = *mpb;

        if (mpConvergenceCriteria->mActualizeRHSIsNeeded == true)
        {
            GetBuilderAndSolver()->BuildRHS(GetScheme(), BaseType::GetModelPart(), b);
        }

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), A, Dx, b);

        KRATOS_CATCH("");
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
        TSystemMatrixType &A = *mpA;
        TSystemVectorType &Dx = *mpDx;
        TSystemVectorType &b = *mpb;

        GetScheme()->CalculateOutputData(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), A, Dx, b);
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        //pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        //set up the system, operation performed just once unless it is required
        //to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
            mReformDofSetAtEachStep == true)
        {
            //setting up the list of the DOFs to be solved
            double setup_dofs_begintime = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double setup_dofs_endtime = OpenMPUtils::GetCurrentTime();
                KRATOS_INFO("setup_dofs_time") << "setup_dofs_time : " << setup_dofs_endtime - setup_dofs_begintime << std::endl;
            }

            //shaping correctly the system
            double setup_system_begin = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double setup_system_end = OpenMPUtils::GetCurrentTime();
                KRATOS_INFO("setup_system_time") << rank << ": setup_system_time : " << setup_system_end - setup_system_begin << std::endl;
            }
        }

        //prints informations about the current time
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            KRATOS_INFO("CurrentTime") << "\nCurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        if (mSolutionStepIsInitialized == false)
        {
            /*typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
            typename TSchemeType::Pointer pScheme = GetScheme();
            int rank = BaseType::GetModelPart().GetCommunicator().MyPID(); */

            //setting up the Vectors involved to the correct size
            double system_matrix_resize_begin = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA, mpDx, mpb, BaseType::GetModelPart().Elements(), BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double system_matrix_resize_end = OpenMPUtils::GetCurrentTime();
                KRATOS_INFO("system_matrix_resize_time") << rank << ": system_matrix_resize_time : " << system_matrix_resize_end - system_matrix_resize_begin << std::endl;
            }
            TSystemMatrixType &A = *mpA;
            TSystemVectorType &Dx = *mpDx;
            TSystemVectorType &b = *mpb;

            //initial operations ... things that are constant over the Solution Step
            pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(), A, Dx, b);

            //initial operations ... things that are constant over the Solution Step
            pScheme->InitializeSolutionStep(BaseType::GetModelPart(), A, Dx, b);
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

        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        TSystemMatrixType &A = *mpA;
        TSystemVectorType &Dx = *mpDx;
        TSystemVectorType &b = *mpb;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), A, Dx, b);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), A, Dx, b);

        //Cleaning memory after the solution
        pScheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            SparseSpaceType::Clear(mpA);
            SparseSpaceType::Clear(mpDx);
            SparseSpaceType::Clear(mpb);

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
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        TSystemMatrixType &A = *mpA;
        TSystemVectorType &Dx = *mpDx;
        TSystemVectorType &b = *mpb;

        //initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        //			BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);
        bool is_converged = false;
        bool ResidualIsUpdated = false;
        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
        is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

        //function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(A);
            TSparseSpace::SetToZero(Dx);
            TSparseSpace::SetToZero(b);

            pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
        }
        else
        {
            TSparseSpace::SetToZero(Dx); //Dx=0.00;
            TSparseSpace::SetToZero(b);

            pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
        }

        // Debugging info
        EchoInfo(iteration_number);

        // Updating the results stored in the database
        UpdateDatabase(A, Dx, b, BaseType::MoveMeshFlag());

        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

        if (is_converged == true)
        {
            //initialisation of the convergence criteria
            mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(b);

                pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
        }

        //Iteration Cicle... performed only for NonLinearProblems
        while (is_converged == false &&
               iteration_number++ < mMaxIterationNumber)
        {
            //setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

            is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(Dx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (GetKeepSystemConstantDuringIterations() == false)
                    {
                        //A = 0.00;
                        TSparseSpace::SetToZero(A);
                        TSparseSpace::SetToZero(Dx);
                        TSparseSpace::SetToZero(b);

                        pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(Dx);
                        TSparseSpace::SetToZero(b);

                        pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(Dx);
                    TSparseSpace::SetToZero(b);

                    pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
                }
            }
            else
            {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            EchoInfo(iteration_number);

            // Updating the results stored in the database
            UpdateDatabase(A, Dx, b, BaseType::MoveMeshFlag());

            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

            ResidualIsUpdated = false;

            if (is_converged == true)
            {

                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(b);

                    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
                    ResidualIsUpdated = true;
                    //std::cout << "mb is calculated" << std::endl;
                }

                is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            MaxIterationsExceeded();
        }

        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (ResidualIsUpdated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(mb);
            //    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
        }

        //calculate reactions if required
        if (mCalculateReactionsFlag == true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), A, Dx, b);
        }

        return is_converged;
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
    TSystemMatrixType &GetSystemMatrix()
    {
        TSystemMatrixType &mA = *mpA;

        return mA;
    }

    /**
     * @brief This method directly sets the input as the LHS
     * @param A The LHS matrix
     */
    void GetDirectSystemMatrix(TSystemMatrixType& A)
    {
        A = *mpA;
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

    typename TLinearSolver::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employed
    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria; /// The pointer to the convergence criteria employed

    TSystemVectorPointerType mpDx; /// The incremement in the solution
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

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    unsigned int mMaxIterationNumber; /// The maximum number of iterations, 30 by default

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mKeepSystemConstantDuringIterations; // Flag to allow keeping system matrix constant during iterations

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief Here the database is updated
     * @param A The LHS matrix of the system of equations
     * @param Dx The incremement in the solution
     * @param b The RHS vector of the system of equations
     * @param MoveMesh The flag that allows to move the mesh
     */

    virtual void UpdateDatabase(
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b,
        const bool MoveMesh)
    {
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        pScheme->Update(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

        // Move the mesh if needed
        if (MoveMesh == true)
        {
            BaseType::MoveMesh();
        }
    }

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber)
    {
        TSystemMatrixType &A = *mpA;
        TSystemVectorType &Dx = *mpDx;
        TSystemVectorType &b = *mpb;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("Dx") << "Solution obtained = " << Dx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << b << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            //             KRATOS_INFO("Solution") << "After first system solution" << std::endl;
            KRATOS_INFO("LHS") << "SystemMatrix = " << A << std::endl;
            KRATOS_INFO("Dx") << "Solution obtained = " << Dx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << b << std::endl;
        }
        else if (this->GetEchoLevel() == 4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), A, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), b);
        }
    }

    /**
     * @brief This method prints information after reach the max number of iterations
     * @todo Replace by logger
     */

    virtual void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        mpConvergenceCriteria->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

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

    ResidualBasedNewtonRaphsonStrategy(const ResidualBasedNewtonRaphsonStrategy &Other){};

    ///@}

}; /* Class ResidualBasedNewtonRaphsonStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY  defined */
