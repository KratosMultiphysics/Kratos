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


#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY )
#define  KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY

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

/// Short class definition.

/// Detail class definition.

//URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

//URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

//URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

//URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


//URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

//URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

//URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

//URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


 
template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
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
     * Constructor.
     */
     
    ResidualBasedNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

        mKeepSystemConstantDuringIterations = false;

        //set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;


        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        //saving the scheme
        mpScheme = pScheme;

        //saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        //setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                 new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (mpLinearSolver)
                             );

        //set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        //tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        //tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        //set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        //by default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    ResidualBasedNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

        mKeepSystemConstantDuringIterations = false;

        //set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;


        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        //saving the scheme
        mpScheme = pScheme;

        //saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        //setting up the default builder and solver
        mpBuilderAndSolver = pNewBuilderAndSolver;

        //set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        //tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        //tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        //set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        //by default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    /** 
     * Destructor.
     */
     
    virtual ~ResidualBasedNewtonRaphsonStrategy()
    {
    }

    /**
     * Destructor.
     */
     
    //Set and Get Scheme ... containing Builder, Update and other

    void SetScheme(typename TSchemeType::Pointer pScheme) 
    {
        mpScheme = pScheme;
    };

    typename TSchemeType::Pointer GetScheme() 
    {
        return mpScheme;
    };

    //Set and Get the BuilderAndSolver

    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver) 
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver() 
    {
        return mpBuilderAndSolver;
    };

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

    //level of echo for the solving strategy
    // 0 -> mute... no echo at all
    // 1 -> printing time and basic informations
    // 2 -> printing linear solver data
    // 3 -> Print of debug informations:
    //		Echo of stiffness matrix, Dx, b...

    virtual void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * Operation to predict the solution ... if it is not called a trivial predictor is used in which the
     * values of the solution step of interest are assumed equal to the old values 
     */
     
    virtual void Predict() override
    {
        KRATOS_TRY
        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
            Initialize();

        //initialize solution step
        if (mSolutionStepIsInitialized == false)
            InitializeSolutionStep();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;


        GetScheme()->Predict(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), mA, mDx, mb);

        //move the mesh if needed
        if (this->MoveMeshFlag() == true) BaseType::MoveMesh();

        KRATOS_CATCH("")
    }

    /**
     * Initialization of member variables and prior operations
     */
     
    virtual void Initialize() override
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false)
        {
            //pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (pScheme->SchemeIsInitialized() == false)
                pScheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (pScheme->ElementsAreInitialized() == false)
                pScheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (pScheme->ConditionsAreInitialized() == false)
                pScheme->InitializeConditions(BaseType::GetModelPart());

            //initialisation of the convergence criteria
            if (mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
                mpConvergenceCriteria->Initialize(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * The problem of interest is solved.
     * This function calls sequentially: Initialize(), InitializeSolutionStep(), Predict(),
     * SolveSolutionStep() and FinalizeSolutionStep().
     * All those functions can otherwise be called separately.
     */
     
    virtual double Solve() override
    {
        Initialize();
        InitializeSolutionStep();
        Predict();
        SolveSolutionStep();
        FinalizeSolutionStep();
        return 0.00;
    }

    /**
     * Clears the internal storage
     */
     
    virtual void Clear() override
    {
        KRATOS_TRY

        SparseSpaceType::Clear(mpA);
        TSystemMatrixType& mA = *mpA;
        SparseSpaceType::Resize(mA, 0, 0);

        SparseSpaceType::Clear(mpDx);
        TSystemVectorType& mDx = *mpDx;
        SparseSpaceType::Resize(mDx, 0);

        SparseSpaceType::Clear(mpb);
        TSystemVectorType& mb = *mpb;
        SparseSpaceType::Resize(mb, 0);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();

        GetScheme()->Clear();

        if (this->GetEchoLevel() > 0)
        {
            if(BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                std::cout << "Newton Raphson strategy Clear function used" << std::endl;
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * This should be considered as a "post solution" convergence check which is useful for coupled 
     * analysis - the convergence criteria used is the one used inside the "solve" step
     */
     
    virtual bool IsConverged() override
    {
        KRATOS_TRY

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;


        if (mpConvergenceCriteria->mActualizeRHSIsNeeded == true)
        {
            GetBuilderAndSolver()->BuildRHS(GetScheme(), BaseType::GetModelPart(), mb);
        }

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), mA, mDx, mb);
        KRATOS_CATCH("")
    }

    /**
     * This operations should be called before printing the results when non trivial results 
     * (e.g. stresses)
     * Need to be calculated given the solution of the step
     * This operations should be called only when needed, before printing as it can involve a non
     * negligible cost
     */
     
    virtual void CalculateOutputData() override
    {
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        GetScheme()->CalculateOutputData(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), mA, mDx, mb);
    }

    /**
     * Performs all the required operations that should be done (for each step) 
     * before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    
    virtual void InitializeSolutionStep() override
    {
        KRATOS_TRY

        //pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

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
                std::cout << "setup_dofs_time : " << setup_dofs_endtime- setup_dofs_begintime << std::endl;
            }

            //shaping correctly the system
            double setup_system_begin = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double setup_system_end = OpenMPUtils::GetCurrentTime();
                std::cout << rank << ": setup_system_time : " << setup_system_end- setup_system_begin << std::endl;
            }
        }

        //prints informations about the current time
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << " " << std::endl;
            std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        if (mSolutionStepIsInitialized == false)
        {
            /*typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
            typename TSchemeType::Pointer pScheme = GetScheme();
            int rank = BaseType::GetModelPart().GetCommunicator().MyPID(); */

            //setting up the Vectors involved to the correct size
            double system_matrix_resize_begin = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->ResizeAndInitializeVectors(mpA, mpDx, mpb, BaseType::GetModelPart().Elements(), BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double system_matrix_resize_end = OpenMPUtils::GetCurrentTime();
                std::cout << rank << ": system_matrix_resize_time : " << system_matrix_resize_end- system_matrix_resize_begin << std::endl;
            }
            TSystemMatrixType& mA = *mpA;
            TSystemVectorType& mDx = *mpDx;
            TSystemVectorType& mb = *mpb;

            //initial operations ... things that are constant over the Solution Step
            pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            //initial operations ... things that are constant over the Solution Step
            pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * Performs all the required operations that should be done (for each step) 
     * after solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    
    virtual void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

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

        KRATOS_CATCH("")
    }

    /**
     * Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    
    virtual bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        //initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        //			BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);
        bool is_converged = false;
        bool ResidualIsUpdated = false;
        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), mA, mDx, mb);

        //function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        else
        {
            TSparseSpace::SetToZero(mDx); //mDx=0.00;
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }

        if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            // 				std::cout << "After first system solution" << std::endl;
            std::cout << "SystemMatrix = " << mA << std::endl;
            std::cout << "solution obtained = " << mDx << std::endl;
            std::cout << "RHS  = " << mb << std::endl;
        }
        if (this->GetEchoLevel() == 4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), mA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), mb);
        }
        
        // Updating the results stored in the database
        UpdateDatabase(mA, mDx, mb, BaseType::MoveMeshFlag());

        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        if (is_converged == true)
        {
            //initialisation of the convergence criteria
            mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), mA, mDx, mb);

            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(mb);

                pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), mA, mDx, mb);
        }


        //Iteration Cicle... performed only for NonLinearProblems
        while (is_converged == false &&
                iteration_number++<mMaxIterationNumber)
        {
            //setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), mA, mDx, mb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(mDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false )
                {
                    if( GetKeepSystemConstantDuringIterations() == false)
                    {
                        //mA = 0.00;
                        TSparseSpace::SetToZero(mA);
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                }
            }
            else
            {
                std::cout << "ATTENTION: no free DOFs!! " << std::endl;
            }

            
            // Updating the results stored in the database
            UpdateDatabase(mA, mDx, mb, BaseType::MoveMeshFlag());

            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            ResidualIsUpdated = false;

            if (is_converged == true)
            {

                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
                    ResidualIsUpdated = true;
                    //std::cout << "mb is calculated" << std::endl;
                }

                is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), mA, mDx, mb);
            }
        }


        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            MaxIterationsExceeded();

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
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
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

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
    }

    void SetKeepSystemConstantDuringIterations(bool value)
    {
        mKeepSystemConstantDuringIterations = value;
    }

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

    typename TSchemeType::Pointer mpScheme;

    typename TLinearSolver::Pointer mpLinearSolver;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;

//     TSystemVectorType mDx;
//     TSystemVectorType mb;
//     TSystemMatrixType mA;
                    
    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb;
    TSystemMatrixPointerType mpA;

    /**
     * Flag telling if it is needed to reform the DofSet at each
     * solution step or if it is possible to form it just once
     * - true  => reforme at each time step
     * - false => form just one (more efficient)
     * Default = false
     */
     
    bool mReformDofSetAtEachStep;

    /**
     * Flag telling if it is needed or not to compute the reaction
     * default = true
     * */
     
    bool mCalculateReactionsFlag;

    bool mSolutionStepIsInitialized;

    ///default = 30
    unsigned int mMaxIterationNumber;

    bool mInitializeWasPerformed;

    //flag to allow keeping system matrix constant during iterations
    bool mKeepSystemConstantDuringIterations;


    ///@}
    ///@name Private Operators
    ///@{
    
    /**
     * Here the database is updated
     */
     
    virtual void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh
    )
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

    //**********************************************************************
    //**********************************************************************

    void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }

    }

    /**
     * function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
     
    virtual int Check() override
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
     
    ResidualBasedNewtonRaphsonStrategy(const ResidualBasedNewtonRaphsonStrategy& Other)
    {
    };


    ///@}

}; /* Class ResidualBasedNewtonRaphsonStrategy */

///@}

///@name Type Definitions 
///@{


///@}

} /* namespace Kratos. */

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY  defined */

