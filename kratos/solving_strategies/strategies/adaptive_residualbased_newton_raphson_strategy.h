//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//



#if !defined(KRATOS_ADAPTIVE_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY )
#define  KRATOS_ADAPTIVE_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY


/* System includes */


/* External includes */


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/matrix_market_interface.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

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
/**   Detail class definition.

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
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class AdaptiveResidualBasedNewtonRaphsonStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( AdaptiveResidualBasedNewtonRaphsonStrategy );

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
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


    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    AdaptiveResidualBasedNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        int MinIterations = 4,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        double ReductionFactor = 0.5,
        double IncreaseFactor = 1.3,
        int NumberOfCycles = 5
    )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY
        //set flags to default values
        SetMaxIterationNumber(MaxIterations);

        mMinIterationNumber = MinIterations;

        mReductionFactor = ReductionFactor;

        mIncreaseFactor = IncreaseFactor;

        mNumberOfCycles = NumberOfCycles;

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
                                 new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(mpLinearSolver)
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

        KRATOS_WATCH("AdaptiveResidualBasedNewtonRaphsonStrategy is chosen");

        KRATOS_CATCH("")
    }

    AdaptiveResidualBasedNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        int MinIterations = 4,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        double ReductionFactor = 0.5,
        double IncreaseFactor = 1.3,
        int NumberOfCycles = 5
    )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY
        //set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;

        mMinIterationNumber = MinIterations;

        mReductionFactor = ReductionFactor;

        mIncreaseFactor = IncreaseFactor;

        mNumberOfCycles = NumberOfCycles;


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

        KRATOS_CATCH("")
    }

    /** Destructor.
    */
    ~AdaptiveResidualBasedNewtonRaphsonStrategy() override {}

    /** Destructor.
    */

    //Set and Get Scheme ... containing Builder, Update and other
    void SetScheme(typename TSchemeType::Pointer pScheme )
    {
        mpScheme = pScheme;
    };
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    //Set and Get the BuilderAndSolver
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver )
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

    void SetMaxIterationNumber( int  MaxIterationNumber)
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
    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/

    /**
    operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
    */
    void Predict() override
    {
        KRATOS_TRY
        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if(mInitializeWasPerformed == false)
            Initialize();

        //initialize solution step
        if (mSolutionStepIsInitialized == false)
            InitializeSolutionStep();

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;


        GetScheme()->Predict(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

        //move the mesh if needed
        if(this->MoveMeshFlag() == true) BaseType::MoveMesh();

        KRATOS_CATCH("")
    }

    //*********************************************************************************
    /**
    the problem of interest is solved
    */
    //**********************************************************************
    bool CyclicSolve()
    {
        KRATOS_TRY


        //pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        //int solstep = pCurrentProcessInfo.GetCurrentSolutionStep();
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if(mInitializeWasPerformed == false)
            Initialize();

        //set up the system, operation performed just once unless it is required
        //to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true )
        {
            //setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme,BaseType::GetModelPart());

            //shaping correctly the system
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
        }

        //prints informations about the current time
        if (this->GetEchoLevel()!=0)
        {
            std::cout << " " << std::endl;
            std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        //updates the database with a prediction of the solution
        Predict();

        //initialize solution step
        if (mSolutionStepIsInitialized == false)
            InitializeSolutionStep();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;



        //initializing the parameters of the Newton-Raphson cicle
        int iteration_number=1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
//			BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);
        bool is_converged = false;
//        bool ResidualIsUpdated = false;
        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
        is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

        //function to perform the building and the solving phase.
        if(BaseType::mRebuildLevel >1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
        }
        else
        {
            TSparseSpace::SetToZero(mDx); //mDx=0.00;
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildRHSAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
        }

        if (this->GetEchoLevel()==3) //if it is needed to print the debug info
        {
// 				std::cout << "After first system solution" << std::endl;
            std::cout << "SystemMatrix = " << mA << std::endl;
            std::cout << "solution obtained = " << mDx << std::endl;
            std::cout << "RHS  = " << mb << std::endl;
        }
        if (this->GetEchoLevel()==4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_"<< BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number  << ".mm";
            WriteMatrixMarketMatrix((char*)(matrix_market_name.str()).c_str() , mA, false);
        }

        //update results
        pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);

        if (is_converged==true)
        {
            //initialisation of the convergence criteria
            mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(mb);

                pBuilderAndSolver->BuildRHS(pScheme,BaseType::GetModelPart(),mb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
        }


        //Iteration Cicle... performed only for NonLinearProblems
        while(	is_converged == false &&
                iteration_number++<mMaxIterationNumber)
        {
            //setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);

            is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(mDx)!=0)
            {
                if(BaseType::mRebuildLevel >1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    //mA = 0.00;
                    TSparseSpace::SetToZero(mA);
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
                }
                else
                {
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHSAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
                }
            }
            else
            {
                std::cout << "ATTENTION: no free DOFs!! " << std::endl;
            }


            //Updating the results stored in the database
            pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);

//            ResidualIsUpdated = false;

            if (is_converged==true)
            {

                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHS(pScheme,BaseType::GetModelPart(),mb);
//                    ResidualIsUpdated = true;
                    //std::cout << "mb is calculated" << std::endl;
                }

                is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
            }
        }

        return is_converged;

        KRATOS_CATCH("")

    }
    //**********************************************************************
    double Solve() override
    {
        KRATOS_TRY
        bool is_converged = false;
        int number_of_cycles = 0;

        while(is_converged == false && number_of_cycles++ < mNumberOfCycles )
        {
            KRATOS_WATCH("****************** inside cyclic_lopp*******************");
            KRATOS_WATCH(mNumberOfCycles);
            KRATOS_WATCH(mMinIterationNumber);
            KRATOS_WATCH(mReductionFactor);
            KRATOS_WATCH(mIncreaseFactor);

            // mSolutionStepIsInitialized = false;
            //Iteration Cicle... performed only for NonLinearProblems
            is_converged = CyclicSolve();

            double iteration_number = BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER];

            if(is_converged == false)
            {
                ModelPart& model_part = BaseType::GetModelPart();
                double CurrentTime = model_part.GetProcessInfo()[TIME];
                KRATOS_WATCH("old delta_time");
                KRATOS_WATCH(model_part.GetProcessInfo()[DELTA_TIME]);
                double delta_time  = model_part.GetProcessInfo()[DELTA_TIME];
                double new_time =  CurrentTime - delta_time * mReductionFactor;


                ReductionStrategy(model_part, new_time);


                //plots a warning if the maximum number of iterations is exceeded
                MaxIterationsExceeded(number_of_cycles);

                KRATOS_WATCH("CurrentTime");
                KRATOS_WATCH(CurrentTime);

                KRATOS_WATCH("new_time");
                KRATOS_WATCH(new_time);

                KRATOS_WATCH("new delta_time");
                KRATOS_WATCH(model_part.GetProcessInfo()[DELTA_TIME]);
            }

            if(iteration_number <=  mMinIterationNumber)
            {
                ModelPart& model_part = BaseType::GetModelPart();
                double DeltaTime = model_part.GetProcessInfo()[DELTA_TIME];
                double New_Dt = DeltaTime * mIncreaseFactor;
                model_part.GetProcessInfo()[DELTA_TIME] = New_Dt;
                //ModelPart.GetProcessInfo().SetCurrentTime(NewTime);
                PrintIncreaseOfDeltaTime(DeltaTime, New_Dt);
            }
        }//end of cyclic loop

        BaseType::GetModelPart().GetProcessInfo()[SCALE] = number_of_cycles;
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;



        //recalculate residual if needed
        // (note that some convergence criteria need it to be recalculated)
        if (mpConvergenceCriteria->GetActualizeRHSflag() == false)
        {
            TSparseSpace::SetToZero(mb);

            mpBuilderAndSolver->BuildRHS(mpScheme,BaseType::GetModelPart(),mb);

            //std::cout << "mb is calculated" << std::endl;
        }

        //calculate reactions if required
        if (mCalculateReactionsFlag ==true)
        {
            mpBuilderAndSolver->CalculateReactions(mpScheme,BaseType::GetModelPart(),mA,mDx,mb);
        }
        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation
        mpScheme->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
        mpBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);

        //Cleaning memory after the solution
        mpScheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            SparseSpaceType::Clear(mpA);
            SparseSpaceType::Clear(mpDx);
            SparseSpaceType::Clear(mpb);

            this->Clear();
        }

        return 0.00;

        KRATOS_CATCH("")

    }

    /**
    this should be considered as a "post solution" convergence check which is useful for coupled analysis
    - the convergence criteria used is the one used inside the "solve" step
    */
    //**********************************************************************
    bool IsConverged() override
    {
        KRATOS_TRY

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;


        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            GetBuilderAndSolver()->BuildRHS(GetScheme(),BaseType::GetModelPart(),mb);
        }

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
        KRATOS_CATCH("")

    }

    //*********************************************************************************
    /**
    this operations should be called before printing the results when non trivial results (e.g. stresses)
    need to be calculated given the solution of the step

    This operations should be called only when needed, before printing as it can involve a non negligible cost
    */
    void CalculateOutputData() override
    {
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
        GetScheme()->CalculateOutputData(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
    }

    //**********************************************************************
    //**********************************************************************
    void Clear() override
    {
        KRATOS_TRY
        std::cout << "Newton Raphson strategy Clear function used" << std::endl;

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Resize(mA,0,0);

        SparseSpaceType::Clear(mpDx);
        SparseSpaceType::Resize(mDx,0);

        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Resize(mb,0);


        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();


        KRATOS_CATCH("");
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

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
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

    typename TSchemeType::Pointer mpScheme;

    typename TLinearSolver::Pointer mpLinearSolver;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;

    /*		TSystemVectorType mDx;
    		TSystemVectorType mb;
    		TSystemMatrixType mA;*/
    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb;
    TSystemMatrixPointerType mpA;

    /**
    Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    - true  => reforme at each time step
    - false => form just one (more efficient)

    Default = false
    */
    bool mReformDofSetAtEachStep;

    /**
    Flag telling if it is needed or not to compute the reactions

    default = true
    */
    bool mCalculateReactionsFlag;

    bool mSolutionStepIsInitialized;

    ///default = 30
    int mMaxIterationNumber;
    int mMinIterationNumber;
    bool mInitializeWasPerformed;

    double mReductionFactor;
    double mIncreaseFactor;

    int mNumberOfCycles;

    /*@} */
    /**@name Private Operators*/
    /*@{ */
    //**********************************************************************
    //**********************************************************************
    void Initialize() override
    {
        KRATOS_TRY

        //pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

        //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
        if (pScheme->SchemeIsInitialized() == false)
            pScheme->Initialize(BaseType::GetModelPart());

        //Initialize The Elements - OPERATIONS TO BE DONE ONCE
        if (pScheme->ElementsAreInitialized() == false)
            pScheme->InitializeElements(BaseType::GetModelPart());

        //initialisation of the convergence criteria
        if (mpConvergenceCriteria->IsInitialized() == false)
            mpConvergenceCriteria->Initialize(BaseType::GetModelPart());


        mInitializeWasPerformed = true;

        KRATOS_CATCH("")
    }


    //**********************************************************************
    //**********************************************************************
    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = GetScheme();



        //setting up the Vectors involved to the correct size
        pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA,mpDx,mpb,BaseType::GetModelPart());

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;


        //initial operations ... things that are constant over the Solution Step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);

        //initial operations ... things that are constant over the Solution Step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);

        mSolutionStepIsInitialized = true;

        KRATOS_CATCH("")
    }


    //**********************************************************************
    //**********************************************************************
    void MaxIterationsExceeded(int cycle_number)
    {
        std::cout << "***************************************************" << std::endl;
        std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
        std::cout << "***************************************************" << std::endl;

        std::cout << "*****************REDUCTION CYCLE********************" << std::endl;
        std::cout <<                     cycle_number                       << std::endl;
        std::cout << "***************************************************" << std::endl;
    }


    void PrintIncreaseOfDeltaTime(double old_dt , double new_dt)
    {
        std::cout << "***************************************************" << std::endl;
        std::cout << "******* ATTENTION: INcreasing dt ********" << std::endl;
        std::cout << "***************************************************" << std::endl;

        std::cout << "*****************OLD DT********************" << std::endl;
        std::cout <<                   old_dt                      << std::endl;
        std::cout << "*****************NEW DT********************" << std::endl;
        std::cout <<                   new_dt                      << std::endl;
    }

    //**********************************************************************
    //**********************************************************************
    void ReductionStrategy(ModelPart& model_part, double NewTime)
    {
        model_part.ReduceTimeStep(model_part, NewTime);

        for(ModelPart::NodeIterator i = model_part.NodesBegin() ;
                i != model_part.NodesEnd() ; ++i)
        {
            (i)->X() = (i)->X0() + i->GetSolutionStepValue(DISPLACEMENT_X , 1);
            (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(DISPLACEMENT_Y , 1);
            (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(DISPLACEMENT_Z , 1);
        }

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
    AdaptiveResidualBasedNewtonRaphsonStrategy(const AdaptiveResidualBasedNewtonRaphsonStrategy& Other);


    /*@} */

}; /* Class ResidualBasedNewtonRaphsonStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY  defined */

