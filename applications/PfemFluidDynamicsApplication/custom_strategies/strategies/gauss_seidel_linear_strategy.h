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
//

#if !defined(KRATOS_GAUSS_SEIDEL_LINEAR_STRATEGY)
#define KRATOS_GAUSS_SEIDEL_LINEAR_STRATEGY

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"
//#include "boost/timer.hpp"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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
    template <class TSparseSpace,
            class TDenseSpace,  //= DenseSpace<double>,
            class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
            >
    class GaussSeidelLinearStrategy
        : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        /** Counted pointer of ClassName */
        KRATOS_CLASS_POINTER_DEFINITION(GaussSeidelLinearStrategy);

        typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

        typedef typename BaseType::TDataType TDataType;

        typedef TSparseSpace SparseSpaceType;

        typedef typename BaseType::TSchemeType TSchemeType;
        typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

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

        //constructor specifying the builder and solver

        GaussSeidelLinearStrategy(
            ModelPart &model_part,
            typename TSchemeType::Pointer pScheme,
            typename TLinearSolver::Pointer pNewLinearSolver,
            typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
            bool ReformDofSetAtEachStep = true,
            bool CalculateNormDxFlag = false)
            : ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part)
        {
            KRATOS_TRY

            /* std::cout<<" GaussSeidelLinearStrategy"<<std::endl; */

            mReformDofSetAtEachStep = ReformDofSetAtEachStep;
            mCalculateNormDxFlag = CalculateNormDxFlag;

            //saving the scheme
            mpScheme = pScheme;

            //saving the linear solver
            mpLinearSolver = pNewLinearSolver;

            //setting up the  builder and solver
            mpBuilderAndSolver = pNewBuilderAndSolver;

            //set flag to start correcty the calculations
            mInitializeWasPerformed = false;

            //tells to the Builder And Solver if the system matrix and vectors need to
            //be reshaped at each step or not
            mReformDofSetAtEachStep = true;
            GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

            //set EchoLevel to the default value (only time is displayed)
            this->SetEchoLevel(1);

            //by default the matrices are rebuilt at each solution step
            BaseType::SetRebuildLevel(1);

            KRATOS_CATCH("")
        }

        /** Destructor.
     */
        virtual ~GaussSeidelLinearStrategy()
        {
            // in trilinos third party library, the linear solver's
            // preconditioner should be freed before the system matrix.
            // we control the deallocation order with Clear().
            this->Clear();
        }

        /** Destructor.
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

        void SetReformDofSetAtEachStepFlag(bool flag)
        {
            mReformDofSetAtEachStep = flag;
            mReformDofSetAtEachStep = true;
            GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
        }

        bool GetReformDofSetAtEachStepFlag()
        {
            return mReformDofSetAtEachStep;
        }

        //level of echo for the solving strategy
        // 0 -> mute... no echo at all
        // 1 -> printing time and basic informations
        // 2 -> printing linear solver data
        // 3 -> Print of debug informations:
        //		Echo of stiffness matrix, Dx, b...

        void SetEchoLevel(int Level) override
        {
            BaseType::SetEchoLevel(Level);
            GetBuilderAndSolver()->SetEchoLevel(Level);
        }

        //*********************************************************************************
        /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/

        /**
    operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
        /* void Predict() */
        /* { */
        /*     KRATOS_TRY */
        /*     //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions */
        /*     //if the operations needed were already performed this does nothing */
        /*     //if(mInitializeWasPerformed == false) */
        /*     //{ */
        /*     //	Initialize(); */
        /*     //	mInitializeWasPerformed = true; */
        /*     //} */

        /*     ////initialize solution step */
        /*     //if (mSolutionStepIsInitialized == false) */
        /*     //	InitializeSolutionStep(); */

        /*     TSystemMatrixType& mA = *mpA; */
        /*     TSystemVectorType& mDx = *mpDx; */
        /*     TSystemVectorType& mb = *mpb; */

        /*     DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet(); */

        /*     this->GetScheme()->Predict(BaseType::GetModelPart(), rDofSet, mA, mDx, mb); */

        /*     KRATOS_CATCH("") */
        /* } */

        //*********************************************************************************
        /**
    the problem of interest is solved
    a double containing norm(Dx) is returned if CalculateNormDxFlag == true, else 0 is returned
     */
        //**********************************************************************

        double Solve() override
        {
            KRATOS_TRY

            //pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

            TSystemMatrixType &mA = *mpA;
            TSystemVectorType &mDx = *mpDx;
            TSystemVectorType &mb = *mpb;

            if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false)
            {
                TSparseSpace::SetToZero(mA);
                TSparseSpace::SetToZero(mDx);
                TSparseSpace::SetToZero(mb);
                // passing smart pointers instead of references here
                // to prevent dangling pointer to system matrix when
                // reusing ml preconditioners in the trilinos tpl
                pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                BaseType::mStiffnessMatrixIsBuilt = true;
            }
            else
            {
                TSparseSpace::SetToZero(mDx);
                TSparseSpace::SetToZero(mb);
                pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
            }

            if (BaseType::GetEchoLevel() == 3) //if it is needed to print the debug info
            {
                std::cout << "SystemMatrix = " << mA << std::endl;
                std::cout << "solution obtained = " << mDx << std::endl;
                std::cout << "RHS  = " << mb << std::endl;
            }

            if (this->GetEchoLevel() == 4) //print to matrix market file
            {
                std::stringstream matrix_market_name;
                matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm";
                TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), mA, false);

                std::stringstream matrix_market_vectname;
                matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
                TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), mb);
            }

            //update results
            DofsArrayType &rDofSet = pBuilderAndSolver->GetDofSet();
            pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            //move the mesh if needed
            /* BaseType::MoveMesh();  */

            //calculate if needed the norm of Dx
            double normDx = 0.00;
            if (mCalculateNormDxFlag == true)
            {
                normDx = TSparseSpace::TwoNorm(mDx);
            }

            //Finalisation of the solution step,
            //operations to be done after achieving convergence, for example the
            //Final Residual Vector (mb) has to be saved in there
            //to avoid error accumulation
            pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
            pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            //Cleaning memory after the solution
            pScheme->Clean();

            return normDx;

            KRATOS_CATCH("")
        }

        TSystemMatrixType &GetSystemMatrix() override
        {
            TSystemMatrixType &mA = *mpA;

            return mA;
        }
        //*********************************************************************************

        double GetResidualNorm() override
        {
            if (TSparseSpace::Size(*mpb) != 0)
                return TSparseSpace::TwoNorm(*mpb);
            else
                return 0.0;
        }

        //*********************************************************************************

        /**
    this operations should be called before printing the results when non trivial results (e.g. stresses)
    need to be calculated given the solution of the step

      This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
        void CalculateOutputData() override
        {
            TSystemMatrixType &mA = *mpA;
            TSystemVectorType &mDx = *mpDx;
            TSystemVectorType &mb = *mpb;

            DofsArrayType &rDofSet = GetBuilderAndSolver()->GetDofSet();
            GetScheme()->CalculateOutputData(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        }

        //*********************************************************************************

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

        /*@} */
        /**@name Inquiry */
        /*@{ */

        /*@} */
        /**@name Friends */
        /*@{ */

        /*@} */

    protected:
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

    private:
        /**@name Static Member Variables */
        /*@{ */

        /*@} */
        /**@name Member Variables */
        /*@{ */

        typename TLinearSolver::Pointer mpLinearSolver;

        typename TSchemeType::Pointer mpScheme;

        typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

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

        //calculates if required the norm of the correction term Dx
        bool mCalculateNormDxFlag;

        /**
    Flag telling if it is needed or not to compute the reactions

      default = true
     */
        bool mInitializeWasPerformed;

        unsigned int mMaxIterationNumber;

        /*@} */
        /**@name Private Operators*/
        /*@{ */
        //**********************************************************************
        //**********************************************************************

        void Initialize() override
        {
            KRATOS_TRY

            if (BaseType::GetEchoLevel() > 2)
                std::cout << "entering in the  Initialize of the GaussSeidelLinearStrategy" << std::endl;

            //pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (pScheme->SchemeIsInitialized() == false)
                pScheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (pScheme->ElementsAreInitialized() == false)
                pScheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (pScheme->ConditionsAreInitialized() == false)
                pScheme->InitializeConditions(BaseType::GetModelPart());

            if (BaseType::GetEchoLevel() > 2)
                std::cout << "exiting the  Initialize of the GaussSeidelLinearStrategy" << std::endl;

            KRATOS_CATCH("")
        }

        //**********************************************************************
        //**********************************************************************

        void InitializeSolutionStep() override
        {
            KRATOS_TRY

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
            typename TSchemeType::Pointer pScheme = GetScheme();
            //int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            /* ProcessInfo& pCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo(); */

            //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
            //if the operations needed were already performed this does nothing
            if (mInitializeWasPerformed == false)
            {
                Initialize();
                mInitializeWasPerformed = true;
            }

            /* std::cout << "Gauss-Seidel VP Strategy,  CurrentTime: " << pCurrentProcessInfo[TIME] << std::endl; */

            //loop to reform the dofset
            // boost::timer system_construction_time;

            //setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());

            //shaping correctly the system
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());

            //setting up the Vectors involved to the correct size
            pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA, mpDx, mpb, BaseType::GetModelPart());
            //if (BaseType::GetEchoLevel() > 1 && rank == 0)
            //  std::cout << "System Construction Time : " << system_construction_time.elapsed() << std::endl;

            TSystemMatrixType &mA = *mpA;
            TSystemVectorType &mDx = *mpDx;
            TSystemVectorType &mb = *mpb;

            //initial operations ... things that are constant over the Solution Step
            pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            //initial operations ... things that are constant over the Solution Step
            pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            KRATOS_CATCH("")
        }

        //**********************************************************************
        //**********************************************************************

        void MaxIterationsExceeded()
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }

        //**********************************************************************
        //**********************************************************************

        void Clear() override
        {
            KRATOS_TRY;
            // if the preconditioner is saved between solves, it
            // should be cleared here.
            GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

            if (mpA != NULL)
                SparseSpaceType::Clear(mpA);
            if (mpDx != NULL)
                SparseSpaceType::Clear(mpDx);
            if (mpb != NULL)
                SparseSpaceType::Clear(mpb);

            //setting to zero the internal flag to ensure that the dof sets are recalculated
            GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
            GetBuilderAndSolver()->Clear();
            GetScheme()->Clear();

            KRATOS_CATCH("");
        }

        /**
     * function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
        int Check() override
        {
            KRATOS_TRY

            BaseType::Check();

            GetBuilderAndSolver()->Check(BaseType::GetModelPart());

            GetScheme()->Check(BaseType::GetModelPart());

            return 0;

            KRATOS_CATCH("")
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
        GaussSeidelLinearStrategy(const GaussSeidelLinearStrategy &Other);

        /*@} */

    }; /* Class GaussSeidelLinearStrategy */

    /*@} */

    /**@name Type Definitions */
    /*@{ */

    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_GAUSS_SEIDEL_LINEAR_STRATEGY  defined */
