// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Long Chen
//

//#if !defined(KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY )
//#define KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY
#if !defined(KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY )
#define  KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY
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
//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

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
    class FormfindingUpdatedReferenceStrategy
        : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        ///@name Type Definitions
        ///@{
        typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

        // Counted pointer of ClassName
        KRATOS_CLASS_POINTER_DEFINITION(FormfindingUpdatedReferenceStrategy);

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

        FormfindingUpdatedReferenceStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            typename TLinearSolver::Pointer pNewLinearSolver,
            typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
            int Maxiterations = 30,
            bool CalculateReactions = false,
            bool ReformDofSetAtEachStep = false,
            bool MoveMeshFlag = false
            )
            : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
        {
            KRATOS_TRY;

            std::cout << "Formfinding strategy created!" << std::endl;

            mKeepSystemConstantDuringIterations = false;

            // set flags to default values
            SetMaxIterationNumber(Maxiterations);
            mCalculateReactionsFlag = CalculateReactions;

            mReformDofSetAtEachStep = ReformDofSetAtEachStep;

            // saving the convergence criteria to be used
            mpConvergenceCriteria = pNewConvergenceCriteria;

            // saving the scheme
            mpScheme = pScheme;

            // saving the linear solver
            mpLinearSolver = pNewLinearSolver;

            // Setting up the default builder and solver
            mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                (
                    new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(mpLinearSolver)
                    );

            // set flags to start correctly the calculations
            mSolutionStepIsInitialized = false;
            mInitializeWasPerformed = false;

            // tells to the builder and solver if the reactions have to be calculated or not
            GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

            // tells to the Builder and Solver if the system matrix and vectors need to be reshaped at each step or not
            GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

            //set EchoLevel to the default value (only time is displayed)
            SetEchoLevel(1);

            //by default the matrices are rebuilt at each iteration
            this->SetRebuildLevel(2);

            mpA = TSparseSpace::CreateEmptyMatrixPointer();
            mpDx = TSparseSpace::CreateEmptyVectorPointer();
            mpb = TSparseSpace::CreateEmptyVectorPointer();

            mHomotopyFactor = 0.0;

            KRATOS_CATCH("");
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
            bool MoveMeshFlag = false
            )
            : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
        {
            KRATOS_TRY;

            std::cout << "Formfinding strategy created!" << std::endl;

            mKeepSystemConstantDuringIterations = false;

            // set flags to default values
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

            mHomotopyFactor = 0.0;

            KRATOS_CATCH("");
        }



        /**
        * Destructor.
        */

        ~FormfindingUpdatedReferenceStrategy() override
        {
        }

        /**
        * Destructor.
        */

        // Set and Get Scheme ... containing Builder, Update and other

        void SetScheme(typename TSchemeType::Pointer SetScheme)
        {
            mpScheme = SetScheme;
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
        }

        // Set and Get InitializePerformedFlag

        void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
        {
            mInitializeWasPerformed = InitializePerformedFlag;
        }

        bool GetInitializePerformedFlag()
        {
            return mInitializeWasPerformed;
        }

        // Set and Get CalculateReactionsFlag

        void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
        {
            mCalculateReactionsFlag = CalculateReactionsFlag;
        }

        bool GetCalculateReactionsFlag()
        {
            return mCalculateReactionsFlag;
        }

        // Set and Get MaxIterationNumber

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

        void SetEchoLevel(int Level) override
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

        void Predict() override
        {
            KRATOS_TRY
            // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
            // if the operations needed were already performed this does nothing
            if (mInitializeWasPerformed == false)
            {
                Initialize();
            }

            // initialize solution step
            if (mSolutionStepIsInitialized == false)
            {
                InitializeSolutionStep();
            }

            TSystemMatrixType& A = *mpA;
            TSystemVectorType& Dx = *mpDx;
            TSystemVectorType& b = *mpb;

            GetScheme()->Predict(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), A, Dx, b);

            //move the mesh if needed
            if (this->MoveMeshFlag() == true)
            {
                BaseType::MoveMesh();
            }

            KRATOS_CATCH("")
        }


        /**
        * Initialization of member variables and prior operations
        */

        void Initialize() override
        {
            KRATOS_TRY;

            if (mInitializeWasPerformed == false)
            {
                // pointers needed in the solution
                typename TSchemeType::Pointer pScheme = GetScheme();
                typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

                // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
                if (pScheme->SchemeIsInitialized() == false)
                {
                    pScheme->Initialize(BaseType::GetModelPart());
                }

                // Initialize The Elements - OPEARATIONS TO BE DONE ONCE
                if (pScheme->ElementsAreInitialized() == false)
                {
                    pScheme->InitializeElements(BaseType::GetModelPart());
                }

                // Initialize The Conditions - OPERATIONS TO BE DONE ONCE
                if (pScheme->ConditionsAreInitialized() == false)
                {
                    pScheme->InitializeConditions(BaseType::GetModelPart());
                }

                // Initialization of the convergence criteria
                if (mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
                {
                    mpConvergenceCriteria->Initialize(BaseType::GetModelPart());
                }

                mInitializeWasPerformed = true;

            }

            KRATOS_CATCH("");
        }

        /**
         * Teh problem of intereset is solved.
         * This function calls sequentially: Initialize(), InitializeSolutionStep(), Predict(),
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
         * Clears the internal storage
         */

        void Clear() override
        {
            KRATOS_TRY

            SparseSpaceType::Clear(mpA);
            TSystemMatrixType& A = *mpA;
            SparseSpaceType::Resize(A, 0, 0);

            SparseSpaceType::Clear(mpDx);
            TSystemVectorType& Dx = *mpDx;
            SparseSpaceType::Resize(Dx, 0);

            SparseSpaceType::Clear(mpb);
            TSystemVectorType& b = *mpb;
            SparseSpaceType::Resize(b, 0);

            //setting to zero the internal flag to ensure that the dof sets are recalculated
            GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
            GetBuilderAndSolver()->Clear();

            GetScheme()->Clear();

            if (this->GetEchoLevel() > 0)
            {
                if (BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
                {
                    std::cout << "Formfinding Updated Reference strategy Clear function used" << std::endl;
                }
            }

            KRATOS_CATCH("");

        }

        /**
         * This should be considered as a "post solution" convergence check which is useful for coupled
         * analysis - the convergence criteria used is the one used inside the "solve" step
         */

        bool IsConverged() override
        {
            KRATOS_TRY;

            TSystemMatrixType& A = *mpA;
            TSystemVectorType& Dx = *mpDx;
            TSystemVectorType& b = *mpb;

            if (mpConvergenceCriteria->mActualizeRHSIsNeeded == true)
            {
                GetBuilderAndSolver()->BuildRHS(GetScheme(), BaseType::GetModelPart(), b);
            }

            return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), A, Dx, b);

            KRATOS_CATCH("");
        }

        /**
         * This operations should be called before printing the results when non trivial results
         * (e.g. stresses)
         * Need to be calculated given the solution of the step
         * This operations should be called only when needed, before printing as it can involve a non
         * negligible cost
         */

        void CalculateOutputData() override
        {
            TSystemMatrixType& A = *mpA;
            TSystemVectorType& Dx = *mpDx;
            TSystemVectorType& b = *mpb;

            GetScheme()->CalculateOutputData(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), A, Dx, b);
        }

        /**
         * Performs all the required operations that should be done (for each step)
         * before solving the solution step.
         * A member variable should be used as a flag to make sure this function is called only once per step.
         */

        void InitializeSolutionStep() override
        {
            KRATOS_TRY;

            // pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
            int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            // set up the system, operation performed just once unless it is required
            // to reform the dof set at each iteration
            if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
            {
                //setting up the list of the DOFs to be solved
                double setup_dofs_begintime = OpenMPUtils::GetCurrentTime();
                pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());
                if (this->GetEchoLevel() > 0 && rank == 0)
                {
                    double setup_dofs_endtime = OpenMPUtils::GetCurrentTime();
                    std::cout << "setup_dofs_time : " << setup_dofs_endtime - setup_dofs_begintime << std::endl;
                }

                //shaping correctly the system
                double setup_system_begin = OpenMPUtils::GetCurrentTime();
                pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
                if (this->GetEchoLevel() > 0 && rank == 0)
                {
                    double setup_system_end = OpenMPUtils::GetCurrentTime();
                    std::cout << rank << ": setup_system_time : " << setup_system_end - setup_system_begin << std::endl;
                }
            }

            // prints information about the current time
            if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                std::cout << " " << std::endl;
                std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
            }

            if (mSolutionStepIsInitialized == false)
            {
                //setting up the Vectors involved to the correct size
                double system_matrix_resize_begin = OpenMPUtils::GetCurrentTime();
                pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA, mpDx, mpb, BaseType::GetModelPart().Elements(), BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
                //pBuilderAndSolver->ResizeAndInitializeVectors(mpA, mpDx, mpb, BaseType::GetModelPart().Elements(), BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
                if (this->GetEchoLevel() > 0 && rank == 0)
                {
                    double system_matrix_resize_end = OpenMPUtils::GetCurrentTime();
                    std::cout << rank << ": system_matrix_resize_time : " << system_matrix_resize_end - system_matrix_resize_begin << std::endl;
                }
                TSystemMatrixType& A = *mpA;
                TSystemVectorType& Dx = *mpDx;
                TSystemVectorType& b = *mpb;

                //initialize operations ... things that are constant over the Solution Step
                pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(), A, Dx, b);

                //initialize operations ... things that are constant over the Solution Step
                pScheme->InitializeSolutionStep(BaseType::GetModelPart(), A, Dx, b);

                ////***************************************************
                //// Update reference configuration for formfinding
                //// Initialize the Elements
                //pScheme->InitializeElements(BaseType::GetModelPart());

                mSolutionStepIsInitialized = true;
            }

            //***************************************************
            // Update reference configuration for formfinding
            // Initialize the Elements
            pScheme->InitializeElements(BaseType::GetModelPart());

            KRATOS_CATCH("");
        }

        /**
         * Performs all the required operations that should be done (for each step)
         * after solving the solution step.
         * A member variable should be used as a flag to make sure this function is called only once per step
         */

        void FinalizeSolutionStep() override
        {
            KRATOS_TRY;

            KRATOS_CATCH("");
        }


        /**
         * Solves the current step. This function returns true if a solution has been found, false otherwise
         */

        bool SolveSolutionStep() override
        {
            // Pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

            TSystemMatrixType& A = *mpA;
            TSystemVectorType& Dx = *mpDx;
            TSystemVectorType& b = *mpb;

            // Initializing the parameters of the Newton-Raphson cycle
            unsigned int iteration_number = 1;
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;      // long: this should be the iteration number of the newton raphson step
            bool is_converged = false;
            bool ResidualIsUpdated = false;
            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
            is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

            // function to perform the building and the solving phase
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

            // Iteration cycle... performed only for NonlinearProblems
            while (is_converged == false &&
                iteration_number++ < mMaxIterationNumber)
            {
                // setting the number of iteration
                BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

                pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

                is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

                // call the linear system solver to find the correction mDx for the system
                // it is not called if there is no system to solve
                if (SparseSpaceType::Size(Dx) != 0)
                {
                    if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                    {
                        if (GetKeepSystemConstantDuringIterations() == false)
                        {
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
                    std::cout << "ATTENTION: no free DOFs!! " << std::endl;
                } // checking if there are free DOFs

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
            } // end Nonlinear iteration cycles

            //plots a warning if the maximum number of iterations is exceeded
            if (iteration_number >= mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                MaxIterationsExceeded();
            }

            // recalculate residual if needed
            // (note that some convergence criteria need it to be recalculated)
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

        } // end SolveSolutionStep


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
        * homotopy factor describing the mapping from original problem to the
        * modified problem
        * 0 ~ only the modified problem is to be solved
        * 1 ~ only the original problem is to be solved
        */

        double mHomotopyFactor;


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

        virtual void EchoInfo(const unsigned int iteration_number)
        {
            TSystemMatrixType& A = *mpA;
            TSystemVectorType& Dx = *mpDx;
            TSystemVectorType& b = *mpb;

            if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
            {
                std::cout << "Solution obtained = " << Dx << std::endl;
                std::cout << "RHS  = " << b << std::endl;
            }
            else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
            {
                //             std::cout << "After first system solution" << std::endl;
                std::cout << "SystemMatrix = " << A << std::endl;
                std::cout << "Solution obtained = " << Dx << std::endl;
                std::cout << "RHS  = " << b << std::endl;
            }
            else if (this->GetEchoLevel() == 4) //print to matrix market file
            {
                std::stringstream matrix_market_name;
                matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
                TSparseSpace::WriteMatrixMarketMatrix((char*)(matrix_market_name.str()).c_str(), A, false);

                std::stringstream matrix_market_vectname;
                matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
                TSparseSpace::WriteMatrixMarketVector((char*)(matrix_market_vectname.str()).c_str(), b);
            }
        }


        //**********************************************************************
        //**********************************************************************

        void MaxIterationsExceeded()
        {
            if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
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

        FormfindingUpdatedReferenceStrategy(const FormfindingUpdatedReferenceStrategy& Other)
        {
        };


        ///@}
    }; /* Class FormfindingUpdatedReferenceStrategy */

       ///@}

       ///@name Type Definitions
       ///@{


       ///@}

} /* namespace Kratos. */

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY defined */
