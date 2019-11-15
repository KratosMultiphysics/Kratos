// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main author:    Manuel Messmer
//

#if !defined(KRATOS_PREBUCKLING_STRATEGY)
#define KRATOS_PREBUCKLING_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"

// Application includes
#include "structural_mechanics_application_variables.h"

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

/// Strategy for linearized prebuckling analysis
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver>
class PrebucklingStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PrebucklingStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType::Pointer SchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> ConvergenceCriteriaType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PrebucklingStrategy(
        ModelPart &rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pEigenSolver,
        BuilderAndSolverPointerType pBuilderAndSolver,
        typename ConvergenceCriteriaType::Pointer pConvergenceCriteria,
        int MaxIteration,
        double InitialStep,
        double SmallStep,
        double BigStep,
        double ConvergenceRatio )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpEigenSolver = pEigenSolver;

        mpBuilderAndSolver = pBuilderAndSolver;

        mpConvergenceCriteria = pConvergenceCriteria;

        mMaxIteration = MaxIteration;

        mInitialStep = InitialStep;

        mSmallStep = SmallStep;

        mBigStep = BigStep;

        mConvergenceRatio = ConvergenceRatio; 

        // Set Eigensolver flags
        mpEigenSolver->SetDofSetIsInitializedFlag(false);
        mpEigenSolver->SetReshapeMatrixFlag(false);
        mpEigenSolver->SetCalculateReactionsFlag(false);
        // Set Static Builder and Solver flags
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->SetCalculateReactionsFlag(false);
        mpBuilderAndSolver->SetReshapeMatrixFlag(false);

        // Set EchoLevel to the default value (only time is displayed)
        this->SetEchoLevel(1);

        // default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);

        //Set Matrices and Vectors to empty pointers
        mpStiffnessMatrix = TSparseSpace::CreateEmptyMatrixPointer();
        mpStiffnessMatrixPrevious = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpRHS = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    PrebucklingStrategy(const PrebucklingStrategy &Other) = delete;

    /// Destructor.
    ~PrebucklingStrategy() override
    {
        // Clear() controls order of deallocation to avoid invalid memory access
        // in some special cases.
        // warning: BaseType::GetModelPart() may be invalid here.
        this->Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetIsInitialized(bool val)
    {
        mInitializeWasPerformed = val;
    }

    bool GetIsInitialized() const
    {
        return mInitializeWasPerformed;
    }

    SchemePointerType &pGetScheme()
    {
        return mpScheme;
    };

    BuilderAndSolverPointerType &pGetEigenSolver()
    {
        return mpEigenSolver;
    };

    BuilderAndSolverPointerType &pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    // SparseMatrixPointerType &pGetMassMatrix()
    // {
    //     return mpMassMatrix;
    // }

    SparseMatrixPointerType &pGetStiffnessMatrixPrevious()
    {
        return mpStiffnessMatrixPrevious;
    }

    SparseMatrixPointerType &pGetStiffnessMatrix()
    {
        return mpStiffnessMatrix;
    }

    SparseVectorPointerType &pGetRHS()
    {
        return mpRHS;
    }

    SparseVectorPointerType &pGetDx()
    {
        return mpDx;
    }

    // SparseMatrixType &GetMassMatrix()
    // {
    //     return *mpMassMatrix;
    // }

    SparseMatrixType &GetStiffnessMatrixPrevious()
    {
        return *mpStiffnessMatrixPrevious;
    }

    SparseMatrixType &GetStiffnessMatrix()
    {
        return *mpStiffnessMatrix;
    }

    SparseVectorType &GetRHS()
    {
        return *mpRHS;
    }

    SparseVectorType &GetDx()
    {
        return *mpDx;
    }

    ConvergenceCriteriaType &GetConvergenceCriteria()
    {
        return mpConvergenceCriteria;
    }

    bool GetSolutionFoundFlag()
    {
        return mSolutionFound;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        this->pBuilderAndSolver()->SetReshapeMatrixFlag(flag);
    }

    bool GetReformDofSetAtEachStepFlag() const
    {
        return this->pBuilderAndSolver()->GetReshapeMatrixFlag();
    }

    /// Set verbosity level of the solving strategy.
    /**
     * - 0 -> mute... no echo at all
     * - 1 -> print time and basic information
     * - 2 -> print linear solver data
     * - 3 -> print debug information
     */
    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        this->pGetEigenSolver()->SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /**
     * Initialization to be performed once before using the strategy.
     */
    void Initialize() override
    {
        KRATOS_TRY

        ModelPart &rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();
        typename ConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            << "Entering Initialize" << std::endl;

        if (mInitializeWasPerformed == false)
        {
            SchemePointerType &pScheme = this->pGetScheme();

            if (pScheme->SchemeIsInitialized() == false)
                pScheme->Initialize(rModelPart);

            if (pScheme->ElementsAreInitialized() == false)
                pScheme->InitializeElements(rModelPart);

            if (pScheme->ConditionsAreInitialized() == false)
                pScheme->InitializeConditions(rModelPart);
        }
        //initialisation of the convergence criteria
        pConvergenceCriteria->Initialize(BaseType::GetModelPart());

        this->SetIsInitialized( true );

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            << "Exiting Initialize" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY

        BuilderAndSolverPointerType &pBuilderAndSolver = this->pGetBuilderAndSolver();
        pBuilderAndSolver->GetLinearSystemSolver()->Clear();
        BuilderAndSolverPointerType &pEigenSolver = this->pGetEigenSolver();
        pEigenSolver->GetLinearSystemSolver()->Clear();

        // if (this->pGetMassMatrix() != nullptr)
        //     this->pGetMassMatrix() = nullptr;

        if (this->pGetStiffnessMatrix() != nullptr)
            this->pGetStiffnessMatrix() = nullptr;

        if (this->pGetStiffnessMatrixPrevious() != nullptr)
            this->pGetStiffnessMatrixPrevious() = nullptr;

        if (this->pGetRHS() != nullptr)
            this->pGetRHS() = nullptr;

        if (this->pGetDx() != nullptr)
            this->pGetDx() = nullptr;

        // Re-setting internal flag to ensure that the dof sets are recalculated
        pBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        pEigenSolver->SetDofSetIsInitializedFlag(false);

        pBuilderAndSolver->Clear();
        pEigenSolver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;
        mSolutionStepIsInitialized = false;


        KRATOS_CATCH("")
    }

    /**
     * Performs all the required operations that should be done (for each step)
     * before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY
        if (!mSolutionStepIsInitialized){
            ModelPart &rModelPart = BaseType::GetModelPart();
            const int rank = rModelPart.GetCommunicator().MyPID();

            KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
                << "Entering InitializeSolutionStep" << std::endl;

            // Solver
            BuilderAndSolverPointerType &pEigenSolver = this->pGetEigenSolver();
            BuilderAndSolverPointerType &pBuilderAndSolver = this->pGetBuilderAndSolver();
            // Scheme
            SchemePointerType &pScheme = this->pGetScheme();
            // Convergence Criteria
            typename ConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
            // Elastic Stiffness
            SparseMatrixPointerType &pStiffnessMatrix = this->pGetStiffnessMatrix();
            SparseMatrixType &rStiffnessMatrix = this->GetStiffnessMatrix();
            // Elastic Stiffness Previous
            SparseMatrixPointerType &pStiffnessMatrixPrevious = this->pGetStiffnessMatrixPrevious();
            SparseMatrixType &rStiffnessMatrixPrevious = this->GetStiffnessMatrixPrevious();
            // Right Hand Side
            SparseVectorPointerType &pRHS = this->pGetRHS();
            SparseVectorType &rRHS = this->GetRHS();
            // Solution Vector
            SparseVectorPointerType &pDx = this->pGetDx();
            SparseVectorType &rDx = this->GetDx();

            // Initialize dummy vectors
            SparseVectorPointerType _pDx = SparseSpaceType::CreateEmptyVectorPointer();
            SparseVectorPointerType _pb = SparseSpaceType::CreateEmptyVectorPointer();
            auto &_rDx = *_pDx;
            auto &_rb = *_pb;

            // Reset solution dofs
            BuiltinTimer system_construction_time;
            if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                pBuilderAndSolver->GetReshapeMatrixFlag() == true)
            {
                // Set up list of dofs
                BuiltinTimer setup_dofs_time;
                
                pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);

                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                // Set global equation ids
                BuiltinTimer setup_system_time;

                pBuilderAndSolver->SetUpSystem(rModelPart);

                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_system_time.ElapsedSeconds() << std::endl;

                // Resize and initialize system matrices
                BuiltinTimer system_matrix_resize_time;

                // Elastic Stiffness matrix; Solution Vector; RHS Vector
                pBuilderAndSolver->ResizeAndInitializeVectors(
                    pScheme, pStiffnessMatrix, pDx, pRHS, rModelPart);
                // Previous Stiffness Matrix 
                pBuilderAndSolver->ResizeAndInitializeVectors(
                    pScheme, pStiffnessMatrixPrevious, _pDx, _pb, rModelPart);

                // Mass Matrix 
                // pBuilderAndSolver->ResizeAndInitializeVectors(
                //     pScheme, pMassMatrix, _pDx, _pb, rModelPart);

                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;

            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;
            
            // Initial operations ... things that are constant over the solution
            // step
            pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),
                                                rStiffnessMatrix, rDx, rRHS);

            // Initial operations ... things that are constant over the solution
            // step
            pScheme->InitializeSolutionStep(BaseType::GetModelPart(), rStiffnessMatrix, rDx, rRHS);

            // Initialisation of the convergence criteria
            // if (pConvergenceCriteria->GetActualizeRHSflag() == true)
            // {
            //     TSparseSpace::SetToZero(rRHS);
            //     pBuilderAndSolver->BuildRHS(pScheme, rModelPart, rRHS);
            // }

            pConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

            // if (pConvergenceCriteria->GetActualizeRHSflag() == true)
            //     TSparseSpace::SetToZero(rRHS);

            mSolutionStepIsInitialized = true;
            KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
                << "Exiting InitializeSolutionStep" << std::endl;
        }
        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        ModelPart &rModelPart = BaseType::GetModelPart();
        ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();

        const int rank = rModelPart.GetCommunicator().MyPID();

        SchemePointerType &pScheme = this->pGetScheme();
        BuilderAndSolverPointerType &pBuilderAndSolver = this->pGetBuilderAndSolver();
        //SparseMatrixType &rMassMatrix = this->GetMassMatrix();
        SparseMatrixType &rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseMatrixType &rStiffnessMatrixPrevious = this->GetStiffnessMatrixPrevious();
        SparseVectorType &rDx = this->GetDx();
        SparseVectorType &rRHS = this->GetRHS();

        typename ConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;
        
        // store inital load condition values
        if( mLoadStepIteration == 0)
        {
            StoreInitialLoadConditions();
        }
        // update load conditions
        // implemented so far are point and surface loads
        // don't change forces in the first load step
        if( mLoadStepIteration > 0)
        {
            UpdateLoadConditions();
        }

        // update loadfactor increment
        double deltaLoadMultiplier = 0.0;
        if( mLoadStepIteration == 1) // inital step
        {
            deltaLoadMultiplier = mInitialStep*(lambda + lambdaPrev);
        }
        else if( mLoadStepIteration % 2 == 1 ) //small step
        {
            deltaLoadMultiplier = mSmallStep*(lambdaPrev );
        }

        BuiltinTimer system_solve_time;
        // initialize nonlinear iteration
        this->pGetScheme()->InitializeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
        pConvergenceCriteria->InitializeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
        is_converged = mpConvergenceCriteria->PreCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

        TSparseSpace::SetToZero(rStiffnessMatrix);
        TSparseSpace::SetToZero(rRHS);
        TSparseSpace::SetToZero(rDx);
        // build system and solve system
        pBuilderAndSolver->BuildAndSolve(pScheme, rModelPart, rStiffnessMatrix,rDx, rRHS);
        // Update internal variables
        pScheme->Update(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
        BaseType::MoveMesh();
        // finalize nonlinear solution step
        pScheme->FinalizeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
        pConvergenceCriteria->FinalizeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
                   
        if (is_converged)
        {
            // if (mpConvergenceCriteria->GetActualizeRHSflag())
            // {
            //     TSparseSpace::SetToZero(rRHS);
            //     this->pGetBuilderAndSolver()->BuildRHS(this->pGetScheme(), rModelPart, rRHS);
            //     //residual_is_updated = true;
            // }

            is_converged = mpConvergenceCriteria->PostCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
        }
        
        // start iteration cycle
        while (is_converged == false &&
            iteration_number++ < mMaxIteration)
        {
            //setting the number of iteration
            rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            // initialize nonlinear iteration step
            pScheme->InitializeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
            pConvergenceCriteria->InitializeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
            is_converged = mpConvergenceCriteria->PreCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

            TSparseSpace::SetToZero(rStiffnessMatrix);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rRHS);
            // build and solve system
            pBuilderAndSolver->BuildAndSolve(pScheme, rModelPart, rStiffnessMatrix,rDx, rRHS);
            
            // update internal variables
            this->pGetScheme()->Update(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
            BaseType::MoveMesh();
            
            // finalize nonlinear iteration step
            this->pGetScheme()->FinalizeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
            pConvergenceCriteria->FinalizeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

            if (is_converged)
            {
                // if (mpConvergenceCriteria->GetActualizeRHSflag())
                // {
                //     TSparseSpace::SetToZero(rRHS);
                //     this->pGetBuilderAndSolver()->BuildRHS(pScheme, rModelPart, rRHS);
                //     //residual_is_updated = true;
                // }

                is_converged = mpConvergenceCriteria->PostCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
            }
        }
        KRATOS_INFO_IF("Nonlinear Loadstep Time: ", BaseType::GetEchoLevel() > 0)
                << system_solve_time.ElapsedSeconds() << std::endl;
        
        if ( !is_converged ) {
                    KRATOS_INFO_IF("Nonlinear Loadstep: ", this->GetEchoLevel() > 0)
            << "Convergence not achieved after ( " << mMaxIteration
            << " ) Iterations !" << std::endl;
        } else {
            KRATOS_INFO_IF("Nonlinear Loadstep: ", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << mMaxIteration << " iterations" << std::endl;
        }

        // vector and matrix are initialized by eigensolver
        DenseVectorType Eigenvalues;
        DenseMatrixType Eigenvectors;
        
        if( mLoadStepIteration % 2 == 0 ) //copy matrices after big step and initial step
        {
            rStiffnessMatrixPrevious = rStiffnessMatrix;
            
            if( mLoadStepIteration > 0) // update lambdaPrev
            {
                lambdaPrev = lambdaPrev + mBigStep*lambda;
            }
        }
        else if( mLoadStepIteration % 2 == 1 ) // evaluate eigenvalue problem after small step
        {
            //#####################################################################################
            //Method 1
            //Poy, Denise Lori-eng; On the Buckling Finite Element Analysis of Beam Structures 
            //#####################################################################################
            // Solve Eigenvalue Problem
            // BuiltinTimer eigen_solve_time;
            // this->pGetEigenSolver()->GetLinearSystemSolver()->Solve(
            //     rStiffnessMatrix,
            //     rStiffnessMatrixPrevious,
            //     Eigenvalues,
            //     Eigenvectors);

            // lambda = 1.0 / ( 1.0 - Eigenvalues(0) ) * deltaLoadMultiplier;
            // for(int i = 0; i < Eigenvalues.size(); i++ )
            // {
            //     Eigenvalues[i] = lambdaPrev + 1.0 / ( 1.0 - Eigenvalues[i] )*deltaLoadMultiplier;
            // }
            // this->AssignVariables(Eigenvalues, Eigenvectors);
            //End Method 1#########################################################################

            //#####################################################################################
            //Method 2: 
            // Jia, X.Mang, H. A.; Assessment of solutions from the consistently linearized eigenproblem 
            // by means of finite difference approximations
            //#####################################################################################
            // Solve Eigenvalue Problem

            // This implementation of the eigenvalue problem relies on the difference of the stiffness matrix 
            // between the current and previous step
            rStiffnessMatrix = rStiffnessMatrixPrevious - rStiffnessMatrix;

            this->pGetEigenSolver()->GetLinearSystemSolver()->Solve(
                rStiffnessMatrixPrevious,
                rStiffnessMatrix,
                Eigenvalues,
                Eigenvectors);

            // Update eigenvalues to loadfactors (Instead of dividing matrix by deltaLoadMultiplier, here eigenvalues are multiplied)
            lambda = Eigenvalues(0)*deltaLoadMultiplier;
            for(int i = 0; i < Eigenvalues.size(); i++ )
            {
                Eigenvalues[i] = lambdaPrev + Eigenvalues[i]*deltaLoadMultiplier;
            }
            
            this->AssignVariables(Eigenvalues, Eigenvectors);
            
            // Reset elastic stiffness matrix for next loadstep
            this->pGetStiffnessMatrix() = nullptr;
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pGetStiffnessMatrix(), pGetDx(), pGetRHS(), rModelPart);

            //End Method 2#########################################################################

            // Convergence criteria for buckling analysis
            if( abs(lambda/lambdaPrev) < mConvergenceRatio )
            {
                mSolutionFound = true;
                KRATOS_INFO_IF("Prebuckling Analysis: ", BaseType::GetEchoLevel() > 0)
                << "Convergence achieved in " << mLoadStepIteration + 1 << " Load Iterations!" << std::endl;
            }   
        }
        
        mLoadStepIteration++;
        
        return true;
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            << "Entering FinalizeSolutionStep" << std::endl;

        typename ConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
        SparseMatrixType &rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseVectorType &rDx = this->GetDx();
        SparseVectorType &rRHS = this->GetRHS();

        pGetBuilderAndSolver()->FinalizeSolutionStep(
            BaseType::GetModelPart(), rStiffnessMatrix, rDx, rRHS);

        pGetScheme()->FinalizeSolutionStep(BaseType::GetModelPart(),
                                           rStiffnessMatrix, rDx, rRHS);

        pConvergenceCriteria->FinalizeSolutionStep(BaseType::GetModelPart(), pGetBuilderAndSolver()->GetDofSet(), rStiffnessMatrix, rDx, rRHS );
        
        //Cleaning memory after the solution
        pGetScheme()->Clean();
        //reset flags for next step
        mSolutionStepIsInitialized = false;

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            << "Exiting FinalizeSolutionStep" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        ModelPart &rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            << "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            << "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

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

    SchemePointerType mpScheme;

    BuilderAndSolverPointerType mpEigenSolver;

    BuilderAndSolverPointerType mpBuilderAndSolver;

    // SparseMatrixPointerType mpMassMatrix;

    SparseMatrixPointerType mpStiffnessMatrix;

    SparseMatrixPointerType mpStiffnessMatrixPrevious;

    SparseVectorPointerType mpRHS;

    SparseVectorPointerType mpDx;

    typename ConvergenceCriteriaType::Pointer mpConvergenceCriteria;
    /// The pointer to the convergence criteria employed

    bool mInitializeWasPerformed = false;

    bool mSolutionStepIsInitialized = false;

    bool mSolutionFound = false;

    unsigned int mLoadStepIteration = 0;

    unsigned int mMaxIteration;

    std::vector<array_1d<double,3>> mpInitialLoads;

    double mInitialStep;
    double mSmallStep;
    double mBigStep;
    double mConvergenceRatio;

    double lambda = 0.0;
    double lambdaPrev = 1.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

/*
Store Initial Loads to ensure that load multiplier always refers to user input
*/
    void StoreInitialLoadConditions()
    {
        ModelPart &rModelPart = BaseType::GetModelPart();
        const auto it_conditions_begin = rModelPart.ConditionsBegin();
        for (int j = 0; j < static_cast<int>(rModelPart.Conditions().size()); j++) 
        {
            auto it_condition = it_conditions_begin + j;
            // Check if condition is point load
            if (it_condition->Has(POINT_LOAD))
            {
                mpInitialLoads.push_back( it_condition->GetValue(POINT_LOAD) );
            }
            // Check if condition is surface load
            if (it_condition->Has(SURFACE_LOAD))
            {
                mpInitialLoads.push_back( it_condition->GetValue(SURFACE_LOAD) );
            }
        } 
    }
/*
Update Load Conditions
*/
    void UpdateLoadConditions()
    {
        ModelPart &rModelPart = BaseType::GetModelPart();
        const auto it_conditions_begin = rModelPart.ConditionsBegin();
        for (int j = 0; j < static_cast<int>(rModelPart.Conditions().size()); j++) 
        {
            auto it_condition = it_conditions_begin + j;
            // Check if condition is point load
            if (it_condition->Has(POINT_LOAD))
            {
                array_1d<double, 3> point_load;
                if( mLoadStepIteration == 1) // Initial step
                {
                    point_load = mpInitialLoads[j] + mInitialStep * mpInitialLoads[j]; 
                }
                else if( mLoadStepIteration % 2 == 0) // Do big step
                {
                    point_load = ( lambdaPrev + mBigStep * lambda ) * mpInitialLoads[j]; 
                }
                else // Do small step
                {
                    point_load = (1 + mSmallStep) * (lambdaPrev) * mpInitialLoads[j];
                }
                // Update Condition
                it_condition->SetValue(POINT_LOAD,point_load);
            }
            // Check if condition is surface load
            if (it_condition->Has(SURFACE_LOAD))
            {
                array_1d<double, 3>& surface_load = it_condition->GetValue(SURFACE_LOAD);
                if( mLoadStepIteration == 1) // initial step
                {
                    surface_load = mpInitialLoads[j] + mInitialStep * mpInitialLoads[j]; 
                }
                else if( mLoadStepIteration % 2 == 0) // Do big step
                {
                    surface_load = ( lambdaPrev + mBigStep * lambda ) * mpInitialLoads[j]; 
                }
                else if( mLoadStepIteration % 2 == 1) // Do small step
                {
                    surface_load = (1 + mSmallStep) * (lambdaPrev) * mpInitialLoads[j];
                }
                // Update Condition
                it_condition->SetValue(SURFACE_LOAD,surface_load);
            }
        }
    }
    // Assign eigenvalues and eigenvectors to kratos variables
    // Copied from eigensolver_strategy
    void AssignVariables(DenseVectorType &rEigenvalues, DenseMatrixType &rEigenvectors )
    {
        ModelPart &rModelPart = BaseType::GetModelPart();
        const std::size_t NumEigenvalues = rEigenvalues.size();

        // store eigenvalues in process info
        rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR] = rEigenvalues;

        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
        {
            ModelPart::NodeType::DofsContainerType &NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();
            
            Matrix &rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
            if (rNodeEigenvectors.size1() != NumEigenvalues || rNodeEigenvectors.size2() != NumNodeDofs)
            {
                rNodeEigenvectors.resize(NumEigenvalues, NumNodeDofs, false);
            }

            // TO BE VERIFIED!! In the current implmentation of Dofs there are nor reordered and only pushec back.
            // // the jth column index of EIGENVECTOR_MATRIX corresponds to the jth nodal dof. therefore,
            // // the dof ordering must not change.
            // if (NodeDofs.IsSorted() == false)
            // {
            //     NodeDofs.Sort();
            // }

            // fill the EIGENVECTOR_MATRIX
            for (std::size_t i = 0; i < NumEigenvalues; i++)
                for (std::size_t j = 0; j < NumNodeDofs; j++)
                {
                    auto itDof = std::begin(NodeDofs) + j;
                    if( !((*itDof)->IsFixed()) )
                    { 
                        rNodeEigenvectors(i, j) = rEigenvectors(i, (*itDof)->EquationId());
                    }
                    else
                    {
                        rNodeEigenvectors(i, j) = 0.0;
                    }
                }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; /* Class PrebucklingStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos */

#endif /* KRATOS_EIGENSOLVER_STRATEGY  defined */
